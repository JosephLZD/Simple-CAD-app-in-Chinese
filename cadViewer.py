#!/usr/bin/env python
#
# Copyright 2020 Doug Blanding (dblanding@gmail.com)
#
# This file is part of cadViewer.
# The latest  version of this file can be found at:
# //https://github.com/dblanding/cadviewer
#
# cadViewer is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# cadViewer is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# if not, write to the Free Software Foundation, Inc.
# 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#


import logging
import math
import sys

from OCC.Core.Visualization import Tesselator
from PyQt5.QtWidgets import QApplication, QMenu, QTreeWidgetItemIterator
from PyQt5.QtGui import QIcon, QPixmap
from OCC.Core.BRep import BRep_Tool
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse
from OCC.Core.BRepBuilderAPI import (BRepBuilderAPI_MakeEdge,
                                     BRepBuilderAPI_MakeFace,
                                     BRepBuilderAPI_MakeWire)
from OCC.Core.BRepFilletAPI import BRepFilletAPI_MakeFillet
from OCC.Core.BRepPrimAPI import (BRepPrimAPI_MakeBox, BRepPrimAPI_MakePrism,
                                  BRepPrimAPI_MakeCylinder, BRepPrimAPI_MakeRevol)
from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_MakeThickSolid
from OCC.Core.gp import gp_Ax1, gp_Ax3, gp_Dir, gp_Pnt, gp_Trsf, gp_Vec, gp_Ax2
from OCC.Core.TColgp import TColgp_Array1OfPnt
from OCC.Core.TopAbs import TopAbs_EDGE
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopoDS import (TopoDS_Vertex, TopoDS_Edge,
                             topods_Edge, topods_Face, topods_Vertex)
from OCC.Core.TopLoc import TopLoc_Location
from OCC.Core.TopTools import TopTools_ListOfShape
from OCCUtils import Topology
import bottle
from mainwindow import MainWindow
import workplane
import scipy
from OCC.Core.GC import GC_MakeCircle
from OCC.Core.TopoDS import topods, topods_Wire, TopoDS_Compound
from OCC.Core.gp import gp_DZ
from OCC.Core.gp import gp_Pnt, gp_OX, gp_Vec, gp_Trsf, gp_DZ, gp_Ax2, gp_Ax3, gp_Pnt2d, gp_Dir2d, gp_Ax2d, gp_DX
from OCC.Core.GC import GC_MakeArcOfCircle, GC_MakeSegment
from OCC.Core.GCE2d import GCE2d_MakeSegment
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeWire, BRepBuilderAPI_MakeFace, \
    BRepBuilderAPI_Transform
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakePrism, BRepPrimAPI_MakeCylinder, BRepPrimAPI_MakeSphere
from OCC.Core.BRepFilletAPI import BRepFilletAPI_MakeFillet
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Fuse
from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_MakeThickSolid, BRepOffsetAPI_ThruSections
from OCC.Core.BRepLib import breplib
from OCC.Core.BRep import BRep_Tool_Surface, BRep_Builder
from OCC.Core.TopoDS import topods, TopoDS_Edge, TopoDS_Compound
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_EDGE, TopAbs_FACE
from OCC.Core.TopTools import TopTools_ListOfShape
from OCC.Core.Geom2d import Geom2d_Ellipse, Geom2d_TrimmedCurve
from OCC.Core.Geom import Geom_CylindricalSurface, Geom_Ellipse, Geom_Curve, Geom_Surface

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG) # set to DEBUG | INFO | ERROR

TOL = 1e-7 # Linear Tolerance
ATOL = TOL # Angular Tolerance
print('TOLERANCE = ', TOL)

#############################################
#
# Workplane creation functions
#
#############################################

def wpBy3Pts(*args):
    """Direction from pt1 to pt2 sets wDir, pt2 is wpOrigin.
    Direction from pt2 to pt3 sets uDir."""
    if win.ptStack:
        # Finish
        p3 = win.ptStack.pop()
        p2 = win.ptStack.pop()
        p1 = win.ptStack.pop()
        wVec = gp_Vec(p1, p2)
        wDir = gp_Dir(wVec)
        origin = p2
        uVec = gp_Vec(p2, p3)
        uDir = gp_Dir(uVec)
        axis3 = gp_Ax3(origin, wDir, uDir)
        wp = workplane.WorkPlane(100, ax3=axis3)
        win.getNewPartUID(wp, typ='w')
        win.clearCallback()
        statusText = "Workplane created."
        win.statusBar().showMessage(statusText)
    else:
        # Initial setup
        win.registerCallback(wpBy3PtsC)
        display.selected_shape = None
        display.SetSelectionModeVertex()
        statusText = "Pick 3 points. Dir from pt1-pt2 sets wDir, pt2 is origin."
        win.statusBar().showMessage(statusText)
        return

def wpBy3PtsC(shapeList, *args):  # callback (collector) for wpBy3Pts
    for shape in shapeList:
        vrtx = topods_Vertex(shape)
        gpPt = BRep_Tool.Pnt(vrtx) # convert vertex to gp_Pnt
        win.ptStack.append(gpPt)
    if len(win.ptStack) == 1:
        statusText = "Now select point 2 (wp origin)."
        win.statusBar().showMessage(statusText)
    elif len(win.ptStack) == 2:
        statusText = "Now select point 3 to set uDir."
        win.statusBar().showMessage(statusText)
    elif len(win.ptStack) == 3:
        wpBy3Pts()

def wpOnFace(*args):
    """ First face defines plane of wp. Second face defines uDir."""
    if not win.faceStack:
        win.registerCallback(wpOnFaceC)
        display.selected_shape = None
        display.SetSelectionModeFace()
        statusText = "Select face for workplane."
        win.statusBar().showMessage(statusText)
        return
    faceU = win.faceStack.pop()
    faceW = win.faceStack.pop()
    wp = workplane.WorkPlane(100, face=faceW, faceU=faceU)
    win.getNewPartUID(wp, typ='w')
    win.clearCallback()
    statusText = "Workplane created."
    win.statusBar().showMessage(statusText)

def wpOnFaceC(shapeList, *args):  # callback (collector) for wpOnFace
    if not shapeList:
        shapeList = []
    for shape in shapeList:
        face = topods_Face(shape)
        win.faceStack.append(face)
    if len(win.faceStack) == 1:
        statusText = "Select face for workplane U direction."
        win.statusBar().showMessage(statusText)
    elif len(win.faceStack) == 2:
        wpOnFace()

def makeWP():   # Default workplane located in X-Y plane at 0,0,0
    wp = workplane.WorkPlane(100)
    win.getNewPartUID(wp, typ='w')
    win.redraw()

#############################################
#
# Create 2d Construction Line functions
#
#############################################

def add_vertex_to_xyPtStack(shapeList):
    """Helper function to convert vertex to gp_Pnt and put on ptStack."""
    wp = win.activeWp
    for shape in shapeList:
        if isinstance(shape, TopoDS_Vertex):  # Guard against wrong type
            vrtx = topods_Vertex(shape)
            pnt = BRep_Tool.Pnt(vrtx) # convert vertex to type <gp_Pnt>
            trsf = wp.Trsf.Inverted()  # New transform. Don't invert wp.Trsf
            pnt.Transform(trsf)
            pt2d = (pnt.X(), pnt.Y())  # 2d point
            win.xyPtStack.append(pt2d)
        else:
            print(f"(Unwanted) shape type: {type(shape)}")

def processLineEdit():
    """pop value from lineEditStack and place on floatStack or ptStack."""

    text = win.lineEditStack.pop()
    if ',' in text:
        try:
            xstr, ystr = text.split(',')
            p = (float(xstr) * win.unitscale, float(ystr) * win.unitscale)
            win.xyPtStack.append(p)
        except:
            print("Problem with processing line edit stack")
    else:
        try:
            win.floatStack.append(float(text))
        except ValueError as e:
            print(f"{e}")

def clineH():   # Horizontal construction line
    if win.xyPtStack:
        wp = win.activeWp
        p = win.xyPtStack.pop()
        win.xyPtStack = []
        wp.hcl(p)
        win.redraw()
    else:
        win.registerCallback(clineHC)
        display.SetSelectionModeVertex()
        win.xyPtStack = []
        win.clearLEStack()
        win.lineEdit.setFocus()
        statusText = "Select point or enter Y-value for horizontal cline."
        win.statusBar().showMessage(statusText)

def clineHC(shapeList, *args):  # callback (collector) for clineH
    add_vertex_to_xyPtStack(shapeList)
    if win.lineEditStack:
        processLineEdit()
    if win.floatStack:
        y = win.floatStack.pop() * win.unitscale
        pnt = (0, y)
        win.xyPtStack.append(pnt)
    if win.xyPtStack:
        clineH()

def clineV():   # Vertical construction line
    if win.xyPtStack:
        wp = win.activeWp
        p = win.xyPtStack.pop()
        win.xyPtStack = []
        wp.vcl(p)
        win.redraw()
    else:
        win.registerCallback(clineVC)
        display.SetSelectionModeVertex()
        win.xyPtStack = []
        win.clearLEStack()
        win.lineEdit.setFocus()
        statusText = "Select point or enter X-value for vertcal cline."
        win.statusBar().showMessage(statusText)

def clineVC(shapeList, *args):  # callback (collector) for clineV
    add_vertex_to_xyPtStack(shapeList)
    if win.lineEditStack:
        processLineEdit()
    if win.floatStack:
        x = win.floatStack.pop() * win.unitscale
        pnt = (x, 0)
        win.xyPtStack.append(pnt)
    if win.xyPtStack:
        clineV()

def clineHV():   # Horizontal + Vertical construction lines
    if win.xyPtStack:
        wp = win.activeWp
        p = win.xyPtStack.pop()
        win.xyPtStack = []
        wp.hvcl(p)
        win.redraw()
    else:
        win.registerCallback(clineHVC)
        display.SetSelectionModeVertex()
        win.xyPtStack = []
        win.clearLEStack()
        win.lineEdit.setFocus()
        statusText = "Select point or enter x,y coords for H+V cline."
        win.statusBar().showMessage(statusText)

def clineHVC(shapeList, *args):  # callback (collector) for clineHV
    add_vertex_to_xyPtStack(shapeList)
    if win.lineEditStack:
        processLineEdit()
    if win.xyPtStack:
        clineHV()

def cline2Pts():
    if len(win.xyPtStack) == 2:
        wp = win.activeWp
        p2 = win.xyPtStack.pop()
        p1 = win.xyPtStack.pop()
        wp.acl(p1, p2)
        win.xyPtStack = []
        win.redraw()
    else:
        win.registerCallback(cline2PtsC)
        display.SetSelectionModeVertex()
        win.xyPtStack = []
        win.clearLEStack()
        win.lineEdit.setFocus()
        statusText = "Select 2 points for Construction Line."
        win.statusBar().showMessage(statusText)

def cline2PtsC(shapeList, *args):  # callback (collector) for cline2Pts
    add_vertex_to_xyPtStack(shapeList)
    if win.lineEditStack:
        processLineEdit()
    if len(win.xyPtStack) == 2:
        cline2Pts()

def clineAng():
    if (win.xyPtStack and win.floatStack):
        wp = win.activeWp
        text = win.floatStack.pop()
        angle = float(text)
        pnt = win.xyPtStack.pop()
        wp.acl(pnt, ang=angle)
        win.xyPtStack = []
        win.redraw()
    else:
        win.registerCallback(clineAngC)
        display.SetSelectionModeVertex()
        win.xyPtStack = []
        win.floatStack = []
        win.lineEditStack = []
        win.lineEdit.setFocus()
        statusText = "在工作平面上选择一个点(或输入x,y坐标)，并输入角度值"
        win.statusBar().showMessage(statusText)

def clineAngC(shapeList, *args):  # callback (collector) for clineAng
    add_vertex_to_xyPtStack(shapeList)
    win.lineEdit.setFocus()
    if win.lineEditStack:
        processLineEdit()
    if (win.xyPtStack and win.floatStack):
        clineAng()

def clineRefAng():
    pass

def clineAngBisec():
    pass

def clineLinBisec():
    if len(win.xyPtStack) == 2:
        wp = win.activeWp
        pnt2 = win.xyPtStack.pop()
        pnt1 = win.xyPtStack.pop()
        wp.lbcl(pnt1, pnt2)
        win.xyPtStack = []
        win.redraw()
    else:
        win.registerCallback(clineLinBisecC)
        display.SetSelectionModeVertex()

def clineLinBisecC(shapeList, *args):  # callback (collector) for clineLinBisec
    add_vertex_to_xyPtStack(shapeList)
    if len(win.xyPtStack) == 2:
        clineLinBisec()

def clinePara():
    pass

def clinePerp():
    pass

def clineTan1():
    pass

def clineTan2():
    pass

def ccirc():
    """Create a c-circle from center & radius or center & Pnt on circle."""
    wp = win.activeWp
    if len(win.xyPtStack) == 2:
        p2 = win.xyPtStack.pop()
        p1 = win.xyPtStack.pop()
        rad = wp.p2p_dist(p1, p2) #计算圆半径
        wp.circle(p1, rad, constr=True)
        win.xyPtStack = []
        win.floatStack = []
        win.redraw()
    elif (win.xyPtStack and win.floatStack):
        pnt = win.xyPtStack.pop()
        rad = win.floatStack.pop() * win.unitscale
        wp.circle(pnt, rad, constr=True)
        win.xyPtStack = []
        win.floatStack = []
        win.redraw()
    else:
        win.registerCallback(ccircC) #刚开始都会调用这个callback函数让用户输入参数
        display.SetSelectionModeVertex()
        win.xyPtStack = []
        win.floatStack = []
        win.lineEditStack = []
        win.lineEdit.setFocus()
        statusText = "选择圆心，输入半径"
        win.statusBar().showMessage(statusText)

def ccircC(shapeList, *args):
    """callback (collector) for ccirc"""
    add_vertex_to_xyPtStack(shapeList)
    win.lineEdit.setFocus()
    if win.lineEditStack:
        processLineEdit()
    if len(win.xyPtStack) == 2:  #两个点，一个圆心，一个圆边上的任一点
        ccirc()
    if (win.xyPtStack and win.floatStack): #一个点（一个圆心），一个数（半径）
        ccirc()

#############################################
#
# Create 2d Edge Profile functions
#
#############################################

def line():
    """Create a profile geometry line between two end points."""
    if len(win.xyPtStack) == 2:
        wp = win.activeWp
        pnt2 = win.xyPtStack.pop()
        pnt1 = win.xyPtStack.pop()
        wp.line(pnt1, pnt2)
        win.xyPtStack = []
        win.redraw()
    else:
        win.registerCallback(lineC)
        display.SetSelectionModeVertex()
        win.xyPtStack = []
        win.lineEdit.setFocus()
        statusText = "选择线段的两个端点."
        win.statusBar().showMessage(statusText)

def lineC(shapeList, *args):
    """callback (collector) for line"""
    add_vertex_to_xyPtStack(shapeList)
    win.lineEdit.setFocus()
    if win.lineEditStack:
        processLineEdit()
    if len(win.xyPtStack) == 2:
        line()

def rect():
    """Create a profile geometry rectangle from two diagonally opposite corners."""
    if len(win.xyPtStack) == 2:
        wp = win.activeWp
        pnt2 = win.xyPtStack.pop()
        pnt1 = win.xyPtStack.pop()
        wp.rect(pnt1, pnt2)
        win.xyPtStack = []
        win.redraw()
    else:
        win.registerCallback(rectC)
        display.SetSelectionModeVertex()
        win.xyPtStack = []
        win.lineEdit.setFocus()
        statusText = "选择矩形的两个对角点."
        win.statusBar().showMessage(statusText)

def rectC(shapeList, *args):
    """callback (collector) for rect"""
    add_vertex_to_xyPtStack(shapeList)
    win.lineEdit.setFocus()
    if win.lineEditStack:
        processLineEdit()
    if len(win.xyPtStack) == 2:
        rect()

def circle():
    """Create a geometry circle from cntr & rad or cntr & pnt on circle."""
    wp = win.activeWp
    if len(win.xyPtStack) == 2:
        p2 = win.xyPtStack.pop()
        p1 = win.xyPtStack.pop()
        rad = wp.p2p_dist(p1, p2)
        wp.circle(p1, rad, constr=False)
        win.xyPtStack = []
        win.floatStack = []
        win.redraw()
    elif (win.xyPtStack and win.floatStack):
        pnt = win.xyPtStack.pop()
        rad = win.floatStack.pop() * win.unitscale
        wp.circle(pnt, rad, constr=False)
        win.xyPtStack = []
        win.floatStack = []
        win.redraw()
    else:
        win.registerCallback(circleC)
        display.SetSelectionModeVertex()
        win.xyPtStack = []
        win.floatStack = []
        win.lineEditStack = []
        win.lineEdit.setFocus()
        statusText = "选择圆心，输入半径(或选择位于圆边的点)"
        win.statusBar().showMessage(statusText)

def circleC(shapeList, *args):
    """callback (collector) for circle"""
    add_vertex_to_xyPtStack(shapeList)
    win.lineEdit.setFocus()
    if win.lineEditStack:
        processLineEdit()
    if len(win.xyPtStack) == 2:
        circle()
    if (win.xyPtStack and win.floatStack):
        circle()

def arcc2p():
    """Create an arc from center pt, start pt and end pt."""
    wp = win.activeWp
    if len(win.xyPtStack) == 3:
        pe = win.xyPtStack.pop()
        ps = win.xyPtStack.pop()
        pc = win.xyPtStack.pop()
        wp.arcc2p(pc, ps, pe)
        win.xyPtStack = []
        win.floatStack = []
        win.redraw()
    else:
        win.registerCallback(arcc2pC)
        display.SetSelectionModeVertex()
        win.xyPtStack = []
        statusText = "选择弧线所在圆的中心，并选择弧线起点和终点"
        win.statusBar().showMessage(statusText)


def arcc2pC(shapeList, *args):
    """callback (collector) for arcc2p"""
    add_vertex_to_xyPtStack(shapeList)
    win.lineEdit.setFocus()
    if win.lineEditStack:
        processLineEdit()
    if len(win.xyPtStack) == 3:
        arcc2p()

def arc3p():
    """Create an arc from start pt, end pt, and 3rd pt on the arc."""
    wp = win.activeWp
    if len(win.xyPtStack) == 3:
        ps = win.xyPtStack.pop()
        pe = win.xyPtStack.pop()
        p3 = win.xyPtStack.pop()
        wp.arc3p(ps, pe, p3)
        win.xyPtStack = []
        win.floatStack = []
        win.redraw()
    else:
        win.registerCallback(arc3pC)
        display.SetSelectionModeVertex()
        win.xyPtStack = []
        statusText = "选择弧线的起点、终点和经过的点"
        win.statusBar().showMessage(statusText)

def arc3pC(shapeList, *args):
    """Callback (collector) for arc3p"""
    add_vertex_to_xyPtStack(shapeList)
    win.lineEdit.setFocus()
    if win.lineEditStack:
        processLineEdit()
    if len(win.xyPtStack) == 3:
        arc3p()

def geom():
    pass

#############################################
#
# 2D Delete functions
#
#############################################

def delCl():
    """Delete selected 2d construction element.

    Todo: Get this working. Able to pre-select lines from the display
    as type <AIS_InteractiveObject> but haven't figured out how to get
    the type <AIS_Line> (or the cline or Geom_Line that was used to make
    it)."""
    wp = win.activeWp
    win.registerCallback(delClC)
    statusText = "选择待删除元素"
    win.statusBar().showMessage(statusText)
    display = win.canva._display.Context
    print(display.NbSelected())  # Use shift-select for multiple lines
    selected_line = display.SelectedInteractive()
    if selected_line:
        print(type(selected_line))  # <AIS_InteractiveObject>
        print(selected_line.GetOwner())  # <Standard_Transient>

def delClC(shapeList, *args):
    """Callback (collector) for delCl"""
    print(shapeList)
    print(args)
    delCl()

def delEl():
    """Delete selected construction element."""
    wp = win.activeWp
    if win.shapeStack:
        while win.shapeStack:
            shape = win.shapeStack.pop()
            if shape in wp.edgeList:
                wp.edgeList.remove(shape)
        win.redraw()
    else:
        win.registerCallback(delElC)
        win.xyPtStack = []
        statusText = "选择待删除元素."
        win.statusBar().showMessage(statusText)

def delElC(shapeList, *args):
    """Callback (collector) for delEl"""
    for shape in shapeList:
        win.shapeStack.append(shape)
    if win.shapeStack:
        delEl()

#############################################
#
# 3D Geometry creation functions
#
#############################################

def sphere_from_vector_and_scalar(X1,Y1,Z1, radius):
    '''
    把创建球体的操作封装起来，且用scipy向量来简化球心位置的赋值
    直接传入球心的位置向量，半径，返回一个球体.
    :param vector: scipy向量，维度为（3，1）
    :param radius: 半径，一个标量
    :return: 球体
    '''
    # X1 = float( vector[0,0] )
    # Y1 = float( vector[1,0] )
    # Z1 = float( vector[2,0] )
    Point = gp_Pnt(X1, Y1, Z1)
    sphere = BRepPrimAPI_MakeSphere(Point, radius)
    return sphere


def luoshuan():
    '''
    螺栓（半径r，长度h）
    '''
    if len(win.floatStack)==2:
        h = win.floatStack.pop() * win.unitscale #注意是栈，先进后出
        r = win.floatStack.pop() * win.unitscale

        name = 'luoshuan'

        coneLocation = gp_Pnt(r,r,0)
        coneAx2 = gp_Ax2(coneLocation, gp_DZ())
        mkCylinder = BRepPrimAPI_MakeCylinder(coneAx2, r, h).Shape()

        #画螺栓顶部（半圆)
        #先画一个box，再画一个直径=边长，圆心=顶边中点的圆，用圆-正方体=半圆
        # neckLocation = gp_Pnt(-r*0.2, -r*0.2, 0)
        # neckAxis = gp_DZ()
        # neckAx2 = gp_Ax2(neckLocation, neckAxis)
        # box = BRepPrimAPI_MakeBox(neckAx2, r*2.4, r*2.4, h).Shape()
        # # 半径
        # radius = r*1.2
        # # 球心位置 创建向量
        # MySphere = sphere_from_vector_and_scalar(r,r,h, radius).Shape()
        # boolean_result = BRepAlgoAPI_Cut(MySphere, box).Shape()

        #画螺栓顶部(六面体 同螺母）
        hole_r = r
        m = h/6
        d = r*2.4
        print("m:",m)
        print("d:",d)
        # 先画个box，再在中间圆柱体，Cut
        corner = gp_Pnt(-d / 2, -d / 2, 0)
        box = BRepPrimAPI_MakeBox(corner, d, d, m).Shape()

        coneLocation = gp_Pnt(0, 0, 0)
        coneAx2 = gp_Ax2(coneLocation, gp_DZ())
        cylinder = BRepPrimAPI_MakeCylinder(coneAx2, hole_r, m).Shape()

        box_cylinder_result = BRepAlgoAPI_Cut(box, cylinder).Shape()

        # 用四个box，旋转之后，cut，生成六面体外形
        import math
        l = math.sqrt(d * d / 2 - d * hole_r + hole_r * hole_r)  # box的边长
        angle = math.acos((d / 2 - hole_r) / l)  # 反余弦函数，计算box旋转角度

        corner1 = gp_Pnt(-d / 2, 0, 0)
        box1 = BRepPrimAPI_MakeBox(corner1, l, l, m).Shape()
        ax1 = gp_Ax1(corner1, gp_Dir(0, 0., 1.))
        aRotTrsf = gp_Trsf()
        aRotTrsf.SetRotation(ax1, angle)
        aTopLoc = TopLoc_Location(aRotTrsf)
        box1.Move(aTopLoc)
        cut_box1_result = BRepAlgoAPI_Cut(box_cylinder_result, box1).Shape()

        corner2 = gp_Pnt(-d / 2, 0, 0)
        box2 = BRepPrimAPI_MakeBox(corner2, l, l, m).Shape()
        aRotTrsf = gp_Trsf()
        aRotTrsf.SetRotation(ax1, -(angle + math.pi / 2))
        aTopLoc = TopLoc_Location(aRotTrsf)
        box2.Move(aTopLoc)
        cut_box2_result = BRepAlgoAPI_Cut(cut_box1_result, box2).Shape()

        corner3 = gp_Pnt(d / 2, 0, 0)
        box3 = BRepPrimAPI_MakeBox(corner3, l, l, m).Shape()
        aRotTrsf = gp_Trsf()
        ax3 = gp_Ax1(corner3, gp_Dir(0, 0., 1.))
        aRotTrsf.SetRotation(ax3, math.pi / 2 - angle)
        aTopLoc = TopLoc_Location(aRotTrsf)
        box3.Move(aTopLoc)
        cut_box3_result = BRepAlgoAPI_Cut(cut_box2_result, box3).Shape()

        corner4 = gp_Pnt(d / 2, 0, 0)
        box4 = BRepPrimAPI_MakeBox(corner4, l, l, m).Shape()
        aRotTrsf = gp_Trsf()
        ax4 = gp_Ax1(corner4, gp_Dir(0, 0., 1.))
        aRotTrsf.SetRotation(ax4, -(math.pi - angle))
        aTopLoc = TopLoc_Location(aRotTrsf)
        box4.Move(aTopLoc)
        cut_box4_result = BRepAlgoAPI_Cut(cut_box3_result, box4).Shape()

        # 将六面体平移到螺栓顶部
        T = gp_Trsf()
        T.SetTranslation(gp_Vec(r, r, h))
        loc = TopLoc_Location(T)
        cut_box4_result.Move(loc)

        luoshuan = BRepAlgoAPI_Fuse(mkCylinder, cut_box4_result).Shape()

        uid = win.getNewPartUID(luoshuan, name=name)
        win.redraw()

    else: #需要提醒用户输入螺栓的半径和高度
        win.registerCallback(luoshuanC)
        win.floatStack = []
        win.lineEditStack = []
        win.lineEdit.setFocus()
        statusText = "输入螺栓的半径:"
        win.statusBar().showMessage(statusText)


def luoshuanC(shapeList, *args):
    #用户需要输入螺栓的半径R和高度l
    if win.lineEditStack: #用户有手动输入半径
        processLineEdit()
    if len(win.floatStack) == 1 :
        statusText = "输入螺栓的高度:"
        win.statusBar().showMessage(statusText)
    elif len(win.floatStack) == 2:
        luoshuan()


def luomu():
    # 螺母
    if len(win.floatStack)==3:
        name = 'luomu'
        #先画个box，再在中间圆柱体，Cut
        hole_r = win.floatStack.pop() * win.unitscale
        m =  win.floatStack.pop() * win.unitscale
        d =  win.floatStack.pop() * win.unitscale
        corner = gp_Pnt(-d/2,-d/2,0)
        box = BRepPrimAPI_MakeBox(corner, d, d, m).Shape()

        coneLocation = gp_Pnt(0, 0, 0)
        coneAx2 = gp_Ax2(coneLocation, gp_DZ())
        cylinder = BRepPrimAPI_MakeCylinder(coneAx2, hole_r, m).Shape()

        box_cylinder_result = BRepAlgoAPI_Cut(box, cylinder).Shape()

        #用四个box，旋转之后，cut，生成六面体外形
        import math
        l = math.sqrt(d*d/2 - d*hole_r + hole_r*hole_r) #box的边长
        angle = math.acos( (d/2-hole_r) / l ) #反余弦函数，计算box旋转角度

        corner1 = gp_Pnt(-d/2, 0, 0)
        box1 = BRepPrimAPI_MakeBox(corner1, l, l, m).Shape()
        ax1 = gp_Ax1(corner1, gp_Dir(0, 0., 1.))
        aRotTrsf = gp_Trsf()
        aRotTrsf.SetRotation(ax1, angle)
        aTopLoc = TopLoc_Location(aRotTrsf)
        box1.Move(aTopLoc)
        cut_box1_result = BRepAlgoAPI_Cut(box_cylinder_result, box1).Shape()

        corner2 = gp_Pnt(-d/2, 0, 0)
        box2 = BRepPrimAPI_MakeBox(corner2, l, l, m).Shape()
        aRotTrsf = gp_Trsf()
        aRotTrsf.SetRotation(ax1, -(angle+math.pi/2))
        aTopLoc = TopLoc_Location(aRotTrsf)
        box2.Move(aTopLoc)
        cut_box2_result = BRepAlgoAPI_Cut(cut_box1_result, box2).Shape()

        corner3 = gp_Pnt(d/2, 0, 0)
        box3 = BRepPrimAPI_MakeBox(corner3, l, l, m).Shape()
        aRotTrsf = gp_Trsf()
        ax3 = gp_Ax1(corner3, gp_Dir(0, 0., 1.))
        aRotTrsf.SetRotation(ax3, math.pi/2-angle)
        aTopLoc = TopLoc_Location(aRotTrsf)
        box3.Move(aTopLoc)
        cut_box3_result = BRepAlgoAPI_Cut(cut_box2_result, box3).Shape()

        corner4 = gp_Pnt(d/2, 0, 0)
        box4 = BRepPrimAPI_MakeBox(corner4, l, l, m).Shape()
        aRotTrsf = gp_Trsf()
        ax4 = gp_Ax1(corner4, gp_Dir(0, 0., 1.))
        aRotTrsf.SetRotation(ax4, -(math.pi-angle))
        aTopLoc = TopLoc_Location(aRotTrsf)
        box4.Move(aTopLoc)
        cut_box4_result = BRepAlgoAPI_Cut(cut_box3_result, box4).Shape()


        uid = win.getNewPartUID(cut_box4_result, name=name)
        win.redraw()

    else:
        win.registerCallback(luomuC)
        win.floatStack = []
        win.lineEditStack = []
        win.lineEdit.setFocus()
        statusText = "输入螺母的长度:"
        win.statusBar().showMessage(statusText)


def luomuC(shapeList, *args):
    #用户需要输入螺母的长度d,厚度m,螺孔半径hole_r
    if win.lineEditStack: #用户有手动输入半径
        processLineEdit()
    if len(win.floatStack) == 1 :
        statusText = "输入螺母的厚度:"
        win.statusBar().showMessage(statusText)
    if len(win.floatStack) == 2 :
        statusText = "输入螺孔半径:"
        win.statusBar().showMessage(statusText)
    elif len(win.floatStack) == 3:
        luomu()




def makeBox():
    name = 'Box'
    myBody = BRepPrimAPI_MakeBox(30, 40, 50).Shape()
    uid = win.getNewPartUID(myBody, name=name)
    win.redraw()
    tess = Tesselator(myBody)
    tess.Compute()
    threejs_string = tess.ExportShapeToThreejsJSONString("json")
    print(threejs_string)

def makeCyl():
    name = 'Cylinder'
    myBody = BRepPrimAPI_MakeCylinder(40, 80).Shape()
    uid = win.getNewPartUID(myBody, name=name)
    win.redraw()

def extrude():
    """Extrude profile on active WP to create a new part."""
    wp = win.activeWp
    if len(win.lineEditStack) == 2:
        name = win.lineEditStack.pop()
        length = float(win.lineEditStack.pop()) * win.unitscale
        wireOK = wp.makeWire()
        if not wireOK:
            print("Unable to make wire.")
            return
        myFaceProfile = BRepBuilderAPI_MakeFace(wp.wire)
        aPrismVec = wp.wVec * length
        myBody = BRepPrimAPI_MakePrism(myFaceProfile.Shape(),
                                       aPrismVec).Shape()
        uid = win.getNewPartUID(myBody, name=name)
        win.redraw()
    else:
        win.registerCallback(extrudeC)
        win.lineEdit.setFocus()
        statusText = "输入拉伸长度，并输入结果组件名称."
        win.statusBar().showMessage(statusText)

def extrudeC(shapeList, *args):
    win.lineEdit.setFocus()
    if len(win.lineEditStack) == 2:
        extrude()

def revolve():
    """Revolve profile on active WP to create a new part."""
    wp = win.activeWp
    if win.lineEditStack and len(win.ptStack) == 2:
        p2 = win.ptStack.pop()
        p1 = win.ptStack.pop()
        name = win.lineEditStack.pop()
        win.clearAllStacks()
        wireOK = wp.makeWire()
        if not wireOK:
            print("Unable to make wire.")
            return
        face = BRepBuilderAPI_MakeFace(wp.wire).Shape()
        revolve_axis = gp_Ax1(p1, gp_Dir(gp_Vec(p1, p2)))
        revolved_shape = BRepPrimAPI_MakeRevol(face, revolve_axis).Shape()
        uid = win.getNewPartUID(revolved_shape, name=name)
        win.statusBar().showMessage('新组件已生成.')
        win.clearCallback()
        win.redraw()
    else:
        win.registerCallback(revolveC)
        win.lineEdit.setFocus()
        statusText = "选择旋转坐标轴上的两个点."
        win.statusBar().showMessage(statusText)

def revolveC(shapeList, *args):
    for shape in shapeList:
        vrtx = topods_Vertex(shape)
        gpPt = BRep_Tool.Pnt(vrtx) # convert vertex to gp_Pnt
        win.ptStack.append(gpPt)
    if len(win.ptStack) == 1:
        statusText = "选择旋转坐标轴上的第二个点."
        win.statusBar().showMessage(statusText)
    elif len(win.ptStack) == 2 and not win.lineEditStack:
        statusText = "输入结果组件名称."
        win.statusBar().showMessage(statusText)
    win.lineEdit.setFocus()
    if win.lineEditStack and len(win.ptStack) == 2:
        revolve()

#############################################
#
# 3D Geometry positioning functons
#
#############################################

def rotateAP():
    ax1 = gp_Ax1(gp_Pnt(0., 0., 0.), gp_Dir(1., 0., 0.))
    aRotTrsf = gp_Trsf()
    angle = math.pi/18 # 10 degrees 每次旋转10度
    aRotTrsf.SetRotation(ax1, angle) #ax1:旋转轴 angle:旋转角度
    aTopLoc = TopLoc_Location(aRotTrsf)
    win.activePart.Move(aTopLoc)
    win.redraw()

def rotateC(*args):
    #输入平移向量
    if win.lineEditStack:
        processLineEdit()
    if len(win.floatStack)==1:
        statusText = "输入y轴偏移量:"
        win.statusBar().showMessage(statusText)
    if len(win.floatStack)==2:
        statusText = "输入z轴偏移量:"
        win.statusBar().showMessage(statusText)
    if len(win.floatStack)==3:
        translation()

def translation():
    if len(win.floatStack)==3:
        z = win.floatStack.pop() * win.unitscale
        y = win.floatStack.pop() * win.unitscale
        x = win.floatStack.pop() * win.unitscale
        T = gp_Trsf()
        T.SetTranslation(gp_Vec(x,y,z))
        loc = TopLoc_Location(T)
        win.activePart.Move(loc)
        win.redraw()
    else:
        win.registerCallback(translationC)
        win.floatStack = []
        win.lineEditStack = []
        win.lineEdit.setFocus()
        statusText = "请输入偏移向量. 首先输入x偏移量:"
        win.statusBar().showMessage(statusText)

def translationC(*args):
    #输入平移向量
    if win.lineEditStack:
        processLineEdit()
    if len(win.floatStack)==1:
        statusText = "输入y偏移量:"
        win.statusBar().showMessage(statusText)
    if len(win.floatStack)==2:
        statusText = "输入z偏移量:"
        win.statusBar().showMessage(statusText)
    if len(win.floatStack)==3:
        translation()


#############################################
#
# 3D Geometry modification functons
#
#############################################

def mill():
    """Mill profile on active WP into active part. 通过二维平面上已有的wire来拉伸，与当前组件相剪切"""
    wp = win.activeWp
    if win.lineEditStack:
        length = float(win.lineEditStack.pop()) * win.unitscale
        wire = wp.wire
        if not wire:
            print("需要先'make wire'.")
            return
        workPart = win.activePart
        wrkPrtUID = win.activePartUID
        punchProfile = BRepBuilderAPI_MakeFace(wire)
        aPrismVec = wp.wVec * length
        tool = BRepPrimAPI_MakePrism(punchProfile.Shape(),
                                       aPrismVec).Shape()
        newPart = BRepAlgoAPI_Cut(workPart, tool).Shape()
        uid = win.getNewPartUID(newPart, ancestor=wrkPrtUID)
        win.statusBar().showMessage('剪切操作完成')
        win.clearCallback()
        win.redraw()
    else:
        win.registerCallback(millC)
        win.lineEdit.setFocus()
        statusText = "输入剪切组件的高度（可为负数）"
        win.statusBar().showMessage(statusText)

def millC(shapeList, *args):
    win.lineEdit.setFocus()
    if win.lineEditStack:
        mill()

def fillet(event=None):
    if (win.lineEditStack and win.edgeStack):
        text = win.lineEditStack.pop()
        filletR = float(text) * win.unitscale
        edges = []
        for edge in win.edgeStack:
            edges.append(edge)
        win.edgeStack = []
        workPart = win.activePart
        wrkPrtUID = win.activePartUID
        mkFillet = BRepFilletAPI_MakeFillet(workPart)
        for edge in edges:
            mkFillet.Add(filletR, edge)
        newPart = mkFillet.Shape()
        win.getNewPartUID(newPart, ancestor=wrkPrtUID)
        win.statusBar().showMessage('圆弧化操作完成')
        win.clearCallback()
    else:
        win.registerCallback(filletC)
        display.SetSelectionModeEdge()
        statusText = "选择需要圆弧化的边，并输入圆角半径"
        win.statusBar().showMessage(statusText)

def filletC(shapeList, *args):  # callback (collector) for fillet
    win.lineEdit.setFocus()
    for shape in shapeList:
        edge = topods_Edge(shape)
        win.edgeStack.append(edge)
    if (win.edgeStack and win.lineEditStack):
        fillet()

def fuse():
    """Fuse two solid shapes together."""
    if win.shapeStack:
        shape = win.shapeStack.pop()
        workpart = win.activePart
        wrkPrtUID = win.activePartUID
        newPart = BRepAlgoAPI_Fuse(workpart, shape).Shape()
        win.getNewPartUID(newPart, ancestor=wrkPrtUID)
        win.statusBar().showMessage('合并操作完成')
        win.clearCallback()
    else:
        win.registerCallback(fuseC)
        statusText = "选择需要与当前组件结合的模型"
        win.statusBar().showMessage(statusText)

def fuseC(shapeList, *args):  # callback (collector) for fuse
    for shape in shapeList:
        win.shapeStack.append(shape)
    if win.shapeStack:
        fuse()

def shell(event=None):
    if (win.lineEditStack and win.faceStack):
        text = win.lineEditStack.pop()
        faces = TopTools_ListOfShape()
        for face in win.faceStack:
            faces.Append(face)
        win.faceStack = []
        workPart = win.activePart
        wrkPrtUID = win.activePartUID
        shellT = float(text) * win.unitscale
        newPart = BRepOffsetAPI_MakeThickSolid(workPart, faces, -shellT, 1.e-3).Shape()
        win.getNewPartUID(newPart, ancestor=wrkPrtUID)
        win.statusBar().showMessage('空心操作完成')
        win.clearCallback()
    else:
        win.registerCallback(shellC)
        display.SetSelectionModeFace()
        statusText = "选择需要去除的面(可多个面)，并输入壳的厚度"
        win.statusBar().showMessage(statusText)

def shellC(shapeList, *args):  # callback (collector) for shell
    win.lineEdit.setFocus()
    for shape in shapeList:
        face = topods_Face(shape)
        win.faceStack.append(face)
    if (win.faceStack and win.lineEditStack):
        shell()





#####################
#                   #
#   Bottle Demo:    #
#                   #
#####################

# Make Bottle step by step
def makePoints(event=None):
    V1, V2, V3, V4, V5, V6 = bottle.makePoints()  # gp_Pnt
    display.DisplayShape(V1.Vertex())
    display.DisplayShape(V2.Vertex())
    display.DisplayShape(V3.Vertex())
    display.DisplayShape(V4.Vertex())
    display.DisplayShape(V5.Vertex())
    display.DisplayShape(V6.Vertex())
    display.FitAll()
    display.EraseAll()
    display.DisplayShape(V1.Vertex())
    display.DisplayShape(V2.Vertex())
    display.DisplayShape(V3.Vertex())
    display.DisplayShape(V4.Vertex())
    display.DisplayShape(V5.Vertex())
    display.Repaint()
    win.statusBar().showMessage('Make Points complete')

def makeLines(event=None):
    e1, e2, e3 = bottle.makeLines()  # TopoDS_Edge
    display.DisplayColoredShape(e1, 'RED')
    display.DisplayColoredShape(e2, 'RED')
    display.DisplayColoredShape(e3, 'RED')
    display.Repaint()
    win.statusBar().showMessage('Make lines complete')

def makeHalfWire(event=None):
    aWire = bottle.makeHalfWire()  # TopoDS_Wire
    display.EraseAll()
    display.DisplayColoredShape(aWire, 'BLUE')
    display.Repaint()
    win.statusBar().showMessage('Make Half Wire complete')

def makeWholeWire(event=None):
    myWireProfile = bottle.makeWholeWire()  # TopoDS_Wire
    display.DisplayColoredShape(myWireProfile, 'BLUE')
    display.Repaint()
    win.statusBar().showMessage('Make whole wire complete')

def makeFace(event=None):
    bottomFace = bottle.makeFace()  # TopoDS_Face
    display.DisplayShape(bottomFace, color='YELLOW', transparency=0.6)
    display.Repaint()
    win.statusBar().showMessage('Make face complete')

def makeBody(event=None):
    myBody = bottle.makeBody()  # TopoDS_Shape
    partName = 'body'
    win.getNewPartUID(myBody, name=partName)
    win.statusBar().showMessage('Bottle body complete')
    win.redraw()

def makeFillets(event=None):
    workPart = win.activePart
    wrkPrtUID = win.activePartUID
    newPrtName = 'bodyWithFillets'
    mkFillet = BRepFilletAPI_MakeFillet(workPart)
    aEdgeExplorer = TopExp_Explorer(workPart, TopAbs_EDGE)
    while aEdgeExplorer.More():
        aEdge = topods_Edge(aEdgeExplorer.Current())
        mkFillet.Add(bottle.thickness / 12., aEdge)
        aEdgeExplorer.Next()
    myBody = mkFillet.Shape()
    win.getNewPartUID(myBody, name=newPrtName, ancestor=wrkPrtUID)
    win.statusBar().showMessage('Bottle with fillets complete')
    win.redraw()

def addNeck(event=None):
    newPrtName = 'bodyWithNeck'
    workPart = win.activePart
    wrkPrtUID = win.activePartUID
    myNeck = bottle.addNeck()  # TopoDS_Shape
    myBody = BRepAlgoAPI_Fuse(workPart, myNeck).Shape()
    win.getNewPartUID(myBody, name=newPrtName, ancestor=wrkPrtUID)
    win.statusBar().showMessage('Add neck complete')
    win.redraw()

#############################################
#
#  Info & Utility functions
#
#############################################

def topoDumpAP():
    Topology.dumpTopology(win.activePart)

def printCurrUID():
    print(win._currentUID)

def printActiveAsyInfo():
    uid = win.activeAsyUID
    treeNode = win.activeAsy
    print(f"Active Assembly Name: {treeNode} \t UID: {uid}")

def printActiveWpInfo():
    print(f"Name: {win.activeWp}")
    print(f"UID: {win.activeWpUID}")

def printActivePartInfo():
    uid = win.activePartUID
    name = win._nameDict.get(uid)
    print(f"Active Part Name: {name} \t UID: {uid}")

def printPartsInActiveAssy():
    asyPrtTree = []
    leafNodes = win.treeModel.leaves(win.activeAsyUID)
    for node in leafNodes:
        pid = node.identifier
        if pid in win._partDict:
            asyPrtTree.append(pid)
    print(asyPrtTree)

def printActPart():
    uid = win.activePartUID
    if uid:
        name = win._nameDict[uid]
        print("Active Part: %s [uid=%i]" % (name, int(uid)))
    else:
        print(None)

def printTreeView():
    """Print 'uid'; 'name'; 'parent' for all items in treeView."""
    iterator = QTreeWidgetItemIterator(win.treeView)
    while iterator.value():
        item = iterator.value()
        name = item.text(0)
        strUID = item.text(1)
        uid = int(strUID)
        pname = None
        parent = item.parent()
        if parent:
            puid = parent.text(1)
            pname = parent.text(0)
        print(f"UID: {uid}; Name: {name}; Parent: {pname}")
        iterator += 1

def printDrawList():
    print("Draw List:", win.drawList)

def printInSync():
    print(win.inSync())

def setUnits_in():
    win.setUnits('in')

def setUnits_mm():
    win.setUnits('mm')



if __name__ == '__main__':
    app = QApplication(sys.argv)
    win = MainWindow()
    win.add_menu('文件')
    win.add_function_to_menu('文件', "加载 STEP", win.loadStep)
    win.add_function_to_menu('文件', "保存 STEP", win.saveStep)
    win.add_function_to_menu('文件', "保存 STEP (仅限于当前组件)", win.saveStepActPrt)
    win.add_menu('工作平面')
    win.add_function_to_menu('工作平面', "依据平面来创建", wpOnFace)
    win.add_function_to_menu('工作平面', "依据三个点来创建", wpBy3Pts)
    win.add_function_to_menu('工作平面', "(默认)底面", makeWP)
    win.add_menu('三维操作')
    win.add_function_to_menu('三维操作', "长方体", makeBox)
    win.add_function_to_menu('三维操作', "圆柱体", makeCyl)
    win.add_function_to_menu('三维操作', "螺栓", luoshuan)
    win.add_function_to_menu('三维操作', "螺母", luomu)
    win.add_function_to_menu('三维操作', "拉伸生成", extrude)
    win.add_function_to_menu('三维操作', "旋转生成", revolve)
    win.add_menu('修改组件')
    win.add_function_to_menu('修改组件', "旋转当前组件", rotateAP)
    win.add_function_to_menu('修改组件', "平移当前组件", translation)
    win.add_function_to_menu('修改组件', "剪切", mill)
    win.add_function_to_menu('修改组件', "圆弧化", fillet)
    win.add_function_to_menu('修改组件', "空心", shell)
    win.add_function_to_menu('修改组件', "合并", fuse)
    # excised dynamic3Dmodification functions
    #win.add_function_to_menu('Modify Active Part', "Lift Face", lift)
    #win.add_function_to_menu('Modify Active Part', "Offset Face", offsetFace)
    #win.add_function_to_menu('Modify Active Part', "Align Face", alignFace)
    #win.add_function_to_menu('Modify Active Part', "Tweak Face", tweakFace)
    #win.add_function_to_menu('Modify Active Part', "Fuse", fuse)
    #win.add_function_to_menu('Modify Active Part', "Remove Face", remFace)
    # win.add_menu('Bottle')
    # win.add_function_to_menu('Bottle', "Step 1: points", makePoints)
    # win.add_function_to_menu('Bottle', "Step 2: lines", makeLines)
    # win.add_function_to_menu('Bottle', "Step 3: half wire", makeHalfWire)
    # win.add_function_to_menu('Bottle', "Step 4: whole wire", makeWholeWire)
    # win.add_function_to_menu('Bottle', "Step 5: face", makeFace)
    # win.add_function_to_menu('Bottle', "Step 6: body", makeBody)
    # win.add_function_to_menu('Bottle', "Step 7: fillets", makeFillets)
    # win.add_function_to_menu('Bottle', "Step 8: neck", addNeck)
    win.add_menu('其他')
    win.add_function_to_menu('其他', "当前组件坐标数据", topoDumpAP)
    win.add_function_to_menu('其他', "当前组件的UID", printCurrUID)
    win.add_function_to_menu('其他', "当前模型结构", printTreeView)
    win.add_function_to_menu('其他', "当前工作平面UID", printActiveWpInfo)
    win.add_function_to_menu('其他', "当前模型信息", printActiveAsyInfo)
    win.add_function_to_menu('其他', "当前组件信息", printActivePartInfo)
    win.add_function_to_menu('其他', "清空输入框", win.clearLEStack)
    # win.add_function_to_menu('其他', "Calculator", win.launchCalc)
    win.add_function_to_menu('其他', "设置单位为in", setUnits_in)
    win.add_function_to_menu('其他', "设置单位为mm", setUnits_mm)

    drawSubMenu = QMenu('绘制')
    win.popMenu.addMenu(drawSubMenu)
    drawSubMenu.addAction('放缩至适应屏幕', win.fitAll)
    drawSubMenu.addAction('重画', win.redraw)
    drawSubMenu.addAction('隐藏全部', win.eraseAll)
    drawSubMenu.addAction('绘制全部', win.drawAll)
    drawSubMenu.addAction('只绘制当前组件', win.drawOnlyActivePart)

    win.treeView.popMenu.addAction('设置为当前', win.setClickedActive)
    win.treeView.popMenu.addAction('设置为透明', win.setTransparent)
    win.treeView.popMenu.addAction('设置为不透明', win.setOpaque)
    win.treeView.popMenu.addAction('重命名', win.editName)

    win.show()
    win.canva.InitDriver()
    display = win.canva._display

    selectSubMenu = QMenu('选择模式')
    win.popMenu.addMenu(selectSubMenu)
    selectSubMenu.addAction('点', display.SetSelectionModeVertex)
    selectSubMenu.addAction('边', display.SetSelectionModeEdge)
    selectSubMenu.addAction('面', display.SetSelectionModeFace)
    selectSubMenu.addAction('模型', display.SetSelectionModeShape)
    selectSubMenu.addAction('不限', display.SetSelectionModeNeutral)
    win.popMenu.addAction('清空输入缓存', win.clearCallback)
    # Construction Line Toolbar buttons
    win.wcToolBar.addAction(QIcon(QPixmap('icons/hcl.gif')), '水平线', clineH)
    win.wcToolBar.addAction(QIcon(QPixmap('icons/vcl.gif')), '垂直线', clineV)
    win.wcToolBar.addAction(QIcon(QPixmap('icons/hvcl.gif')), '水平 + 垂直', clineHV)
    win.wcToolBar.addAction(QIcon(QPixmap('icons/tpcl.gif')), '两点连线', cline2Pts)
    win.wcToolBar.addAction(QIcon(QPixmap('icons/acl.gif')), '倾斜线', clineAng)
    #win.wcToolBar.addAction(QIcon(QPixmap('icons/refangcl.gif')), 'Ref-Ang', clineRefAng)
    #win.wcToolBar.addAction(QIcon(QPixmap('icons/abcl.gif')), 'Angular Bisector', clineAngBisec)
    win.wcToolBar.addAction(QIcon(QPixmap('icons/lbcl.gif')), '平分线', clineLinBisec)
    #win.wcToolBar.addAction(QIcon(QPixmap('icons/parcl.gif')), 'Parallel', clinePara)
    #win.wcToolBar.addAction(QIcon(QPixmap('icons/perpcl.gif')), 'Perpendicular', clinePerp)
    #win.wcToolBar.addAction(QIcon(QPixmap('icons/cltan1.gif')), 'Tangent to circle', clineTan1)
    #win.wcToolBar.addAction(QIcon(QPixmap('icons/cltan2.gif')), 'Tangent 2 circles', clineTan2)
    win.wcToolBar.addAction(QIcon(QPixmap('icons/ccirc.gif')), '圆', ccirc)
    #win.wcToolBar.addAction(QIcon(QPixmap('icons/cc3p.gif')), 'Circle by 3Pts', ccirc)
    #win.wcToolBar.addAction(QIcon(QPixmap('icons/cccirc.gif')), 'Concentric Circle', ccirc)
    #win.wcToolBar.addAction(QIcon(QPixmap('icons/cctan2.gif')), 'Circ Tangent x2', ccirc)
    #win.wcToolBar.addAction(QIcon(QPixmap('icons/cctan3.gif')), 'Circ Tangent x3', ccirc)
    win.wcToolBar.addSeparator()
    #win.wcToolBar.addAction(QIcon(QPixmap('icons/del_cel.gif')), 'Delete Constr', delCl)
    # Profile Line Toolbar buttons
    win.wgToolBar.addAction(QIcon(QPixmap('icons/line.gif')), '线段', line)
    win.wgToolBar.addAction(QIcon(QPixmap('icons/rect.gif')), '矩形', rect)
    #win.wgToolBar.addAction(QIcon(QPixmap('icons/poly.gif')), 'Polygon', geom)
    #win.wgToolBar.addAction(QIcon(QPixmap('icons/slot.gif')), 'Slot', geom)
    win.wgToolBar.addAction(QIcon(QPixmap('icons/circ.gif')), '圆', circle)
    win.wgToolBar.addAction(QIcon(QPixmap('icons/arcc2p.gif')), '两点弧线', arcc2p)
    win.wgToolBar.addAction(QIcon(QPixmap('icons/arc3p.gif')), '三点弧线', arc3p)
    win.wgToolBar.addSeparator()
    #win.wgToolBar.addAction(QIcon(QPixmap('icons/translate.gif')), 'Translate Profile', geom)
    #win.wgToolBar.addAction(QIcon(QPixmap('icons/rotate.gif')), 'Rotate Profile', geom)
    win.wgToolBar.addAction(QIcon(QPixmap('icons/del_el.gif')), '删除边', delEl)

    win.raise_() # bring the app to the top
    app.exec_()
