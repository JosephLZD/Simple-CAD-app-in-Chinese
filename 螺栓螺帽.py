import logging
import math
import sys
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

def luoshuan():

    name = 'luoshuan'

    point = scipy.array([0.0, 0.0, 0.0])
    point = scipy.reshape(point, (3, 1))
    rad = 50
    circle=GC_MakeCircle(gp_Ax2(gp_Pnt(0,0,0), gp_Dir(0,0,0)), 6).Value()

    # 创建用于体拉伸 sweep 的面
    myFaceProfile = BRepBuilderAPI_MakeFace(circle)

    # 我们把面沿Z轴拉伸到指定高度
    height = 35
    aPrismVec = gp_Vec(0, 0, height)
    myBody = BRepPrimAPI_MakePrism(myFaceProfile.Face(), aPrismVec)

    uid = win.getNewPartUID(myBody, name=name)
    win.redraw()



if __name__ == '__main__':
    app = QApplication(sys.argv)
    win = MainWindow()
    win.add_menu('draw')
    win.add_function_to_menu('draw', "螺栓", luoshuan)
    # win.add_function_to_menu('draw', "螺帽", luomao)
    win.raise_() # bring the app to the top
    app.exec_()