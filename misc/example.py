import sys
from math import pi

from OCC.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCC.AIS import *
from OCC.Quantity import *
from OCC.Display.SimpleGui import init_display
from OCC.TopoDS import *
from OCC.gp import *
from OCC.TopLoc import *
from OCC.Geom import *
from OCC.BRep import BRep_Tool_Surface
from OCC.GCE2d import *
from OCC.BRepBuilderAPI import *

display, start_display, add_menu, add_function_to_menu = init_display()


def geom_plane_from_face(aFace):
    """
    Returns the geometric plane entity from a planar surface
    """
    return Handle_Geom_Plane.DownCast(OCC.BRep.BRep_Tool_Surface(aFace)).GetObject()

def redraw(event=None):
    # display with crisp edges and transpaarency
    context = display.Context
    context.RemoveAll()
    context.SetAutoActivateSelection(False)
    aisShape = AIS_Shape(Box)
    h_aisShape = aisShape.GetHandle()
    context.Display(h_aisShape)
    context.SetTransparency(h_aisShape, .6)
    context.HilightWithColor(h_aisShape, OCC.Quantity.Quantity_NOC_BLACK)
    display.FitAll()
    
def makeBox(event=None):
    global Box
    # Make a box
    Box = BRepPrimAPI_MakeBox(60, 60, 50).Shape()
    redraw()
    
def rotateBox():
    aisShape = AIS_Shape(Box)
    ax1 = gp_Ax1(gp_Pnt(0., 0., 0.), gp_Dir(1., 0., 0.))
    aRotTrsf = gp_Trsf()
    angle = pi/6
    aRotTrsf.SetRotation(ax1, angle)
    aTopLoc = TopLoc_Location(aRotTrsf)
    Box.Move(aTopLoc)
    redraw()

def enableFaceSelect(event=None):
    display.selected_shape = None
    display.SetSelectionModeFace()

def makeSqProfile(size, surface):
    # points and segments need to be in CW sequence to get W pointing along Z
    aPnt1 = gp_Pnt2d(-size, size)
    aPnt2 = gp_Pnt2d(size, size)
    aPnt3 = gp_Pnt2d(size, -size)
    aPnt4 = gp_Pnt2d(-size, -size)
    aSegment1 = GCE2d_MakeSegment(aPnt1, aPnt2)
    aSegment2 = GCE2d_MakeSegment(aPnt2, aPnt3)
    aSegment3 = GCE2d_MakeSegment(aPnt3, aPnt4)
    aSegment4 = GCE2d_MakeSegment(aPnt4, aPnt1)
    print 'Next is where something crashes'
    aEdge1 = BRepBuilderAPI_MakeEdge(aSegment1.Value(),
                                     Handle_Geom_Surface(surface))
    aEdge2 = BRepBuilderAPI_MakeEdge(aSegment2.Value(),
                                     Handle_Geom_Surface(surface))
    aEdge3 = BRepBuilderAPI_MakeEdge(aSegment3.Value(),
                                     Handle_Geom_Surface(surface))
    aEdge4 = BRepBuilderAPI_MakeEdge(aSegment4.Value(),
                                     Handle_Geom_Surface(surface))
    print "Doesn't get here (with rotated box)"
    aWire = BRepBuilderAPI_MakeWire(aEdge1.Edge(),
                                    aEdge2.Edge(),
                                    aEdge3.Edge(),
                                    aEdge4.Edge())

    myWireProfile = aWire.Wire()
    return myWireProfile # TopoDS_Wire

def wireProfileOnFace(event=None):
    aShape = display.GetSelectedShape()
    shapes = display.GetSelectedShapes()
    face = None
    if aShape:
        face = topods_Face(aShape)
        print "A shape found:"
    elif shapes:
        aShape = shapes[0]
        face = topods_Face(aShape)
        print len(shapes), "Shapes found"
    if face:
        surface = geom_plane_from_face(face)
        wireProfile = makeSqProfile(50, surface)
        display.DisplayShape(wireProfile)
    else:
        print 'no face'

def exit(event=None):
    sys.exit()

if __name__ == '__main__':
    add_menu('operations')
    add_function_to_menu('operations', makeBox)
    add_function_to_menu('operations', rotateBox)
    add_function_to_menu('operations', enableFaceSelect)
    add_function_to_menu('operations', wireProfileOnFace)
    add_function_to_menu('operations', exit)
    start_display()
