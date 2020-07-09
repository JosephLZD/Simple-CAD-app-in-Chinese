from OCC.Display.WebGl import x3dom_renderer
from OCC.Core.BRep import BRep_Builder
from OCC.Core.TopoDS import TopoDS_Shape
from OCC.Core.BRepTools import breptools_Read
from OCC.Core.BRepPrimAPI import (BRepPrimAPI_MakeBox, BRepPrimAPI_MakePrism,
                                  BRepPrimAPI_MakeCylinder, BRepPrimAPI_MakeRevol)

# loads brep shape
# cylinder_head = TopoDS_Shape()
# builder = BRep_Builder()
# breptools_Read(cylinder_head, './models/cylinder_head.brep', builder)
myBody = BRepPrimAPI_MakeBox(60, 60, 50).Shape()

# render cylinder head in x3dom
my_renderer = x3dom_renderer.X3DomRenderer()
my_renderer.DisplayShape(myBody)