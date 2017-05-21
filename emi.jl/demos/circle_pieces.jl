include("../emi.jl")
using EMI
using EMI.Draw, EMI.Gmsh

loop = Loop([CircleArc(Point(1, 0), Point(0, 0), Point(0, 1)),
             CircleArc(Point(0, 1), Point(0, 0), Point(-1, 0)),
             CircleArc(Point(-1, 0), Point(0, 0), Point(0, -1)),
             CircleArc(Point(0, -1), Point(0, 0), Point(1, 0))])

canvas = Canvas()
canvas = canvas + loop

set_bbox!(canvas, 0.2, 0.3)
gmsh_script(canvas, 0.2)
