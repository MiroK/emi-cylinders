include("../emi.jl")
using EMI
using EMI.Draw, EMI.Gmsh

const R = 1
const C = R/sqrt(2)

loop = Loop([CircleArc(Point(-R, 0), Point(0, 0), Point(0, R)),
             CircleArc(Point(0, R), Point(0, 0), Point(C, C)),
             Line(Point(C, C), Point(0., 0.)),
             Line(Point(0., 0.), Point(C, -C)),
             CircleArc(Point(C, -C), Point(0., 0.), Point(0, -R)),
             CircleArc(Point(0, -R), Point(0, 0), Point(-R, 0))])

canvas = Canvas()
canvas = canvas + loop

set_bbox!(canvas, 0.2, 0.2)
gmsh_script(canvas, 0.2)
