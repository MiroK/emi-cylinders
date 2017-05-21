include("../emi.jl")
using EMI
using EMI.Draw, EMI.Gmsh

loop = Loop([Line(Point(0, 0), Point(1, 0)), 
             CircleArc(Point(1, 0), Point(0, 0), Point(0, 1)),
             Line(Point(0, 1), Point(0, 0))])

canvas = Canvas()
canvas = canvas + loop

set_bbox!(canvas, 0.2, 0.3)
println(gmsh_script(canvas, 0.2))
