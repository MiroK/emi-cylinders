include("emi.jl")
using EMI
using EMI.Draw, EMI.Gmsh

# SOMETHING like this
#   ^
#  / \
# [   ]
#  \ /
#   U

const M, m, dx, dy = 1, 0.8, 0.2, 0.3
@assert M > m
loop = Loop([Line(Point(0, dx+dy+m), Point(0, dy+m)),
             Line(Point(0, dy+m), Point(dy, dy+m)),
             EllipseArc(Point(dy, dy+m), Point(dy+M, dy+m), Point(dy+M, dy)),
             Line(Point(dy+M, dy), Point(dy+M, 0)),
             Line(Point(dy+M, 0), Point(dx+dy+M, 0)),
             Line(Point(dx+dy+M, 0), Point(dx+dy+M, dy)),
             EllipseArc(Point(dx+dy+2*M, dy+m), Point(dx+dy+M, dy+m), Point(dx+dy+M, dy)),
             Line(Point(dx+dy+2*M, dy+m), Point(dx+2*dy+2*M, dy+m)),
             Line(Point(dx+2*dy+2*M, dy+m), Point(dx+2*dy+2*M, dx+dy+m)),
             Line(Point(dx+2*dy+2*M, dx+dy+m), Point(dx+dy+2*M, dx+dy+m)),
             EllipseArc(Point(dx+dy+2*M, dx+dy+m), Point(dx+dy+M, dx+dy+m), Point(dx+dy+M, dx+dy+2*m)),
             Line(Point(dx+dy+M, dx+dy+2*m), Point(dx+dy+M, dx+2*dy+2*m)),
             Line(Point(dx+dy+M, dx+2*dy+2*m), Point(dy+M, dx+2*dy+2*m)),
             Line(Point(dy+M, dx+2*dy+2*m), Point(dy+M, dx+dy+2*m)),
             EllipseArc(Point(dy, dx+dy+m), Point(dy+M, dx+dy+m), Point(dy+M, dx+dy+2*m)),
             Line(Point(dy, dx+dy+m), Point(0, dx+dy+m))])

const R, dx, dy = 1, 0.2, 0.3
loop = Loop([Line(Point(0, R+dy+dx), Point(0, R+dy)),
             Line(Point(0, R+dy), Point(dy, R+dy)),
             CircleArc(Point(dy, R+dy), Point(R+dy, R+dy), Point(R+dy, dy)),
             Line(Point(R+dy, dy), Point(R+dy, 0)),
             Line(Point(R+dy, 0), Point(R+dy+dx, 0)),
             Line(Point(R+dy+dx, 0), Point(R+dy+dx, dy)),
             CircleArc(Point(R+dy+dx, dy), Point(R+dy+dx, R+dy), Point(2R+dy+dx, R+dy)),
             Line(Point(2R+dy+dx, R+dy), Point(2R+2*dy+dx, R+dy)),
             Line(Point(2R+2*dy+dx, R+dy), Point(2R+2*dy+dx, R+dy+dx)),
             Line(Point(2R+2*dy+dx, R+dy+dx), Point(2R+dy+dx, R+dy+dx)),
             CircleArc(Point(2R+dy+dx, R+dy+dx), Point(R+dy+dx, R+dy+dx), Point(R+dy+dx, 2R+dy+dx)),
             Line(Point(R+dy+dx, 2R+dy+dx), Point(R+dy+dx, 2R+2*dy+dx)),
             Line(Point(R+dy+dx, 2R+2*dy+dx), Point(R+dy, 2R+2*dy+dx)),
             Line(Point(R+dy, 2R+2*dy+dx), Point(R+dy, 2R+dy+dx)),
             CircleArc(Point(R+dy, 2R+dy+dx), Point(R+dy, R+dy+dx), Point(dy, R+dy+dx)),
             Line(Point(dy, R+dy+dx), Point(0, R+dy+dx))])

canvas = Canvas()
#canvas = canvas + Circle(Point(0, 0), 1)
#canvas = canvas + Ellipse(Point(4, 0), 1, 0.5)
canvas = canvas + loop

set_bbox!(canvas, 0.2, 0.3)
println(gmsh_script(canvas))
