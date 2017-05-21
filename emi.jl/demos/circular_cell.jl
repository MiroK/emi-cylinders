include("../emi.jl")
using EMI
using EMI.Draw, EMI.Gmsh


# Something like this with circular
#   ^
#  / \
# [   ]
#  \ /
#   U


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
canvas = canvas + loop

set_bbox!(canvas, 0.2, 0.3)
gmsh_script(canvas, 0.2)
