include("../emi.jl")
using EMI
using EMI.Draw, EMI.Gmsh

# from .. run as emi-gmsh demos/FILE.JL demos/FILE
# Something like this with ellipse
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

canvas = Canvas()
canvas = canvas + loop

set_bbox!(canvas, 0.2, 0.3)
gmsh_script(canvas, 0.2)
