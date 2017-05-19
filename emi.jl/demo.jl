include("emi.jl")
using EMI
using EMI.Draw

# SOMETHING like this
#   ^
#  / \
# [   ]
#  \ /
#   U

path = EMI.Draw.make_path([
                           Line(Point(0, 4), Point(0, 3)),
                           Line(Point(0, 3), Point(1, 3)),
                           CircleArc(Point(1, 3), Point(3, 3), Point(3, 1)),
                           Line(Point(3, 1), Point(3, 0))
                          ])

println(path)
