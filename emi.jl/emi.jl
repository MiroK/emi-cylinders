module EMI

module Draw
include("emi_draw.jl")

export Point
export Curve, Line, CircleArc
export Shape, Circle, Ellipse, Rectangle, Square, ClosedPolygon, CompositeLoop
export BoundingBox
export Canvas, set_bbox!

end  # Draw

module Gmsh
using EMI.Draw
include("emi_gmsh.jl")

export gmsh_script

end  # Gmsh 

end # EMI
