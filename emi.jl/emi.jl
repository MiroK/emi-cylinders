module EMI

module Draw
include("emi_draw.jl")

export Point, Shape, Circle, Ellipse, Rectangle, Square, ClosedPolygon, BoundingBox
export Canvas, set_bbox!

end  # Draw

module Gmsh
using EMI.Draw
include("emi_gmsh.jl")

export geo_code

end  # Gmsh 

end # EMI
