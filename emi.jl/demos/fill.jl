include("../emi.jl")
using EMI
using EMI.Draw, EMI.Gmsh

"""Make nxm grid of shapes which are separated by given spacing"""
function fill_canvas(shape::Shape, counts::Tuple{Int, Int}, spacing::Tuple{Real, Real})
    # Compute the transform
    bb = BoundingBox(shape)
    dx = (first(bb.ur) - first(bb.ll)) + first(spacing)
    dy = (last(bb.ur) - last(bb.ll)) + last(spacing)
    shift_x = Point(dx, 0)
    shift_y = Point(0, dy)
   
    n, m = counts
    model = shape
    shapes = [shape]
    for j in 1:m
        # Fill horizontal o -> o -> o 
        for i in 1:n-1
            push!(shapes, model + i*shift_x)
        end
        # Fill vertical
        model = model + shift_y
        j < m && push!(shapes, model)
    end
    Canvas(shapes)
end

# shape = ClosedPolygon([Point(0, 0), Point(1, 0), Point(1, 1)])

const R, dx, dy = 1, 0.2, 0.3
shape = Loop([Line(Point(0, R+dy+dx), Point(0, R+dy)),
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

canvas = fill_canvas(shape, (3, 5), (0.2, 0.1))

set_bbox!(canvas, 0.2, 0.3)
gmsh_script(canvas, 0.2)
