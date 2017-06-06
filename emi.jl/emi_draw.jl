import Base: ==, abs, dot, -, +, cross, first, last, isempty, *, size

#########################
# Drawing simple shapes #
#########################

"""Point in R^2"""
immutable Point
    x::Real
    y::Real
end
Point(p::Tuple{Real, Real}) = Point(p[1], p[2])
==(p::Point, q::Point) = (p.x == q.x) && (p.y == q.y)
-(p::Point, q::Point) = Point(p.x-q.x, p.y-q.y)
+(p::Point, q::Point) = Point(p.x+q.x, p.y+q.y)
*(p::Point, a::Number) = Point(a*p.x, a*p.y)
*(a::Number, p::Point) = Point(a*p.x, a*p.y)

dot(p::Point, q::Point) = p.x*q.x + p.y*q.y
abs(p::Point) = sqrt(dot(p, p))
cross(p::Point, q::Point) = p.x*q.y - p.y*q.x

first(p::Point) = p.x
last(p::Point) = p.y

# ----------------------------------------------------------------------------------------

# We are after making closed polygons. These are some primitives which constitute them
abstract Curve

immutable Line <: Curve
    p0::Point
    p1::Point

    function Line(p0, p1)
        @assert !(p0 == p1)
        new(p0, p1)
    end
end
+(l::Line, vec::Point) = Line(l.p0+vec, l.p1+vec)
first(l::Line) = l.p0
last(l::Line) = l.p1

immutable CircleArc <: Curve
    p0::Point
    center::Point
    p1::Point
    # NOTE: order of p0, p1 translates to orientation
    function CircleArc(p0, center, p1)
        @assert !(p0 == center) && !(p1 == center) && !(p0 == p1)
        new(p0, center, p1)
    end
end
+(c::CircleArc, vec::Point) = CircleArc(c.p0+vec, c.center+vec, c.p1+vec)
first(c::CircleArc) = c.p0
last(c::CircleArc) = c.p1

immutable EllipseArc <: Curve
    p0::Point
    center::Point
    p1::Point
    # NOTE: this is a defintion similar to GMSH with the condition that start is at major
    # point; so p0(start), center, p0, p1(end) is how you'd do it in GMSH. Beware that this
    # convenction enforces orientation
    function EllipseArc(p0, center, p1)
        @assert !(p0 == center) && !(p1 == center) && !(p0 == p1)
        @assert abs(p0-center) > abs(p1-center)
        new(p0, center, p1)
    end
end
+(c::EllipseArc, vec::Point) = EllipseArc(c.p0+vec, c.center+vec, c.p1+vec)
first(c::EllipseArc) = c.p0
last(c::EllipseArc) = c.p1

# ----------------------------------------------------------------------------------------

"""Shapes of cells"""
abstract Shape

immutable Circle <: Shape
    center::Point
    radius::Real

    function Circle(center, radius)
        @assert radius > 0
        new(center, radius)
    end
end
+(shape::Circle, vec::Point) = Circle(shape.center+vec, shape.radius)

immutable Ellipse <: Shape
    center::Point
    size_x::Real
    size_y::Real

    function Ellipse(center, x_axis, y_axis)
        @assert x_axis > 0 && y_axis > 0
        new(center, x_axis, y_axis)
    end
end
+(shape::Ellipse, vec::Point) = Ellipse(shape.center+vec, shape.size_x, shape.size_y)

immutable Rectangle <: Shape
    ll::Point
    size_x::Real
    size_y::Real

    function Rectangle(ll, size_x, size_y)
        @assert size_x > 0 && size_y > 0
        new(ll, size_x, size_y)
    end
end
+(shape::Rectangle, vec::Point) = Rectangle(shape.ll+vec, shape.size_x, shape.size_y)

immutable Square <: Shape
    ll::Point
    size::Real

    function Square(ll, size)
        @assert size > 0
        new(ll, size)
    end
end
+(shape::Square, vec::Point) = Square(shape.ll+vec, shape.size)

immutable ClosedPolygon <: Shape
    points::Vector{Point}

    function ClosedPolygon(points)
        @assert length(points) > 2
        @assert first(points) != last(points)
        new(points)
    end
end
+(shape::ClosedPolygon, vec::Point) = ClosedPolygon(map(p -> p+vec, shape.points))

Triangle(v0::Point, v1::Point, v2::Point) = ClosedPolygon([v0, v1, v2])

# This is a equilateral N-sided polygon with given volume
function NGon(N::Int, V::Real, x0::Real=0, y0::Real=0)
    @assert N > 2 && V > 0
    R = sqrt(2*V/N/sin(2*pi/N))
    angles = 2*pi/N*(1:N)
    ClosedPolygon(map(Point, zip(x0 + R*cos(angles), y0 + R*sin(angles))))
end

function NGonR(N::Int, R::Real, x0::Real=0, y0::Real=0)
    @assert N > 2 && R > 0
    angles = 2*pi/N*(1:N)
    ClosedPolygon(map(Point, zip(x0 + R*cos(angles), y0 + R*sin(angles))))
end

function make_path{T<:Curve}(curves::Vector{T})
    curve = first(curves)
    orientation = [1]
    path = Vector{Curve}([curve])

    for index in 2:length(curves)
        next = curves[index]
        linked = false
        # The first curve is orientated and the curves should follow each other
        # So either they agree on orientation
        link = (last(orientation) == 1) ? last(curve) : first(curve)

        if link == first(next)
            push!(orientation, 1)
            linked = true
        end
        # Or they don't
        if link == last(next)
            push!(orientation, -1)
            linked = true
        end
        # But they must be ordered
        @assert linked "Linking $(index): $(last(curve)), $(first(next)), $(last(next)) | $(next)"
        push!(path, next)
        curve = next
    end
    (path, orientation)
end

immutable Loop <: Shape
    curves::Vector{Curve}
    orientation::Vector{Int}

    function Loop(curves)
        @assert !isempty(curves)
        
        curves, orientation = make_path(curves)
        # Now see if it is closed
        last(orientation) == 1 && @assert first(first(curves)) == last(last(curves))
        last(orientation) == -1 && @assert first(first(curves)) == first(last(curves))

        new(curves, orientation)
    end
end
+(shape::Loop, vec::Point) = Loop(map(p -> p+vec, shape.curves))

# The points should be unique
function is_degenerate(poly::ClosedPolygon)
    points = poly.points
    n = length(poly.points)
    for i in 1:n
        pointi = points[i]
        for j in i+1:n-1
            pointj = points[j]
            pointi == pointj && return true
        end
    end
    false
end

# Checking convexity (not used anywhere, just for fun)
function is_convex(poly::ClosedPolygon)
    points = poly.points
    n = length(points)

    for k in 1:n
        prev = (k == 1) ? n : k-1
        next = (k == n) ? 1 : k+1

        outflow = points[next]-points[k]
        inflow = points[k]-points[prev]

        sign(cross(outflow, inflow)) > 0 && return false
    end
    true
end

# ----------------------------------------------------------------------------------------

"""Smallest rectangle that contains the shape"""
immutable BoundingBox
    ll::Point
    ur::Point
    
    function BoundingBox(ll, ur)
        @assert first(ll) <= first(ur) 
        @assert last(ll) <= last(ur)
        new(ll, ur)
    end
end
size(bbox::BoundingBox, dim::Int) = (dim == 1) ? (bbox.ur.x-bbox.ll.x) : (bbox.ur.y-bbox.ll.y)

BoundingBox(l::Line) = BoundingBox([l.p0, l.p1])
BoundingBox(c::CircleArc) = BoundingBox([c.p0, c.p1])
BoundingBox(c::EllipseArc) = BoundingBox([c.p0, c.p1])

BoundingBox(shape::Circle) = BoundingBox(shape.center - Point(shape.radius, shape.radius),
                                         shape.center + Point(shape.radius, shape.radius))

BoundingBox(shape::Ellipse) = BoundingBox(shape.center - Point(shape.size_x, shape.size_y), 
                                          shape.center + Point(shape.size_x, shape.size_y))
                                          
BoundingBox(shape::Rectangle) = BoundingBox(shape.ll,
                                            shape.ll + Point(shape.size_x, shape.size_y))

BoundingBox(shape::Square) = BoundingBox(shape.ll,
                                         shape.ll + Point(shape.size, shape.size))

BoundingBox(points::Vector{Point}) = BoundingBox(Point(minimum(map(first, points)),
                                                      minimum(map(last, points))),
                                                Point(maximum(map(first, points)),
                                                      maximum(map(last, points))))

BoundingBox(shape::ClosedPolygon) = BoundingBox(shape.points)

BoundingBox(loop::Loop) = BoundingBox(map(BoundingBox, loop.curves))

BoundingBox{T<:Shape}(shapes::Vector{T}) = BoundingBox(map(BoundingBox, shapes))

BoundingBox(b::BoundingBox) = b
BoundingBox(b::BoundingBox, B::BoundingBox) = BoundingBox([b.ll, b.ur, B.ll, B.ur])

function BoundingBox(boxes::Vector{BoundingBox})
    length(boxes) == 1 && return first(boxes)
    # Otherwise, don't want recursion
    BoundingBox([[b.ll for b in boxes]..., [b.ur for b in boxes]...])
end

# Collision between bounding boxes
function collides(b::BoundingBox, B::BoundingBox, tol=1E-13)
    const ll, ur, LL, UR = b.ll, b.ur, B.ll, B.ur
    !(((ur.x+tol < LL.x) || (UR.x + tol < ll.x)) || ((ur.y+tol < LL.y) || (UR.y + tol < ll.y)))
end

collides{T<:Shape}(b::BoundingBox, shape::T) = collides(b, BoundingBox(shape))
collides{T<:Shape}(shape::T, B::BoundingBox) = collides(BoundingBox(shape), B)
collides{T<:Shape, S<:Shape}(b::T, B::S) = collides(BoundingBox(b), BoundingBox(B))

# ----------------------------------------------------------------------------------------

"""Canvas is where we draw the shapes"""
type Canvas
    shapes::Vector{Shape}
    bbox::BoundingBox
end
Canvas{T<:Shape}(shapes::Vector{T})=Canvas(shapes, BoundingBox(shapes))
Canvas() = Canvas(Vector{Shape}(), BoundingBox(Point(-Inf, -Inf), Point(Inf, Inf)))

collides{T<:Shape}(canvas::Canvas, shape::T) = collides(canvas.bbox, BoundingBox(shape))
isempty(canvas::Canvas) = length(canvas.shapes) == 0

"""Add shape - modifies canvas, recomputed bbox. Only allowed for noncolliding"""
function +(canvas::Canvas, shape::Shape)
    if !isempty(canvas)
        @assert !collides(canvas, shape)
        canvas.bbox = BoundingBox(canvas.bbox, BoundingBox(shape))
    else
        canvas.bbox = BoundingBox(shape)
    end
    push!(canvas.shapes, shape)
    canvas
end

function +(canvas::Canvas, shape::ClosedPolygon)
    if !isempty(canvas)
        @assert !is_degenerate(shape)
        @assert !collides(canvas, shape)
        canvas.bbox = BoundingBox(canvas.bbox, BoundingBox(shape))
    else
        canvas.bbox = BoundingBox(shape)
    end
    push!(canvas.shapes, shape)
    canvas
end

BoundingBox(c::Canvas, C::Canvas) = BoundingBox(c.bbox, C.bbox)

function +(c::Canvas, C::Canvas)
    isempty(C) && return c
    # A copy
    if isempty(c)
        c.shapes = C.shapes
        c.bbox = C.bbox
        return c
    end

    c.bbox = BoundingBox(c, C)
    c.shapes = vcat(c.shapes, C.shapes)
    c
end

"""For later drawing increase the bounding rectangle"""
function set_bbox!(canvas::Canvas, ll::Point, ur::Point)
    @assert ll.x < canvas.bbox.ll.x && ll.y < canvas.bbox.ll.y
    @assert ur.x > canvas.bbox.ur.x && ur.y > canvas.bbox.ur.y
    canvas.bbox = BoundingBox(ll, ur)
end

function set_bbox!(canvas::Canvas, padx::Real, pady::Real)
    @assert padx > 0 && pady > 0
    ll = canvas.bbox.ll - Point(padx, pady)
    ur = canvas.bbox.ur + Point(padx, pady)
    canvas.bbox = BoundingBox(ll, ur)
end

#############################
# Drawing composites/tissue #
#############################

