import Base: ==, abs, dot, -, +, cross, first, last, isempty, *

#########################
# Drawing simple shapes #
#########################

"""Point in R^2"""
immutable Point
    x::Real
    y::Real
end
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

    function CircleArc(p0, center, p1)
        @assert !(p0 == center) && !(p1 == center) && !(p0 == p1)
        new(p0, center, p1)
    end
end
+(c::CircleArc, vec::Point) = CircleArc(c.p0+vec, c.center+vec, c.p1+vec)
first(c::CircleArc) = c.p0
last(c::CircleArc) = c.p1

# FIXME: Ellipse here

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
        new(center, size_x, size_y)
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

immutable CompositeLoop <: Shape
    curves::Vector{Curve}

    function CompositeLoop(curves)
        @assert !isempty(curves)
        # The loop must be closed
        length(curves) == 1 && @assert first(first(curves)) == last(first(curves))
        @assert first(first(curves)) == last(last(curves))

        # The line form a connected path, first of this is last of prev
        for k in 2:length(curves)
            @assert first(curves[k]) == last(curves[k-1])
        end
        
        new(curves)
    end
end

# The points should be unique
function is_degenerate(poly::ClosedPolygon)
    points = poly.points
    n = length(poly.points)
    for i in 1:n
        pointi = points[i]
        for j in i:n-1
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

# FIXME: CircleArc, EllipseArc, form a closed loop curve CompositeLoop

# ----------------------------------------------------------------------------------------

"""Smallest rectangle that contains the shape"""
immutable BoundingBox
    ll::Point
    ur::Point
    
    function BoundingBox(ll, ur)
        @assert first(ll) < first(ur) 
        @assert last(ll) < last(ur)
        new(ll, ur)
    end
end

BoundingBox(shape::Circle) = BoundingBox(shape.center - Point(shape.radius, shape.radius),
                                         shape.center + Point(shape.radius, shape.radius))

BoundingBox(shape::Ellipse) = BoundingBox(shape.center - Point(shape.size_x, shape.size_y), 
                                          shape.center + Point(shape.size_x, shape.size_y))
                                          
BoundingBox(shape::Rectangle) = BoundingBox(shape.ll,
                                            shape.ll + Point(shape.size_x, shape.size_y))

BoundingBox(shape::Square) = BoundingBox(shape.ll,
                                         shape.ll + Point(shape.size, shape.size))

BoundingBox(shape::ClosedPolygon) = BoundingBox(Point(minimum(map(first, shape.points)),
                                                      minimum(map(last, shape.points))),
                                                Point(maximum(map(first, shape.points)),
                                                      maximum(map(last, shape.points))))

"""Collision between bounding boxes"""
function collides(b::BoundingBox, B::BoundingBox, tol=1E-13)
    const ll, ur, LL, UR = b.ll, b.ur, B.ll, B.ur
    !(((ur.x+tol < LL.x) || (UR.x + tol < ll.x)) || ((ur.y+tol < LL.y) || (UR.y + tol < ll.y)))
end

collides{T<:Shape}(b::BoundingBox, shape::T) = collides(b, BoundingBox(shape))
collides{T<:Shape}(shape::T, B::BoundingBox) = collides(BoundingBox(shape), B)
collides{T<:Shape, S<:Shape}(b::T, B::S) = collides(BoundingBox(b), BoundingBox(B))

BoundingBox(b::BoundingBox, B::BoundingBox) = BoundingBox(Point(min(b.ll.x, B.ll.x),
                                                                min(b.ll.y, B.ll.y)),
                                                          Point(max(b.ur.x, B.ur.x),
                                                                max(b.ur.y, B.ur.y)))

function BoundingBox{T<:Shape}(shapes::Vector{T})
    length(shapes) == 1 && return BoundingBox(first(shapes))
    length(shapes) == 2 && return BoundingBox(map(BoundingBox, shapes)...)
    # Otherwise
    BoundingBox(BoundingBox(shapes[1:2]), BoundingBox(shapes[3:end]))
end

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
    @assert ll.x < canvas.ll.x && ll.y < canvas.ll.y
    @assert ur.x > canvas.ur.x && ur.y > canvas.ur.y
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
    println(shapes)
    Canvas(shapes)
end

# FIXME: concentric circle, concentric triangles, tiling
