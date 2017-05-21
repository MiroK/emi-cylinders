import Base: length, first, last, getindex

abstract IndexRepr

immutable LineRepr <: IndexRepr
    id::NTuple{2, Int}
end

immutable CircleArcRepr <: IndexRepr
    id::NTuple{3, Int}
end

immutable EllipseArcRepr <: IndexRepr
    id::NTuple{3, Int}
end

for (T, size) in zip((:LineRepr, :CircleArcRepr, :EllipseArcRepr), (2, 3, 3))
    eval(quote
            first(r::$T) = first(r.id)
            last(r::$T) = last(r.id)
            getindex(r::$T, i::Int) = r.id[i]
         end)
end

# Convert to primitives  -----------------------------------------------------------------

function primitives(l::Line)
    points = [l.p0, l.p1]
    indices = LineRepr((1, 2))
    (points, indices)
end

function primitives(c::CircleArc)
    points = [c.p0, c.center, c.p1]
    indices = CircleArcRepr((1, 2, 3))
    (points, indices)
end

function primitives(c::EllipseArc)
    points = [c.p0, c.center, c.p1]
    indices = EllipseArcRepr((1, 2, 3))
    (points, indices)
end

function primitives(shape::BoundingBox)
    dx = shape.ur.x - shape.ll.x
    dy = shape.ur.y - shape.ll.y

    points = [shape.ll, shape.ll + Point(dx, 0), shape.ur, shape.ll + Point(0, dy)]
    primitives(ClosedPolygon(points))
end

function primitives(shape::Square)
    points = [shape.ll,
              shape.ll + Point(shape.size, 0),
              shape.ll + Point(shape.size, shape.size),
              shape.ll + Point(0, shape.size)]

    primitives(ClosedPolygon(points))
end

function primitives(shape::Rectangle)
    points = [shape.ll,
              shape.ll + Point(shape.size_x, 0),
              shape.ll + Point(shape.size_x, shape.size_y),
              shape.ll + Point(0, shape.size_y)]

    primitives(ClosedPolygon(points))
end

# Physical points, curves defined in terms of LOCAL point numbering, orientation
function primitives(shape::ClosedPolygon)
    # Point part
    points = shape.points

    npoints = length(points)
    pindices = 1:npoints  # Local
    curves = map(LineRepr, zip([pindices[2:end]..., first(pindices)],
                               [pindices[1:end-1]..., last(pindices)]))

    (points, [(c, 1) for c in curves])
end

function primitives(shape::Circle)
    points = [shape.center,
              shape.center + Point(shape.radius, 0),
              shape.center + Point(0, shape.radius),
              shape.center + Point(-shape.radius, 0),
              shape.center + Point(0, -shape.radius)]

    curves = map(CircleArcRepr, [(2, 1, 3), (3, 1, 4), (4, 1, 5), (5, 1, 2)])

    (points, [(c, 1) for c in curves])
end

function primitives(shape::Ellipse)
    if shape.size_x > shape.size_y
        points = [shape.center,
                  shape.center + Point(shape.size_x, 0),
                  shape.center + Point(0, shape.size_y),
                  shape.center + Point(-shape.size_x, 0),
                  shape.center + Point(0, -shape.size_y)]
    else
        points = [shape.center,
                  shape.center + Point(0, shape.size_y),
                  shape.center + Point(shape.size_x, 0),
                  shape.center + Point(0, -shape.size_y),
                  shape.center + Point(-shape.size_x, 0)]
    end
    curves = map(EllipseArcRepr, [(2, 1, 3), (4, 1, 3), (4, 1, 5), (2, 1, 5)])
    orientation = [1, -1, 1, -1]

    (points, collect(zip(curves, orientation)))
end

function primitives(shape::Loop)
    curves = shape.curves
    orientation = shape.orientation
    @assert first(orientation) == 1

    points, indices = primitives(first(curves))
    # First curve for sure can add all its points
    all_points = Vector{Point}(points)
    # And also its indices
    all_indices = Vector{IndexRepr}([indices])
    # The connection over last point

    link = last(points)
    link_index = last(indices)
    for (curve, sign) in zip(curves[2:end-1], orientation[2:end-1])
        counter = length(all_points)
        # Connection over first
        points, indices = primitives(curve)
        new_indices = [counter+i for i in 1:length(points)-1]
        if sign == 1
            all_points = vcat(all_points, points[2:end])
            repr = typeof(indices)(tuple(link_index, new_indices...))
            link = last(curve)
            link_index = last(repr)
        else
            all_points = vcat(all_points, points[1:end-1])
            repr = typeof(indices)(tuple(new_indices..., link_index))
            link = first(curve)
            link_index = first(repr)
        end
        push!(all_indices, repr)
    end

    # FIXME: Not sure about this part
    # The last curve closes the loop so if 1/-1 oriented we have seen its last/first point
    # That point was seen as 1, last is given by counter
    counter = length(all_points)
    curve, sign = last(curves), last(orientation)
    points, indices = primitives(curve)
    
    all_points = vcat(all_points, points[2:end-1])  # Only interior points
    new_indices = [counter+i for i in 1:length(points)-2]
    if sign == 1
        repr = typeof(indices)(tuple(link_index, new_indices..., 1))
    else
        repr = typeof(indices)(tuple(1, new_indices..., link_index))
    end
    push!(all_indices, repr)

    (all_points, collect(zip(all_indices, orientation)))
end

# Convert to gmsh ------------------------------------------------------------------------

#function mesh(canvas::Canvas, file::AbstractString, size::Real)
#    base, ext = splitext(file)
#
#    geo_file = "$(base).geo"
#    if ext == ".msh"
#        @assert gmsh_script(canvas, "$(geo_file)", size) > 0
#        run(`gmsh -2 $(geo_file) -o $(file)`)
#    end
#end


function gmsh_script(canvas::Canvas, size::Real)
    count = write(STDOUT, gmsh_code(canvas, size))
    "//$(count)"
end

function gmsh_script(canvas::Canvas, file::AbstractString, size::Real)
    @assert last(splitext(file)) == ".geo" 
    write(open(file, "w"), gmsh_code(canvas, size))
end

function gmsh_code(canvas::Canvas, size::Real)
    counter = 1  # Of any GMSH entities

    code = Vector{AbstractString}(["SIZE = $(size);\n"])
    loops = Vector{Int}()  # Each shape is a closed loop, their collection together with the
                           # canvas loop defines the outer domain

    # We first write shapes
    for (tag, shape) in enumerate(canvas.shapes)
        points, curves = primitives(shape)
        # Will need to convert entities from local to global
        local_to_global = collect(counter:counter+length(points)-1)
        # Code for points
        points_gmsh, counter = gmsh_code(points, local_to_global, counter)
        # In case of lines we want a loop as well; for loops and to define physical ...
        curves_gmsh, shape_loop, counter = gmsh_code(curves, local_to_global, counter)
        # Line Loop is another entity
        loop_gmsh = "Line Loop($(counter)) = {$(join(map(x -> "$(x)", shape_loop), ", "))};"
        counter += 1
        # Plane surface in terms of line loop
        surf_gmsh = "Plane Surface($(counter)) = {$(counter-1)};"
        # Physical, not entities
        physurf_gmsh = "Physical Surface($(tag)) = {$(counter)};"
        physline_gmsh = "Physical Line($(tag)) = {$(join(map(x-> "$(x)", shape_loop), ", "))};"

        loops = vcat(loops, counter-1)  # Add to 'inner' boundary
        # Prepare for next round
        counter += 1  # Counter so that it is ready to use
        shape_code = join([points_gmsh,
                           curves_gmsh,
                           loop_gmsh,
                           surf_gmsh,
                           physurf_gmsh, 
                           physline_gmsh,
                           "//--------------------\n"], "\n")  # For readability
        push!(code, shape_code)
    end
    # Now the canvas
    tag = length(canvas.shapes) + 1

    points, curves = primitives(canvas.bbox)
    local_to_global = collect(counter:counter+length(points)-1)
    
    points_gmsh, counter = gmsh_code(points, local_to_global, counter)
    
    curves_gmsh, shape_loop, counter = gmsh_code(curves, local_to_global, counter)
    loop_gmsh = "Line Loop($(counter)) = {$(join(map(x -> "$(x)", shape_loop), ", "))};"
    counter += 1

    surf_gmsh = "Plane Surface($(counter)) = {$(counter-1), $(join(map(x -> "$(x)", loops), ", "))};"
    physurf_gmsh = "Physical Surface($(tag)) = {$(counter)};"
    physline_gmsh = "Physical Line($(tag)) = {$(join(map(x-> "$(x)", shape_loop), ", "))};"

    canvas_code = join([points_gmsh,
                        curves_gmsh,
                        loop_gmsh,
                        surf_gmsh,
                        physurf_gmsh, 
                        physline_gmsh,
                        "//--------------------\n"], "\n")
    push!(code, canvas_code)

    # Done
    join(code, "\n")
end

function gmsh_code(points::Vector{Point}, local_to_global::Vector{Int}, counter::Int)
    # Code for points
    points_gmsh = Vector{AbstractString}()
    
    for (local_index, point) in enumerate(points)
        str = "Point($(local_to_global[local_index])) = {$(point.x), $(point.y), 0., SIZE};"
        push!(points_gmsh, str)
        counter +=1
    end
    points_gmsh = join(points_gmsh, "\n")

    (points_gmsh, counter)
end

function gmsh_code(line::Tuple{LineRepr, Int}, local_to_global::Vector{Int}, counter::Int)
    line, orientation = line
    v0, v1 = local_to_global[first(line)], local_to_global[last(line)]
    line_gmsh = "Line($(counter)) = {$(v0), $(v1)};"
    loop_piece = orientation*counter

    counter += 1
    (line_gmsh, loop_piece, counter)
end

function gmsh_code(curve::Tuple{CircleArcRepr, Int}, local_to_global::Vector{Int}, counter::Int)
    curve, orientation = curve
    v0, v1 = local_to_global[first(curve)], local_to_global[last(curve)]
    center = local_to_global[curve[2]]
    curve_gmsh = "Circle($(counter)) = {$(v0), $(center), $(v1)};"
    loop_piece = orientation*counter

    counter += 1
    (curve_gmsh, loop_piece, counter)
end

function gmsh_code(curve::Tuple{EllipseArcRepr, Int}, local_to_global::Vector{Int}, counter::Int)
    curve, orientation = curve
    v0, center, v1 = local_to_global[curve[1]], local_to_global[curve[2]], local_to_global[curve[3]]
    curve_gmsh = "Ellipse($(counter)) = {$(v0), $(center), $(v0), $(v1)};"
    loop_piece = orientation*counter

    counter += 1
    (curve_gmsh, loop_piece, counter)
end

function gmsh_code{T<:IndexRepr}(curves::Vector{Tuple{T, Int}}, local_to_global::Vector{Int}, counter::Int)
    curves_gmsh = Vector{AbstractString}()
    shape_loop = Vector{Int}() 
    
    for curve in curves
        curve_gmsh, loop_piece, counter = gmsh_code(curve, local_to_global, counter)
        push!(curves_gmsh, curve_gmsh)
        push!(shape_loop, loop_piece)
    end
    curves_gmsh = join(curves_gmsh, "\n")

    (curves_gmsh, shape_loop, counter)
end
