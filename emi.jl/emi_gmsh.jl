# Convert to primitives  -----------------------------------------------------------------
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

# Physical points, and curves defined in terms of LOCAL point numbering
function primitives(shape::ClosedPolygon)
    # Point part
    points = shape.points

    npoints = length(points)
    pindices = 1:npoints  # Local
    lines = [zip(pindices[2:end], pindices[1:end-1])..., (first(pindices), last(pindices))]
    lines = Dict(:Line => lines)

    (points, lines)
end

function primitives(shape::Circle)
    points = [shape.center,
              shape.center + Point(shape.radius, 0),
              shape.center + Point(0, shape.radius),
              shape.center + Point(-shape.radius, 0),
              shape.center + Point(0, -shape.radius)]
    lines = Dict(:Circle => [(2, 1, 3), (3, 1, 4), (4, 1, 5), (5, 1, 2)])

    (points, lines)
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
    # NOTE: first, center, last, orientation
    lines = Dict(:Ellipse => [(2, 1, 3, 1), (4, 1, 3, -1), (4, 1, 5, 1), (2, 1, 5, -1)])

    (points, lines)
end

# Convert to gmsh ------------------------------------------------------------------------

function gmsh_script(canvas::Canvas, size::Real=1.)
    counter = 1  # Of any GMSH entities

    gmsh_code = Vector{AbstractString}(["SIZE = $(size);\n"])
    loops = Vector{Int}()  # Each shape is a closed loop, their collection together with the
                           # canvas loop defines the outer domain

    # We first write shapes
    for (tag, shape) in enumerate(canvas.shapes)
        points, lines = primitives(shape)
        # Will need to convert entities from local to global
        local_to_global = collect(counter:counter+length(points)-1)
        # Code for points
        points_gmsh, counter = gmsh_script(points, local_to_global, counter)
        # In case of lines we want a loop as well; for loops and to define physical ...
        lines_gmsh, shape_loop, counter = gmsh_script(lines, local_to_global, counter)
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
                           lines_gmsh,
                           loop_gmsh,
                           surf_gmsh,
                           physurf_gmsh, 
                           physline_gmsh,
                           "//--------------------\n"], "\n")  # For readability
        push!(gmsh_code, shape_code)
    end

    tag = length(canvas.shapes) + 1
    # Now the canvas
    points, lines = primitives(canvas.bbox)
    local_to_global = collect(counter:counter+length(points)-1)
    
    points_gmsh, counter = gmsh_script(points, local_to_global, counter)
    
    lines_gmsh, shape_loop, counter = gmsh_script(lines, local_to_global, counter)
    loop_gmsh = "Line Loop($(counter)) = {$(join(map(x -> "$(x)", shape_loop), ", "))};"
    counter += 1

    surf_gmsh = "Plane Surface($(counter)) = {$(counter-1), $(join(map(x -> "$(x)", loops), ", "))};"
    physurf_gmsh = "Physical Surface($(tag)) = {$(counter)};"
    physline_gmsh = "Physical Line($(tag)) = {$(join(map(x-> "$(x)", shape_loop), ", "))};"

    canvas_code = join([points_gmsh,
                        lines_gmsh,
                        loop_gmsh,
                        surf_gmsh,
                        physurf_gmsh, 
                        physline_gmsh,
                        "//--------------------\n"], "\n")
    push!(gmsh_code, canvas_code)

    # Done
    join(gmsh_code, "\n")
end

function gmsh_script(points::Vector{Point}, local_to_global::Vector{Int}, counter::Int)
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

function gmsh_script(lines::Dict, local_to_global::Vector{Int}, counter::Int)
    # Lines are written differently based on straight, circle, etc..
    lines_gmsh = Vector{AbstractString}()
    shape_loop = Vector{Int}() 
    for line_type in keys(lines)
        # Straight, 2 points
        if line_type == :Line
            for line in lines[line_type]
                v0, v1 = local_to_global[first(line)], local_to_global[last(line)]
                str = "Line($(counter)) = {$(v0), $(v1)};"
                push!(lines_gmsh, str)
                push!(shape_loop, counter)
                counter +=1
            end
        end

        # Circle, 3 points as begin, center, end
        if line_type == :Circle
            for line in lines[line_type]
                v0, v1 = local_to_global[first(line)], local_to_global[last(line)]
                center = local_to_global[line[2]]
                str = "Circle($(counter)) = {$(v0), $(center), $(v1)};"
                push!(lines_gmsh, str)
                push!(shape_loop, counter)
                counter +=1
            end
        end

        # Ellipse, 4 points as begin, center, end; sign
        if line_type == :Ellipse
            for line in lines[line_type]
                v0, center, v1 = local_to_global[line[1]], local_to_global[line[2]], local_to_global[line[3]]
                sgn = line[4]

                str = "Ellipse($(counter)) = {$(v0), $(center), $(v0), $(v1)};"
                push!(lines_gmsh, str)
                push!(shape_loop, sgn*counter)
                counter +=1
            end
        end
    end
    lines_gmsh = join(lines_gmsh, "\n")

    (lines_gmsh, shape_loop, counter)
end

# FIXME: Obligue
# FIXME: do we want cotrol over sizes?
