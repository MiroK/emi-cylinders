import Base: length, first, last, getindex

abstract IndexRepr

immutable LineRepr <: IndexRepr
    id::NTuple{2, Int}
end
interior_indices(::LineRepr) = ()

immutable CircleArcRepr <: IndexRepr
    id::NTuple{3, Int}
end
interior_indices(::CircleArcRepr) = (2,)

immutable EllipseArcRepr <: IndexRepr
    id::NTuple{3, Int}
end
interior_indices(::EllipseArcRepr) = (2,)

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
    # NOTE: loop as it comes here has curves defined in terms of points which can be not can
    # coincide. Here we try to keep only the uniqe points - the check for duplicates,
    # hoever, involves only the interior points of curves.

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
        # First, interior, Last
        if sign == 1
            # First, is known 
            global_indices = [link_index]
            # Need to check interior points
            for interior_id in interior_indices(indices)
                point = points[interior_id]
                index = findfirst(q -> q == point, all_points)
                # Not seen
                if index == 0
                    push!(all_points, point)

                    counter += 1
                    index = counter
                end
                # Seen with index
                push!(global_indices, index)
            end
            # The last gut is new
            link = last(points)
            push!(all_points, link)
            counter += 1

            link_index = counter
            push!(global_indices, link_index)

            repr = typeof(indices)(tuple(global_indices...))
        else
            # Last, is known so remember to add it later
            # First is new
            link = first(points)
            push!(all_points, link)
            counter +=1

            global_indices = [counter]
            # Need to check interior points
            for interior_id in interior_indices(indices)
                point = points[interior_id]
                index = findfirst(q -> q == point, all_points)
                # Not seen
                if index == 0
                    push!(all_points, point)

                    counter += 1
                    index = counter
                end
                # Seen with index
                push!(global_indices, index)
            end
            # The last point is there already as link_index
            push!(global_indices, link_index)

            repr = typeof(indices)(tuple(global_indices...))
            link_index = first(repr)
        end
        push!(all_indices, repr)
    end

    # The last curve closes the loop so if 1/-1 oriented we have seen its last/first point
    # That point was seen as 1, last is given by link_index
    counter = length(all_points)
    curve, sign = last(curves), last(orientation)
    points, indices = primitives(curve)
   
    global_indices = []
    fi, li = (sign == 1) ? (link_index, 1) : (1, link_index)
        
    push!(global_indices, fi)
    # Interior guys
    for interior_id in interior_indices(indices)
        point = points[interior_id]
        index = findfirst(q -> q == point, all_points)
        # Not seen
        if index == 0
            push!(all_points, point)

            counter += 1
            index = counter
        end
        # Seen with index
        push!(global_indices, index)
    end
    push!(global_indices, li)
    
    repr = typeof(indices)(tuple(global_indices...))
    push!(all_indices, repr)

    (all_points, collect(zip(all_indices, orientation)))
end

# Convert to gmsh ------------------------------------------------------------------------

function mesh(canvas::Canvas, size::Vector, file::AbstractString)
    base, ext = splitext(file)
    for (index, size_i) in enumerate(size)
        file_i = base * "$(index)" * ext
        mesh(canvas, size_i, file_i)
    end
end

function mesh(canvas::Canvas, size::Real, file::AbstractString)
    base, ext = splitext(file)

    geo_file = "$(base).geo"
    if ext == ".msh"
        @assert gmsh_script(canvas, "$(geo_file)", size)
        gmsh = success(`gmsh -2 $(geo_file) -o $(file)`)
        return gmsh && success(`rm $(geo_file)`)
    end
    # XML mesh, mesh itself & physical & facet
    msh_file = "$(base).msh"
    if ext == ".xml"
        @assert mesh(canvas, size, msh_file)
        convert = success(`dolfin-convert $(msh_file) $(file)`)
        return convert && success(`rm $(msh_file)`)
    end
    # HDF5
    xml_file = "$(base).xml"
    if ext == ".h5"
        @assert mesh(canvas, size, xml_file)

const python = "
from dolfin import Mesh, MeshFunction, mpi_comm_world, HDF5File;
import os;

xml_mesh = '$(xml_file)'
base, ext = os.path.splitext(xml_mesh);
assert os.path.exists(xml_mesh);

xml_volumes = '$(base)_physical_region.xml';
assert os.path.exists(xml_volumes);

xml_facets = '$(base)_facet_region.xml';
has_facet_regions = os.path.exists(xml_facets);

mesh = Mesh(xml_mesh);
cell_f = MeshFunction('size_t', mesh, xml_volumes);

h5_file = '$(base).h5';
out = HDF5File(mesh.mpi_comm(), h5_file, 'w');
out.write(mesh, '/mesh');
out.write(cell_f, '/cell_markers');

if has_facet_regions: out.write(MeshFunction('size_t', mesh, xml_facets), '/facet_markers')
#"
        convert = success(`python -c "$(python)"`)
        xmls = [xml_file, "$(base)_physical_region.xml", "$(base)_facet_region.xml"]
        return convert && all(success(`rm $(file)`) for file in filter(isfile, xmls))
    end
end

gmsh_script(canvas::Canvas, size::Real) = write(STDOUT, gmsh_code(canvas, size))

function gmsh_script(canvas::Canvas, file::AbstractString, size::Real)
    @assert last(splitext(file)) == ".geo" 
    f = open(file, "w")
    count = write(f, gmsh_code(canvas, size))
    close(f)
    count > 0
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
        # NOTE, all stuff below is what gmsh_code of a curve should do .... Anyways, there
        # will be rewrite of this because of 3d
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
