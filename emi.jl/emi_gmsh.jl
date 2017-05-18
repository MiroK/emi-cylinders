# Flat -------------------------------------------------------------------------------
function geo_code(shape::Square, start::Int=1, tag::Int=0)
    points = [shape.ll,
              shape.ll + Point(shape.size, 0),
              shape.ll + Point(shape.size, shape.size),
              shape.ll + Point(0, shape.size)]
    geo_code(ClosedPolygon(points), start, tag)
end

function geo_code(shape::Rectangle, start::Int=1, tag::Int=0)
    points = [shape.ll,
              shape.ll + Point(shape.size_x, 0),
              shape.ll + Point(shape.size_x, shape.size_y),
              shape.ll + Point(0, shape.size_y)]
    geo_code(ClosedPolygon(points), start, tag::Int)
end

"""Convert the polygon to code in gmsh scripting language"""
function geo_code(shape::ClosedPolygon, start::Int=1, tag::Int=0)
    # NOTE: start is a counter for gmsh entities. Every shape we convert is a closed
    # loop and we might want to label the created loop(surface) as well as its boundary

    # Point part
    points = shape.points
    npoints = length(points)
    pindices = start:(start+npoints-1)

    points = ["Point($(index)) = {$(point.x), $(point.y), 0., SIZE};"
               for (index, point) in zip(pindices, points)]
    points = join(points, "\n")

    # Line part
    start = last(pindices) + 1
    lines = [zip(pindices[2:end], pindices[1:end-1])..., (first(pindices), last(pindices))]
    nlines = length(lines)
    lindices = start:(start+nlines-1)

    lines = ["Line($(index)) = {$(first(line)), $(last(line))};"
               for (index, line) in zip(lindices, lines)]
    lines = join(lines, "\n")

    # The loop
    start = last(lindices) + 1
    loop = "Line Loop($(start)) = {$(join(map(x-> "$(x)", lindices), ", "))};"

    # Physical groups
    physical_surface = "Physical Surface($(tag)) = {$(start)};"
    physical_line = "Physical Line($(tag)) = {$(join(map(x-> "$(x)", lindices), ", "))};"

    code = join([points, lines, loop, physical_surface, physical_line], "\n")
    start += 1  # The next index suitable for numbering

    (code, start)
end

# FIXME: do we want cotrol over sizes?
# FIXME: Obligue
# FIXME: Canvas
