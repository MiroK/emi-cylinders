include("emi.jl")
using EMI
using EMI.Draw, EMI.Gmsh

s = Square(Point(0, 0), 1)
code, start = geo_code(s, 1)

println(code)
println(start)
