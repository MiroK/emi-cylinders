include("emi.jl")
using EMI
using EMI.Draw, EMI.Gmsh

s = Square(Point(0, 0), 1)
s1 = Square(Point(2, 0), 1)
s2 = Square(Point(4, 0), 1)

c = Canvas()
c = c + s
c = c + s1
c = c + s2
set_bbox!(c, 0.2, 0.3)

println(gmsh_script(c))

#println(code)
#println(start)
