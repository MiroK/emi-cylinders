include("emi.jl")
using EMI
using EMI.Draw, EMI.Gmsh

# s = Square(Point(0, 0), 1)
# s1 = Square(Point(2, 0), 1)
# s2 = Square(Point(4, 0), 1)
# 
# c = Canvas()
# c = c + s
# c = c + s1
# c = c + s2
# set_bbox!(c, 0.2, 0.3)

#canvas = Canvas()
#for N in 3:10
#    ngon = NGon(N, 1) + (N-3)*Point(2, 0)
#
#    # println(ngon)
#    canvas = canvas + ngon 
#end

canvas = Canvas()
# canvas = canvas + Circle(Point(0, 0), 1)
canvas = canvas + Ellipse(Point(0, 0), 1, 0.5)

set_bbox!(canvas, 0.2, 0.3)
println(gmsh_script(canvas))

#println(code)
#println(start)
