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

canvas = Canvas()
#for N in 3:10
#    ngon = NGon(N, 1) + (N-3)*Point(2, 0)
#
#    # println(ngon)
#    canvas = canvas + ngon 
#end
loop = Loop([Line(Point(0, 0), Point(1, 0)),
             Line(Point(1, 0), Point(1, 1)),
             Line(Point(1, 1), Point(0, 0))])


#canvas = canvas + Circle(Point(0, 0), 1)
#canvas = canvas + Ellipse(Point(4, 0), 1, 0.5)
canvas = canvas + loop

set_bbox!(canvas, 0.2, 0.3)
println(gmsh_script(canvas))

#println(code)
#println(start)
