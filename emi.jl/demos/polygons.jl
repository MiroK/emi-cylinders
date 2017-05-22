include("../emi.jl")
using EMI
using EMI.Draw, EMI.Gmsh

canvas = Canvas()
for N in 3:10
    ngon = NGon(N, 1) + (N-3)*Point(2, 0)  # Fixed size shit
    canvas = canvas + ngon 
end

set_bbox!(canvas, 0.2, 0.3)
gmsh_script(canvas, 0.2)
