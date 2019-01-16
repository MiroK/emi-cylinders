"""
This example solves the monodomain model with the parsimonious cell model,
for the thin sheet geometry we will use in the EMI/bidomain comparison.
"""



# Import the cbcbeat module
from cbcbeat import *
from parsimonious import *

import numpy as np

import matplotlib.pyplot as plt

# Turn on FFC/FEniCS optimizations
parameters["form_compiler"]["representation"] = "uflacs"
parameters["form_compiler"]["cpp_optimize"] = True
flags = ["-O3", "-ffast-math", "-march=native"]
parameters["form_compiler"]["cpp_optimize_flags"] = " ".join(flags)
parameters["form_compiler"]["quadrature_degree"] = 3


# Define the computational domain
"""
Mesh dimensions 20mm*4.4mm*0.022mm
"""

tilesize = (0.100,0.025,0.025) #mm**3

#create a mesh of 200x200x1 cells
L = 200*tilesize[0]
W = 200*tilesize[1]
H = tilesize[2]
dx = 0.1 #mm
n_L = int(L/dx)
n_W = int(W/dx)
n_H = 1
mesh = BoxMesh(Point(0,0,0),Point(L,W,H),n_L,n_W,n_H)
print("Number of nodes: %g, number of elements: %g" %(mesh.num_cells(), mesh.num_vertices()))


time = Constant(0.0)

# scaled conductivities:
# g = g_0/(C_m*chi), with M_0 being the physical cond.
C_m = 0.01 #microF*mm**-2
chi = 200 #mm**-1

g_il = 0.2/(C_m*chi)
g_it = 0.02/(C_m*chi)
g_el = 0.8/(C_m*chi)
g_et = 0.2/(C_m*chi)

M_e = Constant((((g_el, 0, 0),
                 (0, g_et, 0),
                 (0, 0, g_et))))

M_i = Constant((((g_il, 0, 0),
                 (0, g_it, 0),
                 (0, 0, g_it))))


# define cell model with zero stimulus
cell_params = Parsimonious.default_parameters()
cell_params['amp'] = 0.0
cell_model = Parsimonious(params=cell_params)

# Define a stimulus current
stim = Expression("t > 5.0 && t < 5.5 && x[0] < length*0.1? 80:0", t=time, length=L, degree=1)

#Electrodes placed symmetrically around x[0] = 10, 4mm apart, 1mm long
e1 = CompiledSubDomain("near(x[1], 0) && x[0] > 7-tol && x[0] < 8+tol && on_boundary",tol = 1e-10)
e2 = CompiledSubDomain("near(x[1], 0) && x[0] > 12-tol && x[0] < 13+tol && on_boundary",tol = 1e-10)

#all currents are scaled with Cm*chi, J_app given in uA*mm**-2:
s1 = Expression("t < 2.0 ? J_s1*C_m*chi:0", t=time,
                J_s1=28, C_m=C_m, chi = chi, degree=1)

s2 = Expression("t > t_s2 && t< t_s2 + 2.0 ? J_s2*C_m*chi:0", t=time,
                t_s2 = 190, J_s2=28, C_m=C_m, chi = chi, degree=1)


boundary_markers = MeshFunction('size_t',mesh,mesh.topology().dim()-1)
boundary_markers.set_all(0)
e1.mark(boundary_markers,1)
e2.mark(boundary_markers,2)
applied_c = Markerwise((-(s1+s2),(s1+s2)),(1,2),boundary_markers)


# Collect this information into the CardiacModel class
cardiac_model = CardiacModel(mesh, time, M_i, M_e, cell_model, stimulus=stim,applied_current=None)

# Customize and create a splitting solver
ps = SplittingSolver.default_parameters()
ps["theta"] = 0.5                        # Second order splitting scheme
"""
ps["pde_solver"] = "monodomain"          # Use Monodomain model for the PDEs
ps["CardiacODESolver"]["scheme"] = "RL1" # 1st order Rush-Larsen for the ODEs
ps["MonodomainSolver"]["linear_solver_type"] = "iterative"
ps["MonodomainSolver"]["algorithm"] = "cg"
ps["MonodomainSolver"]["preconditioner"] = "petsc_amg"
"""
ps["pde_solver"] = "bidomain"          # Use Monodomain model for the PDEs
ps["CardiacODESolver"]["scheme"] = "RL1" # 1st order Rush-Larsen for the ODEs
ps["BidomainSolver"]["linear_solver_type"] = "iterative"
ps["BidomainSolver"]["algorithm"] = "cg"
ps["BidomainSolver"]["preconditioner"] = "petsc_amg"

solver = SplittingSolver(cardiac_model, params=ps)

# Extract the solution fields and set the initial conditions
(vs_, vs, vur) = solver.solution_fields()
vs_.assign(cell_model.initial_conditions())

# Time stepping parameters
dt = 0.1
T = 5.0
interval = (0.0, T)


"""
Plot time series in a set of points, save full solution 
to file:
"""
no_points = 5
x = [np.array([x,W/2,H/2]) for x in np.linspace(0,L,no_points)]
time = [interval[0]]
v_series = [[vs_(x[i])[0] for i in range(no_points)]]
u_series = [[0.0 for i in range(no_points)]]

v_file = File('solutions/membrane.pvd')
u_file = File('solutions/extra.pvd')


timer = Timer("XXX Forward solve") # Time the total solve
for (timestep, fields) in solver.solve(interval, dt):
    print("(t_0, t_1) = (%g, %g)" %timestep)

    # Extract the components of the field (vs_ at previous timestep,
    # current vs, current vur)
    (vs_, vs, vur) = fields

    time.append(timestep[1])
    v_series.append([vs(x[i])[0] for i in range(no_points)])
    u_series.append([vur(x[i])[1] for i in range(no_points)])

    #plot every 5th millisecond:
    if (timestep[-1]+0.5*dt)%5 < dt:
        v_file << vs.split()[0]
        u_file << vur.split()[1]
    # Print memory usage (just for the fun of it)
    #print memory_usage()

timer.stop()
list_timings(TimingClear.keep, [TimingType.user])

"""
code below is for plotting time series of v and/or u, in a 
a series of points defined above.
"""
outfile = open('time_series.txt','w')
outfile.write(str(time))
outfile.write('\n')

v_values = list(zip(*v_series))
u_values = list(zip(*u_series))
for i in range(no_points):
    plt.plot(time, v_values[i],label = 'v, x = %g' %(x[i][0]))
    outfile.write(str(v_values[i]))
    outfile.write('\n')
    #plt.plot(time,u_values[i],label = 'u, x= %g' %(x[i][0]) )

outfile.close()

plt.legend()
plt.xlabel('Time (ms)')
plt.ylabel('Membrane potential (mV)')
plt.savefig('time_series.png')
plt.show()

print("Success!")

