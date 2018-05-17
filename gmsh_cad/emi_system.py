from dolfin import *

parameters['form_compiler']['representation'] = 'uflacs'
parameters['form_compiler']['cpp_optimize'] = True
parameters['form_compiler']['cpp_optimize_flags'] = '-O3 -ffast-math -march=native'
parameters['ghost_mode'] = 'shared_facet'

mesh_file = 'cell_grid_2d.h5'

comm = mpi_comm_world()
h5 = HDF5File(comm, mesh_file, 'r')
mesh = Mesh()
h5.read(mesh, 'mesh', False)
# The mesh comes in micro meters. Below it is more convenient to work in cm
mesh.coordinates()[:] *= 1E-4
# Facets in the mesh have tags 0, 1, 2. One is for interfaces between
# cells and cells and the exterior. Two is used for marking boundary facets
# of the domain - this is where typically zero DirichletBCs are applied
# for the potential
surfaces = MeshFunction('size_t', mesh, mesh.topology().dim()-1)
h5.read(surfaces, 'facet')
# The domain is split into 2 subdomains marked as 1 and 2 (cell interior,
# cell exterior). These differ by conductivities
volumes = MeshFunction('size_t', mesh, mesh.topology().dim())
h5.read(volumes, 'physical')

cell = mesh.ufl_cell()
# We have 3 spaces S for sigma = -kappa*grad(u)   [~electric field]
#                  U for potential u
#                  Q for transmebrane potential p
Sel = FiniteElement('RT', cell, 1)
Vel = FiniteElement('DG', cell, 0)
Qel = FiniteElement('Discontinuous Lagrange Trace', cell, 0)

W = FunctionSpace(mesh, MixedElement([Sel, Vel, Qel]))
sigma, u, p = TrialFunctions(W)
tau, v, q = TestFunctions(W)

# Grounding for potential
bcs = [DirichletBC(W.sub(2), Constant(0), surfaces, 2)]

# Make measures aware of subdomains
dx = Measure('dx', domain=mesh, subdomain_data=volumes)
dS = Measure('dS', domain=mesh, subdomain_data=surfaces)
ds = Measure('ds', domain=mesh, subdomain_data=surfaces)
    
# Normal fo the INTERIOR surface. Note that 1, 2 marking of volume makes
# 2 cells the '+' cells w.r.t to surface and n('+') would therefore be their
# outer normal (that is an outer normal of the outside). ('-') makes the orientation
# right
n = FacetNormal(mesh)('-')

# Now onto the weak form
# Electric properties of membrane and interior/exterior
C_m = Constant(1)         # 1 mu F / cm^2                        
cond_int = Constant(5)    # 5 mS / cm
cond_ext = Constant(20)   # 20 mS / cm
# Time step
dt_fem = Constant(1E-3)   # ms

# The source term as a function Q is coming from ODE solver. Here it is
# just some random function
Q = FunctionSpace(mesh, Qel)  
p0 = interpolate(Constant(1), Q)
# And additional source on the boundary is the ionic current. For simplicity
I_ion = p0

# The system
a = ((1/cond_int)*inner(sigma, tau)*dx(1)+(1/cond_ext)*inner(sigma, tau)*dx(2)
     - inner(div(tau), u)*dx(1) - inner(div(tau), u)*dx(2)
     + inner(p('+'), dot(tau('+'), n))*dS(1)
     - inner(div(sigma), v)*dx(1) - inner(div(sigma), v)*dx(2)
     + inner(q('+'), dot(sigma('+'), n))*dS(1)
     - (C_m/dt_fem)*inner(q('+'), p('+'))*dS(1))

L = inner(q('+'), I_ion('+')-(C_m/dt_fem)*p0('+'))*dS(1)

# Additional terms to set to zero the dofs of W.sub(2) which are not on
# the interfaces
a -= inner(p('+'), q('+'))*dS(0) + inner(p, q)*ds(2)
L -= inner(Constant(0)('+'), q('+'))*dS(0) + inner(Constant(0), q)*ds(2)

A, b = PETScMatrix(), PETScVector()
assemble_system(a, L, bcs, A_tensor=A, b_tensor=b)

#import numpy as np
#for i in range(A.size(0)):
#    cols, vals = A.getrow(i)
#    assert np.linalg.norm(vals, 1) > 0
print A.size(0)
