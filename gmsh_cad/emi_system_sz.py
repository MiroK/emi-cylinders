import petsc4py, sys

# I like to be in control of the command line
petsc4py.init(sys.argv)

from dolfin import *
from petsc4py import PETSc
from mpi4py import MPI

comm = MPI.COMM_WORLD

# Create stages for logging
st_mesh = PETSc.Log.Stage("Mesh SetUp")
st_ass  = PETSc.Log.Stage("Assembly")
st_pc   = PETSc.Log.Stage("BDDC SetUp")
st_ksp  = PETSc.Log.Stage("KSP Solve")

# Dolfin global parameters
opts = PETSc.Options()
parameters['mesh_partitioner'] = opts.getString('meshpartitioner','ParMETIS')
parameters['form_compiler']['representation'] = 'uflacs'
parameters['form_compiler']['cpp_optimize'] = True
parameters['form_compiler']['cpp_optimize_flags'] = '-O3 -ffast-math -march=native'
parameters['ghost_mode'] = 'shared_facet'

mesh_file = opts.getString('meshfile','Tiles/tile_1_narrow_2_2.h5')
#mesh_file = '2Dtest/cell_grid_2d.h5'
#mesh_file = '500Kdofs/cell_grid.h5'
#mesh_file = '8Mdofs/cell_grid.h5'
#mesh_file = 'cell_grid_2d.h5'

with st_mesh :
  h5 = HDF5File(comm, mesh_file, 'r')
  mesh = Mesh()
  h5.read(mesh, 'mesh', False)
  # The mesh comes in micro meters. Below it is more convenient to work in cm
  mesh.coordinates()[:] *= 1E-4

  # Facets in the mesh have tags 0, 1. One is for interfaces between
  # cells and cells and the exterior. The domain is split into 2 subdomains
  # marked as 1 and 0 (cell interior, cell exterior). These differ by conductivities

  ext_tag, int_tag = 0, 1
  not_iface_tag, iface_tag = 0, 1

  surfaces = MeshFunction('size_t', mesh, mesh.topology().dim()-1, 0)
  h5.read(surfaces, 'surfaces')

  volumes = MeshFunction('size_t', mesh, mesh.topology().dim(), 0)
  h5.read(volumes, 'volumes')

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

# W.sub(0) bcs set strongly correspond to insulation. Skipping (tau.n)*u on
# bdry means potential there is weakly zero, grounding.

if opts.getBool('view_volumes',False) :
  file = File("Volumes.pvd")
  file << volumes

if opts.getBool('view_surfaces',False) :
  file = File("Surfaces.pvd")
  file << surfaces

# Make measures aware of subdomains
dx = Measure('dx', domain=mesh, subdomain_data=volumes)
dS = Measure('dS', domain=mesh, subdomain_data=surfaces)
ds = Measure('ds', domain=mesh, subdomain_data=surfaces)

# Normal fo the INTERIOR surface. Note that 1, 2 marking of volume makes
# 2 cells the '+' cells w.r.t to surface and n('+') would therefore be their
# outer normal (that is an outer normal of the outside). ('-') makes the orientation
# right (set later)
n = FacetNormal(mesh)

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
a = ((1/cond_int)*inner(sigma, tau)*dx(int_tag)+(1/cond_ext)*inner(sigma, tau)*dx(ext_tag)
     - inner(div(tau), u)*dx(int_tag) - inner(div(tau), u)*dx(ext_tag)
     + inner(p('+'), dot(tau('+'), n('-')))*dS(iface_tag)
     - inner(div(sigma), v)*dx(int_tag) - inner(div(sigma), v)*dx(ext_tag)
     + inner(q('+'), dot(sigma('+'), n('-')))*dS(iface_tag)
     - (C_m/dt_fem)*inner(q('+'), p('+'))*dS(iface_tag))

L = inner(q('+'), I_ion('+')-(C_m/dt_fem)*p0('+'))*dS(iface_tag)

# NOTE: in general the cell interfaces might end up on the surface boundary
a += (inner(p, dot(tau, n))*ds(iface_tag)
      + inner(q, dot(sigma, n))*ds(iface_tag)
      - (C_m/dt_fem)*inner(q, p)*ds(iface_tag))

L += inner(q, I_ion-(C_m/dt_fem)*p0)*ds(iface_tag)

# Additional terms to set to zero the dofs of W.sub(2) which are not on
# the interfaces
a -= inner(p('+'), q('+'))*dS(not_iface_tag) + inner(p, q)*ds(not_iface_tag)
L -= inner(Constant(0)('+'), q('+'))*dS(not_iface_tag) + inner(Constant(0), q)*ds(not_iface_tag)

# In the presence of internal facet terms, the local to global maps have ghost cells (through shared facets)
# However, only one process insert values there -> we need to prune empty local rows/columns
# from the other processes. This can be done programmatically from PETSc, but it is currently
# in a development branch
opts.setValue("-test_matis_fixempty", None)

# Assemble system: matrix in unassembled (MATIS) format
A, b = PETScMatrix(), PETScVector()
as_backend_type(A).mat().setOptionsPrefix("test_")
as_backend_type(A).mat().setType("is")
as_backend_type(A).mat().setFromOptions()
#assemble_system(a, L, bcs, A_tensor=A, b_tensor=b)
with st_ass : assemble_system(a, L, A_tensor=A, b_tensor=b)

## test unassembled format
#Aij, bij = PETScMatrix(), PETScVector()
#assemble_system(a, L, A_tensor=Aij, b_tensor=bij)
#A1 = PETSc.Mat()
#PETSc.Sys.Print("CONVERTING")
#as_backend_type(A).mat().convert('aij',A1)
#A1.axpy(-1,as_backend_type(Aij).mat())
#PETSc.Sys.Print('Norm assembly: %f' % A1.norm(PETSc.NormType.INFINITY))
#del A1, Aij

## test reassembly
#A1 = as_backend_type(A).mat().duplicate(copy=True)
#as_backend_type(A).mat().zeroEntries()
#assemble_system(a, L, A_tensor=A, b_tensor=b)
#
#A1.axpy(-1,as_backend_type(A).mat())
#PETSc.Sys.Print('Norm reassembly: %f' % A1.norm(PETSc.NormType.INFINITY))
#A1.viewFromOptions("-diff_view")

as_backend_type(A).mat().viewFromOptions("-my_view")

# Create PETSc Krylov solver (from petsc4py)
ksp = PETSc.KSP()
ksp.create(PETSc.COMM_WORLD)
ksp.setOptionsPrefix("test_")

# Set the Krylov solver type and set tolerances
# We can use CG since DRT dofs will never be part of the interface between subdomains
# This is because DRT dofs are only coupled to RT dofs and not cell dofs
# Note: If DRT dofs were coupled to DG dofs, we could have put all the DRT dofs (see setBDDCPrimalVerticesIS) that lie on the
#       interface between subdomains in the primal space
#       The coarse matrix will be diagonal on those DRT dofs (most of them) that are not on
#       the intra/extra cellular interface
ksp.setType("cg")
ksp.setTolerances(rtol=1.0e-8, atol=1.0e-12, divtol=1.0e10, max_it=300)
ksp.setOperators(as_backend_type(A).mat())
pc = ksp.getPC()
pc.setType("bddc")
#pc.setBDDCPrimalVerticesIS(PETSc.IS().createGeneral(W.sub(2).dofmap().dofs(),PETSc.COMM_WORLD))

# Options
opts.setValue("-test_ksp_view", None)
opts.setValue("-test_ksp_converged_reason", None)
opts.setValue("-test_ksp_monitor_singular_value", None)
opts.setValue("-test_ksp_norm_type", "natural")

# Don't turn these off
opts.setValue("-test_pc_bddc_detect_disconnected", None)
opts.setValue("-test_pc_bddc_detect_disconnected_filter", None)
opts.setValue("-test_pc_bddc_use_faces", None)
opts.setValue("-test_pc_bddc_benign_trick", None)
opts.setValue("-test_pc_bddc_nonetflux", None)
opts.setValue("-test_pc_bddc_schur_exact", None)
opts.setValue("-test_pc_bddc_use_deluxe_scaling", None)
opts.setValue("-test_pc_bddc_deluxe_zerorows", None)
opts.setValue("-test_pc_bddc_use_local_mat_graph", "0")
opts.setValue("-test_pc_bddc_adaptive_userdefined", None)

# Better off you have MUMPS or SUITESPARSE for the factorizations

# Sometimes MUMPS fails with error -9 (increase Schur workspace....this is annoying)
opts.setValue("-test_sub_schurs_mat_mumps_icntl_14",500)
opts.setValue("-mat_mumps_icntl_14",500)

# Local solvers (MUMPS)
opts.setValue("-test_pc_bddc_dirichlet_pc_type","cholesky") # This is actually LDL^T
opts.setValue("-test_pc_bddc_neumann_pc_type","cholesky") # This is actually LDL^T
opts.setValue("-test_pc_bddc_dirichlet_pc_factor_mat_solver_type","mumps")
opts.setValue("-test_pc_bddc_neumann_pc_factor_mat_solver_type","mumps")

# Alternative local factorizations with SUITESPARSE
#opts.setValue("-test_pc_bddc_dirichlet_pc_type","lu")
#opts.setValue("-test_pc_bddc_neumann_pc_type","lu")
#opts.setValue("-test_pc_bddc_dirichlet_pc_factor_mat_solver_type","umfpack")
#opts.setValue("-test_pc_bddc_neumann_pc_factor_mat_solver_type","umfpack")

# Alternative redundant coarse factorization with SUITESPARSE
#opts.setValue("-test_pc_bddc_coarse_pc_type","redundant")
#opts.setValue("-test_pc_bddc_coarse_redundant_pc_factor_mat_solver_type","umfpack")

# Number of additional levels : 0 means standard 2-level BDDC
nlevels = opts.getInt('nlevels',0)
coarsening = opts.getInt('coarsening',2)

# debugs check from level
lcheck = opts.getInt('check_from_level',nlevels+1)
dcheck = opts.getInt('check_from_level_dbg',1)

# Coarse solver (MUMPS or BDDC)
opts.setValue("-test_pc_bddc_coarse_pc_factor_mat_solver_type","mumps")
if nlevels < 1:
  opts.setValue("-test_pc_bddc_coarse_pc_type","cholesky") # This is actually LDL^T

if lcheck <= 0 :
  opts.setValue("-test_pc_bddc_check_level", dcheck)

if lcheck <= 1 :
  opts.setValue("-test_pc_bddc_coarse_pc_bddc_check_level", dcheck)

opts.setValue("-test_pc_bddc_levels",nlevels)
opts.setValue("-test_pc_bddc_coarsening_ratio",coarsening)
opts.setValue("-test_pc_bddc_coarse_sub_schurs_mat_mumps_icntl_14",500)
opts.setValue("-test_pc_bddc_coarse_pc_bddc_detect_disconnected", None)
opts.setValue("-test_pc_bddc_coarse_pc_bddc_detect_disconnected_filter", None)
opts.setValue("-test_pc_bddc_coarse_pc_bddc_use_deluxe_scaling",None)
opts.setValue("-test_pc_bddc_coarse_pc_bddc_deluxe_zerorows",None)
opts.setValue("-test_pc_bddc_coarse_pc_bddc_schur_exact",None)
opts.setValue("-test_pc_bddc_coarse_pc_bddc_use_local_mat_graph",None)
opts.setValue("-test_pc_bddc_coarse_check_ksp_monitor",None)
opts.setValue("-test_pc_bddc_coarse_check_ksp_converged_reason",None)
opts.setValue("-test_pc_bddc_coarse_check_ksp_type","cg")
opts.setValue("-test_pc_bddc_coarse_check_ksp_norm_type","natural")
if opts.getBool('cheby',False) :
  opts.setValue("-test_pc_bddc_coarse_ksp_type","chebyshev")
  opts.setValue("-test_pc_bddc_use_coarse_estimates",None)
  opts.setValue("-test_pc_bddc_coarse_pc_bddc_use_coarse_estimates",None)

for j in range(1, nlevels):
  if lcheck < j+2 :
    opts.setValue("-test_pc_bddc_coarse_l" + str(j) + "_pc_bddc_check_level", dcheck)

  opts.setValue("-test_pc_bddc_coarse_l" + str(j) + "_sub_schurs_mat_mumps_icntl_14",500)
  opts.setValue("-test_pc_bddc_coarse_l" + str(j) + "_pc_bddc_detect_disconnected",None)
  opts.setValue("-test_pc_bddc_coarse_l" + str(j) + "_pc_bddc_detect_disconnected_filter",None)
  opts.setValue("-test_pc_bddc_coarse_l" + str(j) + "_pc_bddc_use_deluxe_scaling",None)
  opts.setValue("-test_pc_bddc_coarse_l" + str(j) + "_pc_bddc_deluxe_zerorows",None)
  opts.setValue("-test_pc_bddc_coarse_l" + str(j) + "_pc_bddc_schur_exact",None)
  opts.setValue("-test_pc_bddc_coarse_l" + str(j) + "_pc_bddc_use_local_mat_graph",None)
  opts.setValue("-test_pc_bddc_coarse_l" + str(j) + "_check_ksp_type","cg")
  opts.setValue("-test_pc_bddc_coarse_l" + str(j) + "_check_ksp_norm_type","natural")
  opts.setValue("-test_pc_bddc_coarse_l" + str(j) + "_check_ksp_monitor",None)
  opts.setValue("-test_pc_bddc_coarse_l" + str(j) + "_check_ksp_converged_reason",None)
  if opts.getBool('cheby',False) :
    opts.setValue("-test_pc_bddc_coarse_l" + str(j) + "_ksp_type","chebyshev")
    opts.setValue("-test_pc_bddc_coarse_l" + str(j) + "_pc_bddc_use_coarse_estimates",None);

# prevent from having a bad coarse solver
opts.setValue("-test_pc_bddc_coarse_l3_redundant_pc_factor_mat_solver_type","mumps");
opts.setValue("-test_pc_bddc_coarse_l3_redundant_pc_type","cholesky")
opts.setValue("-test_pc_bddc_coarse_l2_redundant_pc_factor_mat_solver_type","mumps");
opts.setValue("-test_pc_bddc_coarse_l2_redundant_pc_type","cholesky")
opts.setValue("-test_pc_bddc_coarse_l1_redundant_pc_factor_mat_solver_type","mumps");
opts.setValue("-test_pc_bddc_coarse_l1_redundant_pc_type","cholesky")
opts.setValue("-test_pc_bddc_coarse_redundant_pc_factor_mat_solver_type","mumps")
opts.setValue("-test_pc_bddc_coarse_redundant_pc_type","cholesky")

#Solve
ksp.setFromOptions()
with st_pc : pc.setUp()

sol = Function(W)

with st_ksp : ksp.solve(as_backend_type(b).vec(), as_backend_type(sol.vector()).vec())

# prevent from deadlocks when garbage collecting
del pc, ksp

# Visualize computed solution
if opts.getBool('view_solution',False) :
  as_backend_type(sol.vector()).update_ghost_values()
  (sols, solv, solq) = sol.split()
  file = File("electric.pvd")
  file << sols
  file = File("potential.pvd")
  file << solv
  #file = File("tpotential.pvd")
  #file << solq
