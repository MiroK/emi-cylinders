import petsc4py, sys
petsc4py.init(sys.argv)
from petsc4py import PETSc

from cbcbeat.cellsolver import CardiacODESolver

from parsimonious import Parsimonious
from emi_pde_tune import pde_components
from probing import ApproxPointProbe, probe_cell_at

from mpi4py import MPI as pyMPI
from dolfin import *
import numpy as np


# Setup 
mesh_file = './tile_1_hein_GMSH307_10_1.h5'
# Get mesh to setup the ode solver
comm = MPI.comm_world  # FIXME!
h5 = HDF5File(comm, mesh_file, 'r')
mesh = Mesh()
h5.read(mesh, 'mesh', False)

# Mesh box dimensions in mm. So mesh is expected to be nn as well.
# NOTE: update by hand - waste of resouces to do it with mpi
ncells_x, ncells_y = 10, 1
Length = 0.1*ncell_x    # In x as multiples of single cell
Width = 0.025*ncells_y   # In y
Height = 0.025

# Here are some sample points for transmembrane potential: On top surface
# of the cylinder
sample_points = np.array([probe_cell_at(m, n=0, level=2, point=p)
                         for m in range(ncell_x) for p in (-1, 1)])

print(sample_points)

exit()

# ODE setup
C_m = 0.01 #microF*mm**-2
chi = 200 #mm**-1

# define cell model with zero stimulus
cell_params = Parsimonious.default_parameters()
cell_params['amp'] = 0.0
cell_model = Parsimonious(params=cell_params)

# Define a stimulus current
stim = Expression("t > 0.0 && t < 0.5 && x[0] < 0.2 ? 200:0", t=0., length=Length, degree=1)

ode_parameters = {'dt': 1E-2,
                  'stimulus': stim}
# ---
ode_solver = CardiacODESolver(mesh,
                              time=Constant(0),
                              model=cell_model,
                              I_s=ode_parameters['stimulus'])

interval = (0.0, 5.0)


# NOTE: a generator; nothing is computed so far
ode_solutions = ode_solver.solve(interval, ode_parameters['dt'])  # Potentials only

# EMI setup
emi_parameters = {'C_int': 3*C_m,
                  'C_ext': 2*C_m,
                  'C_m': C_m,
                  'dt': 0.1,
                  # Potentials as P1 function on the mesh are wired up with the EMI = p0
                  'p0': ode_solver.vs[0],
                  'I_ion': Constant(0.0)}  # FIXME: what should this be?

# Pieces for emi are a, L, W, bcs (so assembly) and transfer operator
emi_pieces = pde_components(mesh, h5, emi_parameters)

# Since these are meant to be simulations with fixed time step, the matrix
# is assembled once
a, L, bcs = (emi_pieces[key] for key in ('a', 'L', 'bcs'))

A, b = PETScMatrix(), PETScVector()

# --------------------------------------------------------------------
# Setup Stefano's solver
# --------------------------------------------------------------------
opts = PETSc.Options()

# In the presence of internal facet terms, the local to global maps have ghost cells (through shared facets)
# However, only one process insert values there -> we need to prune empty local rows/columns
# from the other processes. This can be done programmatically from PETSc, but it is currently
# in a development branch
opts.setValue("-test_matis_fixempty", None)

Amat = A.mat()
Amat.setOptionsPrefix("test_")
Amat.setType("is")
Amat.setFromOptions()

# Assembly
emi_assembler = SystemAssembler(a, L, bcs)
emi_assembler.assemble(A)

Amat.viewFromOptions("-my_view")

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
ksp.setOperators(Amat)
pc = ksp.getPC()
pc.setType("bddc")

# Options
opts.setValue("-test_ksp_view", None)
opts.setValue("-test_ksp_converged_reason", None)
opts.setValue("-test_ksp_monitor_singular_value", None)
opts.setValue("-test_ksp_norm_type", "natural")

# Don't turn these off
opts.setValue("-test_pc_bddc_detect_disconnected", None)
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

nlevels = 1
# Coarse solver (MUMPS or BDDC)
opts.setValue("-test_pc_bddc_coarse_pc_factor_mat_solver_type","mumps")
if nlevels < 1:
    opts.setValue("-test_pc_bddc_coarse_pc_type","cholesky") # This is actually LDL^T

opts.setValue("-test_pc_bddc_levels",nlevels)
opts.setValue("-test_pc_bddc_coarse_sub_schurs_mat_mumps_icntl_14",500)
opts.setValue("-test_pc_bddc_coarse_pc_bddc_use_deluxe_scaling",None)
opts.setValue("-test_pc_bddc_coarse_pc_bddc_deluxe_zerorows",None)
opts.setValue("-test_pc_bddc_coarse_pc_bddc_schur_exact",None)
opts.setValue("-test_pc_bddc_coarse_pc_bddc_use_local_mat_graph",None)
opts.setValue("-test_pc_bddc_coarse_check_ksp_monitor",None)
opts.setValue("-test_pc_bddc_coarse_check_ksp_converged_reason",None)
opts.setValue("-test_pc_bddc_coarse_check_ksp_type","cg")
opts.setValue("-test_pc_bddc_coarse_check_ksp_norm_type","natural")

for j in range(nlevels):
    opts.setValue("-test_pc_bddc_coarse_l" + str(j) + "_sub_schurs_mat_mumps_icntl_14",500)
    opts.setValue("-test_pc_bddc_coarse_l" + str(j) + "_pc_bddc_use_deluxe_scaling",None)
    opts.setValue("-test_pc_bddc_coarse_l" + str(j) + "_pc_bddc_deluxe_zerorows",None)
    opts.setValue("-test_pc_bddc_coarse_l" + str(j) + "_pc_bddc_schur_exact",None)
    opts.setValue("-test_pc_bddc_coarse_l" + str(j) + "_pc_bddc_use_local_mat_graph",None)
    opts.setValue("-test_pc_bddc_coarse_l" + str(j) + "_check_ksp_type","cg")
    opts.setValue("-test_pc_bddc_coarse_l" + str(j) + "_check_ksp_norm_type","natural")
    opts.setValue("-test_pc_bddc_coarse_l" + str(j) + "_ksp_type","chebyshev")
    opts.setValue("-test_pc_bddc_coarse_l" + str(j) + "_check_ksp_monitor",None)
    opts.setValue("-test_pc_bddc_coarse_l" + str(j) + "_check_ksp_converged_reason",None)
    opts.setValue("-test_pc_bddc_coarse_l" + str(j) + "_check_ksp_type","cg")
    opts.setValue("-test_pc_bddc_coarse_l" + str(j) + "_check_ksp_norm_type","natural")

# prevent to have a bad coarse solver
opts.setValue("-test_pc_bddc_coarse_l3_redundant_pc_factor_mat_solver_package","mumps");
opts.setValue("-test_pc_bddc_coarse_l2_redundant_pc_factor_mat_solver_package","mumps");
opts.setValue("-test_pc_bddc_coarse_l1_redundant_pc_factor_mat_solver_package","mumps");
opts.setValue("-test_pc_bddc_coarse_redundant_pc_factor_mat_solver_package","mumps")

ksp.setFromOptions()

# Solution loop consists of ODE solves in interval and every emi_dt/ode_dt
# there is an exchange between transmembrane potentials
fem_ode_sync = int(emi_parameters['dt']/ode_parameters['dt'])
# The EMI is solver for electric field, potential, transm. potential
W = emi_pieces['W']
wh = Function(W)  
# The transmembrane potentials in ode solver are in
P1 = ode_solver.VS.sub(0).collapse() 
p_ode = Function(P1)
# Meanwhile transmembrane potentials from DLT are in
P0 = W.sub(2).collapse()
p_emi = Function(P0)

# We will the transfer operator W(2)-> P0 -> P1 -> ODE(0)
# W(2)-> P0
to_P0_from_EMI = FunctionAssigner(P0, W.sub(2))
# P0 -> P1
P0_to_P1 = emi_pieces['toODE']
# P1 -> ODE(0)
to_ODE_from_P1 = FunctionAssigner(ode_solver.VS.sub(0), P1)


# solver = LUSolver()
# solver.set_operator(A)
# FIXME: how to get it then, config via PETSc?
#solver.parameters['reuse_factorization'] = True

# For now the output consists of point values along x through domain center
no_points = 5
probe_points = np.c_[np.linspace(0, Length, no_points),
                     0.5*Width*np.ones(no_points),
                     0.5*Height*np.ones(no_points)]
# NOTE: it it is not possible to sample wh because of DLT elements. So
# we need one more assigner for the potential
V = W.sub(1).collapse()
potential_assigner = FunctionAssigner(V, W.sub(1))

uh = Function(V)
potential_assigner.assign(uh, wh.sub(1))

# probes = PointProbe(uh, probe_points)
# Initial reading
# probe_values = probes.sample(uh)

# The idea is to have columns of t and potentials readings
# if pyMPI.COMM_WORLD.rank == 0: table = np.r_[interval[0], probe_values[:, 0]]

# PETS.Vec s for solver
b_vec = b.vec()
x_vec = as_backend_type(wh.vector()).vec()
    
step_count = 0
for ((t0, t1), ode_solution) in ode_solutions:
    step_count += 1

    if step_count == fem_ode_sync:
        step_count = 0
        # L is wired up with ode so that is updated. Optionally update
        # values of boundary conditions
        for expr in electrode_values:
            hasattr(expr, 't') and setattr(expr, 't', float(t1))

        # The new rhs
        emi_assembler.assemble(b)
      
        # New (sigma, u, p) ...
        info('\tSolving linear system of size %d' % A.size(0))
        
        ksp.solve(b_vec, x_vec)

        # Update emi potential to standalone DLT function. To emi from wh(2)
        to_P0_from_EMI.assign(p_emi, wh.sub(2))
        # Map DLT function to standalone P1 function (potential for ODE)
        P0_to_P1.mult(p_emi.vector(), p_ode.vector())
        as_backend_type(p_ode.vector()).update_ghost_values()
        # Finally assign that to first component of ODE solution; there
        # other components are states. To (0) from ...
        to_ODE_from_P1.assign(ode_solution.sub(0), p_ode)

        # Check
        print(np.any(np.isnan(wh.vector().get_local())), np.any(np.isinf(wh.vector().get_local())))

        # IO
        potential_assigner.assign(uh, wh.sub(1))
        # probe_values = probes.sample(uh)

        # if pyMPI.COMM_WORLD.rank == 0:
        #     table = np.c_[table, np.r_[float(t1), probe_values[:, 0]]]

# Final dump
# if pyMPI.COMM_WORLD.rank == 0:
#     np.savetxt('probe_readings.txt', table,
#                header='First row it time, remaining are probe readings of potentials')

del pc, ksp
