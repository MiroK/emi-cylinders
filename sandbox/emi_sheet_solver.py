from cbcbeat.cellsolver import CardiacODESolver

from parsimonious import Parsimonious
from emi_pde import pde_components
from probing import PointProbe

from mpi4py import MPI as pyMPI
from dolfin import *
import numpy as np


# Setup 
mesh_file = './tile_1_hein_GMSH307_1_1.h5'
# Get mesh to setup the ode solver
comm = MPI.comm_world  # FIXME!
h5 = HDF5File(comm, mesh_file, 'r')
mesh = Mesh()
h5.read(mesh, 'mesh', False)

# Mesh box dimensions. NOTE: update by hand - waste of resouces to do it with
# mpi
Length = 0.1*1    # In x as multiples of single cell
Width = 0.025*1   # In y
Height = 0.025

# ODE setup
C_m = 0.01 #microF*mm**-2
chi = 200 #mm**-1

# define cell model with zero stimulus
cell_params = Parsimonious.default_parameters()
cell_params['amp'] = 0.0
cell_model = Parsimonious(params=cell_params)

# Define a stimulus current
stim = Expression("t > 5.0 && t < 5.5 && x[0] < length*0.1? 80:0", t=0., length=Length, degree=1)

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
# Electrodes placed symmetrically around x[0] = 10, 4mm apart, 1mm long (with 200 cells, 20mm)
#                                              mid, Length/5, Lenth/20
electrode_length = Length/20.
electrode_distance = Length/5.
e1 = CompiledSubDomain("near(x[1], 0) && x[0] > (c-l)-tol && x[0] < (c+l)+tol && on_boundary",
                       tol = 1e-10,
                       c=Length/2.-electrode_distance/2.-electrode_length/2.,
                       l=electrode_length/2.)

e2 = CompiledSubDomain("near(x[1], 0) && x[0] > (c-l)-tol && x[0] < (c+l)+tol && on_boundary",
                       tol = 1e-10,
                       c=Length/2.+electrode_distance/2.+electrode_length/2.,
                       l=electrode_length/2.)

# All currents are scaled with Cm*chi, J_app given in uA*mm**-2:
# Value must be vector A*n where A is magniture and n is the normal (outer,
# matching that of e0_domain). Here it normal to y == 0
s1 = Expression(('0', "t < 2.0 ? -J_s1*C_m*chi:0", '0'), t=0.0,
                J_s1=28, C_m=C_m, chi = chi, degree=1)

s2 = Expression(('0', "t > t_s2 && t< t_s2 + 2.0 ? -J_s2*C_m*chi:0", '0'), t=0.0,
                t_s2 = 190, J_s2=28, C_m=C_m, chi = chi, degree=1)

# List value for convenience of update (if they are t-dep)
electrode_values = [s1, s2]
# Domain listed for symmetry
electrode_domains = [e1, e2]

emi_parameters = {'C_int': 3*C_m,
                  'C_ext': 2*C_m,
                  'C_m': C_m,
                  'dt': 0.1,
                  # Electrodes [(where, value)]
                  'bc_pairs': list(zip(electrode_domains, electrode_values)),
                  # Potentials as P1 function on the mesh are wired up with the EMI = p0
                  'p0': ode_solver.vs[0],
                  'I_ion': Constant(0.0)}  # FIXME: what should this be?

# Pieces for emi are a, L, W, bcs (so assembly) and transfer operator
emi_pieces = pde_components(mesh, h5, emi_parameters)

# Since these are meant to be simulations with fixed time step, the matrix
# is assembled once
a, L, bcs = (emi_pieces[key] for key in ('a', 'L', 'bcs'))

A, b = PETScMatrix(), PETScVector()

emi_assembler = SystemAssembler(a, L, bcs)
emi_assembler.assemble(A)


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


solver = LUSolver()
solver.set_operator(A)
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

probes = PointProbe(uh, probe_points)
# Initial reading
probe_values = probes.sample(uh)

# The idea is to have columns of t and potentials readings
if pyMPI.COMM_WORLD.rank == 0: table = np.r_[interval[0], probe_values[:, 0]]
    
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
        solver.solve(wh.vector(), b)

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
        probe_values = probes.sample(uh)

        if pyMPI.COMM_WORLD.rank == 0:
            table = np.c_[table, np.r_[float(t1), probe_values[:, 0]]]

# Final dump
if pyMPI.COMM_WORLD.rank == 0:
    np.savetxt('probe_readings.txt', table,
               header='First row it time, remaining are probe readings of potentials')
