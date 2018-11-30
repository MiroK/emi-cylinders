from cbcbeat.cellsolver import CardiacODESolver
from parsimonious import Parsimonious
from dolfin import *


parameters["form_compiler"]["representation"] = "uflacs"
parameters["form_compiler"]["cpp_optimize"] = True
flags = ["-O3", "-ffast-math", "-march=native"]
parameters["form_compiler"]["cpp_optimize_flags"] = " ".join(flags)
parameters["form_compiler"]["quadrature_degree"] = 3


mesh = UnitSquareMesh(32, 32)

cell_params = Parsimonious.default_parameters()
cell_params['amp'] = 0.0
cell_model = Parsimonious(params=cell_params)

stimulus = Constant(1)

# FIXME: Set time step 
ode_solver = CardiacODESolver(mesh, time=Constant(0), model=cell_model,
                              I_s=stimulus)

dt_fem = 1E-2
dt_ode = 1E-3

interval = (0.0, 1.0)
    
fem_ode_sync = int(dt_fem/dt_ode)
# NOTE: a generator; nothing is computed so far
ode_solutions = ode_solver.solve(interval, dt_ode)  # Potentials only

step_count = 0
for ((t0, t1), ode_solution) in ode_solutions:
    step_count += 1
    info('Time is (%g, %g)' % (t0, t1))
    if step_count == fem_ode_sync:
        step_count = 0
        print '\t', ode_solution.vector().norm('l2')
        # ODE -> p0_neuron
        #    p0_neuron.assign(ode_solution)
        #    # Upscale p0_neuron->p0
        #    assign_toQ_fromQ_neuron(p0, p0_neuron)
        
        #    # We could have changing in time simulation
        #    if 't' in site_current.user_parameters:
        #        site_current.t = float(t1)
        #    # Assemble right-hand side (changes with time, so need to reassemble)                
        #    assembler.assemble(b)  # Also applies bcs
      
        #    # New (sigma, u, p) ...
        #    info('\tSolving linear system of size %d' % A.size(0))
        #    info('\tNumber of true unknowns %d' % (A.size(0) - ncstr_dofs))
        #    la_solver.solve(w.vector(), b)

        #    # Update u_out and current_out for output
        #    toV_fromW1.assign(u_out, w.sub(1))
        #    # NOTE: the current_form is points to w which has been updated by solve
        #    w_aux.vector()[:] = assemble(current_form)
        #    toQ_fromW2.assign(current_aux, w_aux.sub(2))
        #    assign_toQ_neuron_fromQ(current_out, current_aux)
        #
        #    yield t1, u_out, current_out, A.size(0)
        #    
        #    # Now transfer the new transm potential down to ode ...
        #    toQ_fromW2.assign(p0, w.sub(2))         # Compt to Q
        #    assign_toQ_neuron_fromQ(p0_neuron, p0)  # To membrane space
        #    ode_solution.assign(p0_neuron)         # As IC for ODE
