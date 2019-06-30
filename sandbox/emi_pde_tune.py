from membrane_interaction import emi_to_ode_operator
from bc_marking import tag_if_not
from mpi4py import MPI as pyMPI
from dolfin import *
import numpy as np


parameters['form_compiler']['representation'] = 'uflacs'
parameters['form_compiler']['cpp_optimize'] = True
parameters['form_compiler']['cpp_optimize_flags'] = '-O3 -ffast-math -march=native'
parameters["form_compiler"]["quadrature_degree"] = 3

parameters['ghost_mode'] = 'shared_facet'

# EMI system is defined by conductivity of
# - [C_int, C_ext, C_m] = intra/extracellular spaces and C_m conducivity,
# - [dt] = time step
#   and the Expression that is their value.
# - [p0] A P1 function that represents the trans
# - [I_ion] The ionic current term
def pde_components(mesh, h5, emi_parameters):
    '''
    Ingredients of EMI system are W, a: WxW->R, L: W->R, bcs and operator for
    maping transmembrane potentials from membrane facets to membrane
    vertices (where ODE needs it).
    '''
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

    # The boundary conditions that are applied in this case are insulation
    # on all but the bottom boundary. The bottom is grounded. Note that
    # x and y extreme boundaries have pieces of cells. TODO: for now we
    # avoid these.
    comm = mesh.mpi_comm()
    x_gmin = [comm.allreduce(xi, op=pyMPI.MIN) for xi in mesh.coordinates().min(axis=0)]
    x_gmax = [comm.allreduce(xi, op=pyMPI.MAX) for xi in mesh.coordinates().max(axis=0)]

    # Now the extemities
    bdry_foos = []
    for i, (xi_min, xi_max) in enumerate(zip(x_gmin, x_gmax)):
        bdry_foos.extend([CompiledSubDomain('near(x[%d], %.16f)' % (i, xi_min)),
                          CompiledSubDomain('near(x[%d], %.16f)' % (i, xi_max))])

    bdry_tags = range(2, 2+len(bdry_foos))
    # Carefully introduce new surfaces
    for tag, bdry in zip(bdry_tags, bdry_foos):
        tag_if_not(tag, bdry, surfaces, 0)

    # W.sub(0) bcs set strongly correspond to insulation. Skipping (tau.n)*u on
    # bdry means potential there is weakly zero, grounding.
    dx = Measure('dx', domain=mesh, subdomain_data=volumes)
    dS = Measure('dS', domain=mesh, subdomain_data=surfaces)
    ds = Measure('ds', domain=mesh, subdomain_data=surfaces)

    n = FacetNormal(mesh)

    # Now onto the weak form
    # Electric properties of membrane and interior/exterior
    C_m = Constant(emi_parameters['C_m'])
    cond_int = Constant(emi_parameters['C_int'])
    cond_ext = Constant(emi_parameters['C_ext'])
    # Time step
    dt = Constant(emi_parameters['dt'])

    # The source term comming from ODE
    p0 = emi_parameters['p0']
    # And additional source on the boundary is the ionic current. For simplicity
    I_ion =emi_parameters['I_ion']

    # The system
    a = ((1/cond_int)*inner(sigma, tau)*dx(int_tag)+(1/cond_ext)*inner(sigma, tau)*dx(ext_tag)
         - inner(div(tau), u)*dx(int_tag) - inner(div(tau), u)*dx(ext_tag)
         + inner(p('+'), dot(tau('+'), n('-')))*dS(iface_tag)
         - inner(div(sigma), v)*dx(int_tag) - inner(div(sigma), v)*dx(ext_tag)
         + inner(q('+'), dot(sigma('+'), n('-')))*dS(iface_tag)
         - (C_m/dt)*inner(q('+'), p('+'))*dS(iface_tag))

    L = inner(q('+'), I_ion('+')-(C_m/dt)*p0('+'))*dS(iface_tag)

    # NOTE: in general the cell interfaces might end up on the surface boundary
    a += (inner(p, dot(tau, n))*ds(iface_tag)
          + inner(q, dot(sigma, n))*ds(iface_tag)
          - (C_m/dt)*inner(q, p)*ds(iface_tag))

    L += inner(q, I_ion-(C_m/dt)*p0)*ds(iface_tag)

    # Additional terms to set to zero the dofs of W.sub(2) which are not on
    # the interfaces
    a -= inner(p('+'), q('+'))*dS(not_iface_tag) + inner(p, q)*ds(not_iface_tag)
    L -= inner(Constant(0)('+'), q('+'))*dS(not_iface_tag) + inner(Constant(0), q)*ds(not_iface_tag)

    bcs = []
    # Insulate all the boundaries; by leaving out 7 we are making it grounded
    bcs = [DirichletBC(W.sub(0), Constant((0, 0, 0)), surfaces, tag) for tag in bdry_tags[:-1]]

    # Constrain LM on these facets
    for tag in bdry_tags:
        a -= inner(p('+'), q('+'))*dS(tag) + inner(p, q)*ds(tag)
        L -= inner(Constant(0)('+'), q('+'))*dS(tag) + inner(Constant(0), q)*ds(tag)

    # And finaly a function that will update p0 like guys
    toODE = emi_to_ode_operator(surfaces, (1, ))

    return {'a': a, 'L': L, 'W': W, 'bcs': bcs, 'toODE': toODE}

# --------------------------------------------------------------------

if __name__ == '__main__':
    from slepc4py import SLEPc
    
    mesh_file = '../gmsh_cad/tile_1_narrow.h5'

    comm = MPI.comm_world
    h5 = HDF5File(comm, mesh_file, 'r')
    mesh = Mesh()
    h5.read(mesh, 'mesh', False)  

    emi_parameters = {'C_int': 1,
                      'C_ext': 1.2,
                      'C_m': 0.5,
                      'dt': 0.1,
                      'p0': Constant(1),
                      'I_ion': Constant(1)}

    comps = pde_components(mesh, h5, emi_parameters)

    A, _ = assemble_system(comps['a'], comps['L'], comps['bcs'])

    import numpy as np

    # e = np.sort(np.abs(np.linalg.eigvals(A.array())))
