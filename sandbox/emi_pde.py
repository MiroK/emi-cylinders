from membrane_interaction import emi_to_ode_operator
from bc_marking import electrode_bcs
from dolfin import *


parameters['form_compiler']['representation'] = 'uflacs'
parameters['form_compiler']['cpp_optimize'] = True
parameters['form_compiler']['cpp_optimize_flags'] = '-O3 -ffast-math -march=native'
parameters["form_compiler"]["quadrature_degree"] = 3

parameters['ghost_mode'] = 'shared_facet'

# EMI system is defined by conductivity of
# - [C_int, C_ext, C_m] = intra/extracellular spaces and C_m conducivity,
# - [dt] = time step
# - [bc_pairs] = pairs of compiled subdomains where bcs on current should be applied
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
    
    # For boundary conditions we want to add tags 2, 3, ... for each
    # expression. NOTE: this leaves 1 (the cell surface) unchanged
    if emi_parameters['bc_pairs']:
        bc_subdomains, bc_expressions = list(zip(*emi_parameters['bc_pairs']))
    
        membrane, tags = electrode_bcs(surfaces, avoid=1, electrodes=bc_subdomains)
        # NOTE: while these are expressions for a scalar current FEniCS requires
        # them in a form v*n where n is the boundary normal, i.e. they must be
        # vector valued
        bcs = [DirichletBC(W.sub(0), bc_expr, surfaces, tag) for tag, bc_expr in zip(tags, bc_expressions)]

        # NOTE: okay at this point electrode_bcs has changes some facets from
        # 0 to tag value. The have to be put back to forms
        for tag in tags:
            a -= inner(p('+'), q('+'))*dS(tag) + inner(p, q)*ds(tag)
            L -= inner(Constant(0)('+'), q('+'))*dS(tag) + inner(Constant(0), q)*ds(tag)
    else:
        bcs = []

    # And finaly a function that will update p0 like guys
    toODE = emi_to_ode_operator(surfaces, (1, ))

    return {'a': a, 'L': L, 'W': W, 'bcs': bcs, 'toODE': toODE}

# --------------------------------------------------------------------

if __name__ == '__main__':
    mesh_file = '../gmsh_cad/tile_1_narrow.h5'

    comm = mpi_comm_world()
    h5 = HDF5File(comm, mesh_file, 'r')
    mesh = Mesh()
    h5.read(mesh, 'mesh', False)  

    emi_parameters = {'C_int': 1,
                      'C_ext': 1.2,
                      'C_m': 0.5,
                      'dt': 0.1,
                      'bc_pairs': [(DomainBoundary(), Expression(('0', '0', '1'), degree=0))],
                      'p0': Constant(1),
                      'I_ion': Constant(1)}

    print pde_components(mesh, h5, emi_parameters)
