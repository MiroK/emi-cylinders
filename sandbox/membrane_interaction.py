from __future__ import print_function
from dolfin import *
import numpy as np

parameters['ghost_mode'] = 'shared_facet'

# We consider EMI formulation using RT1-DG0-DLT0 elements where the transmembrane
# potential lives on ALL facets of the mesh. Meanwhile the ODE problem,
# whose input is the transmembrane potential is solved in the vertices,
# more precisely, ODE ouput lives on P1 space. However, the P1 space serves
# as an input to the EMI step. A P1 function is okay input for EMI and
# in this sense ODE -> EMI step is easy. Here we handle the case from
# EMI -> ODE
def emi_to_ode_operator(membrane, ids, mean='arithmetic'):
    '''
    Produce A such that for u a vector of coefs of DLT0 function A*u 
    gives coefs of P1 function.
    '''
    mesh = membrane.mesh()

    V = FunctionSpace(mesh, 'CG', 1)  # ODE
    Q = FunctionSpace(mesh, 'Discontinuous Lagrange Trace', 0)

    p = TrialFunction(Q)
    v = TestFunction(V)

    K = FacetArea(mesh)
        
    a = 0
    # Take advantage of the fact u*q over facet is the same as u(mid)*facet_area
    n = Constant(mesh.topology().dim())
    if mean == 'arithmetic':
        for tag in ids:
            dMeas = dS(subdomain_data=membrane, subdomain_id=tag)
            a += (n/K('+'))*inner(v('+'), p('+'))*dMeas
            
            dMeas = ds(subdomain_data=membrane, subdomain_id=tag)
            a += (n/K)*inner(v, p)*dMeas
    # Clement
    else:
        for tag in ids:
            dMeas = dS(subdomain_data=membrane, subdomain_id=tag)
            a += n*inner(v('+'), p('+'))*dMeas
            
            dMeas = ds(subdomain_data=membrane, subdomain_id=tag)
            a += n*inner(v, p)*dMeas
        
    A = assemble(a)

    mat = as_backend_type(A).mat()
    x = mat.getRowSum()
    x.reciprocal()
    mat.diagonalScale(x)

    return A


def face_code(axis, value, bounds):
    '''x[axis] == value and other are within bounds'''
    code = ['near(x[%d], %g)' % (axis, value)]
    axes = list(range(len(bounds)+1))
    # Now the code for remaining axis
    for i, interval in zip((axes[:axis] + axes[axis+1:]), bounds):
        face = (interval[0], i, i, interval[1])
        code.append('(%g - tol < x[%d]) && (x[%d] < %g + tol)' % face)
        # All this must true for the face to be marked
    code = ' && '.join(code)
    
    return code


def hypercube_bdry(intervals):
    '''CompiledSubDomain code marking boundary of I0 x I1 ... box'''
    dim = len(intervals)
    assert all(l < h for l, h in intervals)
    
    axes = list(range(dim))

    bdry_code = []
    for axis, interval in enumerate(intervals):
        for bound in interval:
            bdry_code.append(face_code(axis, bound, intervals[:axis]+intervals[axis+1:]))
    # Any face is a boudary
    bdry_code = ' || '.join(['(%s)' % c for c in bdry_code])

    return bdry_code


def check_convergence(meshes, membrane, f, mean):
    '''Convergence on meshes'''
    hs, errors = [], []
    for mesh in meshes:
        facet_f = MeshFunction('size_t', mesh, mesh.topology().dim()-1, 0)
        membrane.mark(facet_f, 1)
        File('f.pvd') << facet_f
        
        V = FunctionSpace(mesh, 'CG', 1)
        Q = FunctionSpace(mesh, 'Discontinuous Lagrange Trace', 0)

        q = interpolate(f, Q)
        u = Function(V)
    
        EMI2ODE = emi_to_ode_operator(facet_f, (1, ), mean=mean)
        EMI2ODE.mult(q.vector(), u.vector())
        as_backend_type(u.vector()).update_ghost_values()  # Keep this

        e = inner(u-f, u-f)*ds(subdomain_data=facet_f, subdomain_id=1)
        e += inner(avg(u-f), avg(u-f))*dS(subdomain_data=facet_f, subdomain_id=1)

        error = sqrt(abs(assemble(e)))
        errors.append(error)

        hmin = MPI.min(mesh.mpi_comm(), mesh.hmin())
        hs.append(hmin)
        
    errors = np.array(errors)
    hs = np.array(hs)
    rates = np.r_[np.nan, np.log(errors[1:]/errors[:-1])/np.log(hs[1:]/hs[:-1])]

    return np.c_[hs, errors, rates], u

# --------------------------------------------------------------------

if __name__ == '__main__':

    mean = 'arithmetic'
    #############
    # Two d
    #############
    meshes = lambda : (UnitSquareMesh(n, 2*n) for n in (4, 8, 16, 32))

    f = Constant(2)
    membrane = CompiledSubDomain(hypercube_bdry([[0.25, 0.75], [0.25, 0.75]]),
                                 tol=1E-10)
    # This is a must
    table, _ = check_convergence(meshes(), membrane, f, mean)
    assert np.all(table[:, 1] < 1E-10), table  # Should be exact
    
    f = Expression('x[0]+x[1]', degree=1)
    membrane = CompiledSubDomain(hypercube_bdry([[0.25, 0.75], [0.25, 0.75]]),
                                 tol=1E-10)
    # More advanced, but still simple surface
    table, _ = check_convergence(meshes(), membrane, f, mean)
    assert table[-1, 2] > 1, table  # Should converge 
    
    # More involved surface
    f = Expression('x[0]+x[1]', degree=1)
    membrane = ' || '.join([hypercube_bdry([[0.25, 0.75], [0.25, 0.75]]),
                            face_code(0, 0.5, [[0.25, 0.75]]),
                            face_code(1, 0.5, [[0.25, 0.75]])])

    membrane = CompiledSubDomain(membrane, tol=1E-10)
    # More advanced, but still simple surface
    table, u = check_convergence(meshes(), membrane, f, mean)
    assert table[-1, 2] > 1, table  # Should converge 

    #############
    # Three d
    #############
    meshes = lambda : (UnitCubeMesh(*(n, )*3) for n in (4, 8, 16, 32))

    f = Constant(2)
    membrane = CompiledSubDomain(hypercube_bdry([[0.25, 0.75], [0.25, 0.75], [0.25, 0.75]]),
                                 tol=1E-10)
    # This is a must
    table, _ = check_convergence(meshes(), membrane, f, mean)
    assert np.all(table[:, 1] < 1E-10)  # Should be exact

    f = Expression('x[0]*x[0]+x[1]*x[1]+x[2]*x[2]', degree=1)
    membrane = CompiledSubDomain(hypercube_bdry([[0.25, 0.75], [0.25, 0.75], [0.25, 0.75]]),
                                 tol=1E-10)
    # More advanced, but still simple surface
    table, u = check_convergence(meshes(), membrane, f, mean)
    assert table[-1, 2] > 1, table  # Should converge 

    # Other with cut
    membrane = ' || '.join([hypercube_bdry([[0.25, 0.75], [0.25, 0.75], [0.25, 0.75]]),
                            face_code(0, 0.5, [[0.25, 0.75], [0.25, 0.75]]),
                            face_code(1, 0.5, [[0.25, 0.75], [0.25, 0.75]]),
                            face_code(2, 0.5, [[0.25, 0.75], [0.25, 0.75]])])

    membrane = CompiledSubDomain(membrane, tol=1E-10)
    # More advanced, but still simple surface
    table, u = check_convergence(meshes(), membrane, f, mean)
    assert table[-1, 2] > 1, table  # Should converge 

    # # What is the difference between these guys
    # if MPI.size(mpi_comm_world()) == 1:
    #     set_log_level(WARNING)
        
    #     comm = mpi_comm_world()
    #     h5 = HDF5File(comm, '../gmsh_cad/tile_1_narrow.h5', 'r')
    #     mesh = Mesh()
    #     h5.read(mesh, 'mesh', False)

    #     surfaces = MeshFunction('size_t', mesh, mesh.topology().dim()-1, 0)
    #     h5.read(surfaces, 'surfaces')

    #     f = Expression('x[0] + 2*x[1] -3*x[2]', degree=1)
    #     V = FunctionSpace(mesh, 'CG', 1)
    #     Q = FunctionSpace(mesh, 'Discontinuous Lagrange Trace', 0)

    #     q = interpolate(f, Q)
        
    #     u_mean = Function(V)
    #     EMI2ODE = emi_to_ode_operator(surfaces, (1, ), mean='arithmetic')
    #     EMI2ODE.mult(q.vector(), u_mean.vector())
    #     as_backend_type(u_mean.vector()).update_ghost_values()

    #     u_clem = Function(V)
    #     EMI2ODE = emi_to_ode_operator(surfaces, (1, ), mean='clement')
    #     EMI2ODE.mult(q.vector(), u_clem.vector())
    #     as_backend_type(u_clem.vector()).update_ghost_values()  # Keep this

    #     errors = []
    #     for u in (u_mean, u_clem):
    #         e = inner(u-f, u-f)*ds(subdomain_data=surfaces, subdomain_id=1)
    #         e += inner(avg(u-f), avg(u-f))*dS(subdomain_data=surfaces, subdomain_id=1)
    #         errors.append(assemble(e))

    #     print('clem', errors[0], 'mean', errors[1])

    #     mesh.init(2)
    #     mesh.init(2, 0)

    #     f2v = mesh.topology()(2, 0)
    #     x = mesh.coordinates()
        
    #     pointwise_error = lambda u, f=f, x=x, surfaces=surfaces: (
    #         abs(u(vertex)-f(vertex))/abs(f(vertex))
    #         for facet in SubsetIterator(surfaces, 1)
    #         for vertex in x[f2v(facet.index())])

    #     # Where is 20% error
    #     for facet in SubsetIterator(surfaces, 1):
    #         if max(abs(u(vertex)-f(vertex))/abs(f(vertex))
    #                for vertex in x[f2v(facet.index())]) > 0.2:  
    #             surfaces[facet] = 2
    #     File('debug.pvd') << surfaces
        
    #     stats = lambda u: [g(pointwise_error(u))
    #                        for g in (min, max, lambda f: np.mean(np.fromiter(f, dtype=float)))]
        
    #     print('clem', stats(u_clem))
    #     print('mean', stats(u_mean))
