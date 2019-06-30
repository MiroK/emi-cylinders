from dolfin import MeshFunction, File
import numpy as np


def tag_if_not(tag, subdomain, f, init_tag):
    '''
    Use subdomain to tag entities of f tag-value but only if the 
    current entity value is init_tag.
    '''
    g = MeshFunction('size_t', f.mesh(), f.dim(), 0)
    subdomain.mark(g, tag)

    g_values = g.array()
    f_values = f.array()

    f_values[:] = np.where(np.logical_and(f_values == init_tag, g_values == tag),
                           g_values,
                           f_values)
    
    return f


def electrode_bcs(facet_f, avoid, electrodes):
    '''
    Update facet function such that tag for electrodes marks its surfaces
    only if that surface is not tagged as avoid.
    '''
    mesh = facet_f.mesh()
    
    array = facet_f.array()
    tags = []
    for t, domain in enumerate(electrodes, avoid+1):
        f = MeshFunction('size_t', mesh, facet_f.dim(), 0)
        domain.mark(f, t)

        a = f.array()
        
        array[:] = np.where(np.logical_and(a == t, array != avoid), a, array)
        tags.append(t)

    return facet_f, tags

# --------------------------------------------------------------------

if __name__ == '__main__':
    from dolfin import *

    if False:
        mesh = UnitSquareMesh(64, 64)

        cell_f = MeshFunction('size_t', mesh, mesh.topology().dim()-1, 0)
        CompiledSubDomain('near(x[0], 0) && x[1] > 0.25 -tol && x[1] < 0.75+tol', tol=1E-10).mark(cell_f, 1)
        CompiledSubDomain('near(x[0], 1.0) && x[1] > 0.25 -tol && x[1] < 0.75+tol', tol=1E-10).mark(cell_f, 1)
        CompiledSubDomain('near(x[1], 0.0)', tol=1E-10).mark(cell_f, 1)
        CompiledSubDomain('near(x[1], 1.0)', tol=1E-10).mark(cell_f, 1)

        electrodes = [CompiledSubDomain('near(x[0], 0)', tol=1E-10),
                      CompiledSubDomain('near(x[0], 1)', tol=1E-10)]

        # Basically want 2 only on x[0] == {0, 1} and x[1] < 0.25 or x[1] > 0.75
        electrode_bcs(cell_f, avoid=1, electrodes=electrodes)

        File('cell_f.pvd') << cell_f

        # The real deal
        comm = MPI.comm_world
        h5 = HDF5File(comm, 'tile_1_narrow.h5', 'r')
        mesh = Mesh()
        h5.read(mesh, 'mesh', False)

        surfaces = MeshFunction('size_t', mesh, mesh.topology().dim()-1, 0)
        h5.read(surfaces, 'surfaces')

        print(mesh.coordinates().min(axis=0), mesh.coordinates().max(axis=0))
        
        xmin = mesh.coordinates().min(axis=0)[0]
        left = CompiledSubDomain('near(x[0], a, tol)', a=xmin, tol=1E-10)

        xmax = mesh.coordinates().max(axis=0)[0]
        right = CompiledSubDomain('near(x[0], a, tol)', a=xmax, tol=1E-10)

        electrode_bcs(surfaces, avoid=1, electrodes=[left, right])
        File('foo.pvd') << surfaces

        assert set(surfaces.array()) == set((0, 1, 2, 3))


    comm = MPI.comm_world
    h5 = HDF5File(comm, 'tile_1_narrow.h5', 'r')
    mesh = Mesh()
    h5.read(mesh, 'mesh', False)

    # At this point we only have tags 0 and 1
    surfaces = MeshFunction('size_t', mesh, mesh.topology().dim()-1, 0)
    h5.read(surfaces, 'surfaces')
    assert comm.size > 1 or set(surfaces.array()) == set((0, 1))

    left = CompiledSubDomain('near(x[0], -54)')
    right = CompiledSubDomain('near(x[0], 54)')
    front = CompiledSubDomain('near(x[1], -12)')
    back = CompiledSubDomain('near(x[1], 12)')
    bottom = CompiledSubDomain('near(x[2], -14)')
    top = CompiledSubDomain('near(x[2], 14)')

    subdomains = (left, right, front, back, bottom, top)
    # Now we want to introduce tags
    for tag, subdomain in enumerate(subdomains, 2):
        # It is only okay to overwrite 0
        tag_if_not(tag, subdomain, surfaces, init_tag=0)

    File('foo.pvd') << surfaces
