from dolfin import MeshFunction, File
import numpy as np


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

    File('xxx.pvd') << facet_f
        
    return facet_f, tags

# --------------------------------------------------------------------

if __name__ == '__main__':
    from dolfin import *
    
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
    comm = mpi_comm_world()
    h5 = HDF5File(comm, '../gmsh_cad/tile_1_narrow.h5', 'r')
    mesh = Mesh()
    h5.read(mesh, 'mesh', False)

    surfaces = MeshFunction('size_t', mesh, mesh.topology().dim()-1, 0)
    h5.read(surfaces, 'surfaces')

    xmin = mesh.coordinates().min(axis=0)[0]
    left = CompiledSubDomain('near(x[0], a, tol)', a=xmin, tol=1E-10)

    xmax = mesh.coordinates().max(axis=0)[0]
    right = CompiledSubDomain('near(x[0], a, tol)', a=xmax, tol=1E-10)
    
    electrode_bcs(surfaces, avoid=1, electrodes=[left, right])
    File('foo.pvd') << surfaces

    assert set(surfaces.array()) == set((0, 1, 2, 3))
