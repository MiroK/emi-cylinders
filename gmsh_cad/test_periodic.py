from dolfin import MeshFunction
import numpy as np


def compute_vertex_periodicity(mesh, master, slave, to_master):
    '''Compute mapping from slave vertices to master vertices'''
    tdim = mesh.topology().dim()
    f = MeshFunction('size_t', mesh, tdim-1, 0)
    master.mark(f, 2)
    slave.mark(f, 3)

    f2v = mesh.topology()(tdim-1, 0)
    master_vertices = list(set(sum((f2v(f.index()).tolist() for f in SubsetIterator(f, 2)), [])))
    slave_vertices = set(sum((f2v(f.index()).tolist() for f in SubsetIterator(f, 3)), []))

    assert len(master_vertices) == len(slave_vertices), (len(master_vertices), len(slave_vertices))

    master_vertex_x = x[master_vertices]

    mapping = {}
    error = 0
    for s in slave_vertices:
        xs = x[s]
        mapped = to_master(xs)
        # Local to master_vertex_x
        dist = np.sqrt(np.sum((master_vertex_x - mapped)**2, axis=1))
        mapped_index = np.argmin(dist)
        # Wrt to vertex numbering
        m = master_vertices[mapped_index]
        mapping[s] = m
        error = max(error, dist[mapped_index])

    return error, mapping
    
# --------------------------------------------------------------------

if __name__ == '__main__':
    from dolfin import *

    mesh_file = 'tile_2x2.h5'
    mesh_file = 'tile_1.h5'

    comm = mpi_comm_world()
    h5 = HDF5File(comm, mesh_file, 'r')
    mesh = Mesh()
    h5.read(mesh, 'mesh', False)

    assert comm.tompi4py().size == 1

    x = mesh.coordinates()
    min_ = np.min(x, axis=0)
    max_ = np.max(x, axis=0)

    tol = 1E-9  # The precision in gmsh isn't great, probably why DOLFIN's
    # periodic boundary computation is not working

    # Check x periodicity
    master = CompiledSubDomain('near(x[0], A, tol) && on_boundary', A=min_[0], tol=tol)
    slave = CompiledSubDomain('near(x[0], A, tol) && on_boundary', A=max_[0], tol=tol)

    shift_x = np.array([max_[0]-min_[0], 0, 0])
    to_master = lambda x, shift=shift_x: x - shift

    error, mapping = compute_vertex_periodicity(mesh, master, slave, to_master)
    assert error < 10*tol, error
    print(error)

    # Check y periodicity
    master = CompiledSubDomain('near(x[1], A, tol) && on_boundary', A=min_[1], tol=tol)
    slave = CompiledSubDomain('near(x[1], A, tol) && on_boundary', A=max_[1], tol=tol)

    shift_x = np.array([0, max_[1]-min_[1], 0])
    to_master = lambda x, shift=shift_x: x - shift

    error, mapping = compute_vertex_periodicity(mesh, master, slave, to_master)
    assert error < 10*tol, error
    print(error)

