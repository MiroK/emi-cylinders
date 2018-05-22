from dolfin import MeshFunction, info
import numpy as np


def compute_entity_periodicity(tdim, mesh, master, slave, to_master):
    '''Mapping from slave entities of tdim to master'''
    assert 0 <= tdim < mesh.topology().dim()

    _, vertex_mapping = compute_vertex_periodicity(mesh, master, slave, to_master)
    # Done for vertices
    if tdim == 0:
        return vertex_mapping

    f = MeshFunction('size_t', mesh, tdim, 0)
    master.mark(f, 2)
    slave.mark(f, 3)

    master_entities = (e.index() for e in SubsetIterator(f, 2))
    slave_entities = (e.index() for e in SubsetIterator(f, 3))
        
    mesh.init(tdim, 0)
    e2v = mesh.topology()(tdim, 0)
    # Define in terms of vertices. Invert for loopup
    master_vertices = {tuple(sorted(e2v(e))): e for e in master_entities}
    # For slave we define via mapped vertices
    slave_entities = {e: tuple(sorted(vertex_mapping[v] for v in e2v(e))) for e in slave_entities}

    assert len(master_vertices) == len(slave_entities)

    mapping = {}
    while slave_entities:
        s, vertices = slave_entities.popitem()
        
        m = master_vertices[vertices]
        mapping[s] = m
        
        master_vertices.pop(vertices)
    assert not slave_entities and not master_vertices
    
    return mapping
    
    
def compute_vertex_periodicity(mesh, master, slave, to_master):
    '''Compute mapping from slave vertices to master vertices'''
    tdim = mesh.topology().dim()
    f = MeshFunction('size_t', mesh, tdim-1, 0)
    master.mark(f, 2)
    slave.mark(f, 3)

    mesh.init(tdim-1, 0)
    f2v = mesh.topology()(tdim-1, 0)
    master_vertices = list(set(sum((f2v(f.index()).tolist() for f in SubsetIterator(f, 2)), [])))
    slave_vertices = set(sum((f2v(f.index()).tolist() for f in SubsetIterator(f, 3)), []))

    assert len(master_vertices) == len(slave_vertices), (len(master_vertices), len(slave_vertices))

    error, mapping = 0., {}
    while slave_vertices:
        s = slave_vertices.pop()
        xs = x[s]
        mapped = to_master(xs)
        # Local to master_vertex_x
        master_vertex_x = x[master_vertices]

        dist = np.sqrt(np.sum((master_vertex_x - mapped)**2, axis=1))
        mapped_index = np.argmin(dist)
        # Wrt to vertex numbering
        m = master_vertices[mapped_index]
        mapping[s] = m
        error = max(error, dist[mapped_index])

        master_vertices.remove(m)
    assert not slave_vertices and not master_vertices
    
    return error, mapping
    
# --------------------------------------------------------------------

if __name__ == '__main__':
    from dolfin import *

    mesh_file = 'tile_2x2.h5'
    #mesh_file = 'tile_1.h5'

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
    # Just checking
    slave_vertices, master_vertices = map(list, zip(*mapping.items()))
    print(np.linalg.norm(np.sqrt(np.sum((to_master(x[slave_vertices]) - x[master_vertices])**2, axis=1)),
                         np.inf))
    # Mapping for higher entities
    compute_entity_periodicity(1, mesh, master, slave, to_master)
    compute_entity_periodicity(2, mesh, master, slave, to_master)

    # Check y periodicity
    master = CompiledSubDomain('near(x[1], A, tol) && on_boundary', A=min_[1], tol=tol)
    slave = CompiledSubDomain('near(x[1], A, tol) && on_boundary', A=max_[1], tol=tol)

    shift_x = np.array([0, max_[1]-min_[1], 0])
    to_master = lambda x, shift=shift_x: x - shift

    error, mapping = compute_vertex_periodicity(mesh, master, slave, to_master)
    assert error < 10*tol, error
    print(error)
    # Just checking
    slave_vertices, master_vertices = map(list, zip(*mapping.items()))
    print(np.linalg.norm(np.sqrt(np.sum((to_master(x[slave_vertices]) - x[master_vertices])**2, axis=1)),
                         np.inf))

    # Mapping for higher entities
    compute_entity_periodicity(1, mesh, master, slave, to_master)
    compute_entity_periodicity(2, mesh, master, slave, to_master)



