from dolfin import CompiledSubDomain, Mesh, MeshEditor, Timer, info
from test_periodic import compute_vertex_periodicity
import numpy as np


def TileMesh(tile, shape, TOL=1E-9):
    '''
    [tile tile;
     tile tile;
     tile tile;
     tile tile]

    The shape is an ntuple describing the number of pieces put next 
    to each other in the i-th axis.
    '''
    # All the axis shapes needs to be power of two
    assert all((((v & (v - 1)) == 0) and v > 0) for v in shape)

    gdim = tile.geometry().dim()
    assert len(shape) <= gdim

    # Do nothing
    if all(v == 1 for v in shape): return tile

    # We want to evolve cells, vertices of the mesh using geometry information
    # and periodicity info
    x = tile.coordinates()
    min_x = np.min(x, axis=0)
    max_x = np.max(x, axis=0)
    shifts = max_x - min_x
    
    shifts_x = []  # Geometry
    vertex_mappings = []  # Periodicity
    # Compute geometrical shift for EVERY direction:
    for axis in range(len(shape)):
        shift = shifts[axis]
        # Vector to shift cell vertices
        shift_x = np.zeros(gdim); shift_x[axis] = shift
        shifts_x.append(shift_x)

        # Compute periodicity in the vertices
        to_master = lambda x, shift=shift_x: x - shift
        # Mapping facets
        master = CompiledSubDomain('near(x[i], A, tol)', i=axis, A=min_x[axis], tol=TOL)
        slave = CompiledSubDomain('near(x[i], A, tol)', i=axis, A=max_x[axis], tol=TOL)

        error, vertex_mapping = compute_vertex_periodicity(tile, master, slave, to_master)
        # Fail when exended direction is no periodic
        assert error < 10*TOL, error
        
        vertex_mappings.append(vertex_mapping)
    # The final piece of information is cells
    cells = tile.cells()
        
    # Evolve
    while shape:
        # Evolve is a bang method on vertex_mappings, shifts_x
        x, cells, shape = evolve(x, cells, vertex_mappings, shifts_x, shape)

    # Done evolving, we can write data
    tdim = tile.topology().dim()
    ctype = str(tile.ufl_cell())

    return make_mesh(x, cells, ctype, tdim, gdim)

        
def evolve(x, cells, vertex_mappings, shifts_x, shape):
    '''Evolve tile along the last exis'''
    axis, gdim = len(shape) - 1, x.shape[1]
    assert gdim > axis >= 0

    # We're done evolving if only one tile is to be plae in the axis dir
    if shape[axis] == 1:
        vertex_mappings.pop()  # No longer needed
        shifts_x.pop()  # Do not touch x and cells
        return x, cells, shape[:-1]

    # Use the axis's periodicity and shifting.
    # NOTE: used only here and discarded
    vertex_mapping, shift_x = vertex_mappings.pop(), shifts_x.pop()

    master_vertices = vertex_mapping.values()
    slave_vertices = vertex_mapping.keys()

    refine = shape[axis]
    while refine > 1:    
        n = len(x)
        # To make the tile piece we add all but the master vertices
        new_vertices = np.fromiter(sorted(set(range(n)) - set(master_vertices)), dtype=int)
        # Verices of the glued tiles
        x = np.r_[x, x[new_vertices] + shift_x]

        # NOTE: using getitem and arrays seems to be on par in efficiency
        # with dicts. So then I keep translate as array because efficiency
        translate = np.arange(n)
        # Offset the free
        translate[new_vertices] = n + np.arange(len(new_vertices))
        # Those at master positions take slave values
        translate[master_vertices] = slave_vertices

        # Cells of the glued tiles
        new_cells = np.zeros_like(cells)
        new_cells.ravel()[:] = map(translate.__getitem__, cells.flat)

        cells = np.r_[cells, new_cells]
        # Update the periodicty mapping - slaves are new
        slave_vertices = map(translate.__getitem__, slave_vertices)
        # For the directions that do not evolve we add the periodic pairs
        for vm in vertex_mappings:
            vm.update(dict(zip(translate[vm.keys()], translate[vm.values()])))
        # Iterate
        refine /= 2
        shift_x *= 2
    # Discard data not needed in next evolution
    return x, cells, shape[:-1]


def make_mesh(vertices, cells, ctype, tdim, gdim):
    '''Mesh by MeshEditor from vertices and cells'''
    mesh = Mesh()
    editor = MeshEditor()
    editor.open(mesh, ctype, tdim, gdim)

    editor.init_vertices(len(vertices))
    editor.init_cells(len(cells))

    for vi, x in enumerate(vertices): editor.add_vertex(vi, x)
    
    for ci, c in enumerate(cells): editor.add_cell(ci, *c)

    editor.close()

    return mesh


# --------------------------------------------------------------------

if __name__ == '__main__':
    from dolfin import File, HDF5File, mpi_comm_world, Mesh, UnitSquareMesh

    if False:
        tile = UnitSquareMesh(1, 1)
    else:
        mesh_file = 'tile_2x2.h5'  #'cell_grid_2d.h5'
    
        comm = mpi_comm_world()
        h5 = HDF5File(comm, mesh_file, 'r')
        tile = Mesh()
        h5.read(tile, 'mesh', False)

    for n in (2, 4, 8):
        t = Timer('x')
        mesh = TileMesh(tile, (n, n))
        info('\nTiling took %g s' % t.stop())
        
    File('test.pvd') << mesh


    
