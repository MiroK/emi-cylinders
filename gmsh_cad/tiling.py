from dolfin import CompiledSubDomain, Mesh, MeshEditor, Timer, info
from test_periodic import compute_vertex_periodicity
import numpy as np


def TileMesh(tile, shape):
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
    # Evolve
    while shape:
        tile, shape = glue(tile, shape)

    return tile

        
def glue(tile, shape):
    '''Evolve tile along the last exis'''
    axis, gdim = len(shape) - 1, tile.geometry().dim()
    assert gdim > axis >= 0

    print 'Glue along', axis
    # We're done evolving if only one tile is to be plae in the axis dir
    if shape[axis] == 1: return tile
    # The first time we seed the algorithm with mesh data
    x = tile.coordinates()
    min_x = np.min(x, axis=0)
    max_x = np.max(x, axis=0)

    min_axis = min_x[axis]
    max_axis = max_x[axis]
    shift = max_axis - min_axis
    # The delta is the same
    shift_x = np.zeros(gdim); shift_x[axis] = shift
    
    to_master = lambda x, shift=shift_x: x - shift
    # This first time we have to compute data from mesh
    TOL = 1E-9
    master = CompiledSubDomain('near(x[i], A, tol)', i=axis, A=min_axis, tol=TOL)
    slave = CompiledSubDomain('near(x[i], A, tol)', i=axis, A=max_axis, tol=TOL)

    # Always from max to min and extend beyond max
    _, vertex_mapping = compute_vertex_periodicity(tile, master, slave, to_master)
    # From mesh
    cells = tile.cells()
    # Master never changes
    master_vertices = vertex_mapping.values()
    refine = shape[axis]
    while refine > 1:
        print refine
        n = len(x)
        # To make the tile piece we add all but the master vertices
        new_vertices = np.fromiter(sorted(set(range(n)) - set(master_vertices)), dtype=int)
        # Verices of the glued tiles
        x = np.r_[x, x[new_vertices] + shift_x]

        translate = np.arange(n)
        # Offset the free
        translate[new_vertices] = n + np.arange(len(new_vertices))
        # Those at master positions take slave values
        slave_vertices = vertex_mapping.keys()
        translate[master_vertices] = slave_vertices

        # Cells of the glued tiles
        new_cells = np.zeros_like(cells)
        new_cells.ravel()[:] = map(translate.__getitem__, cells.flat)

        cells = np.r_[cells, new_cells]

        # Update the periodicty mapping - slaves are new
        vertex_mapping = dict(zip(master_vertices, map(translate.__getitem__, slave_vertices)))

        refine /= 2

    # Done evolving, we can write data
    tdim = tile.topology().dim()
    ctype = str(tile.ufl_cell())

    # We have a new tile
    tile = make_mesh(x, cells, ctype, tdim, gdim)
    # to be evolved in the new shape
    shape = shape[:-1]

    return tile, shape


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
        tile = UnitSquareMesh(2, 2)
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


    
