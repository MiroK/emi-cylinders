from tiling_cpp import (fill_mesh, fill_mesh_function, fill_mesh_valuecollection,
                        fill_mf_from_mvc)
from test_periodic import compute_vertex_periodicity

from dolfin import (CompiledSubDomain, Mesh, MeshEditor, MeshFunction, MeshValueCollection,
                    Timer, info, SubsetIterator)
from collections import defaultdict
from itertools import izip
import numpy as np

from test import make_tile, powers2
# FIXME: 
#
# mesh connectivity to identify entities takes the most time
#
# if make_mesh stored MeshValueCollections and there was way to 
# make mvc to MeshFunctions then a lot of space could be saved


def merge(x, x_cells, shift, master_vertices, slave_vertices, y=None, y_cells=None):
    
    if y is None and y_cells is None:
        y, y_cells = x, x_cells

    n = len(y)
    # To make the tile piece we add all but the master vertices
    new_vertices = set(range(n))
    new_vertices.difference_update(master_vertices)
    new_vertices = np.fromiter(new_vertices, dtype=int)
    assert np.all(np.diff(new_vertices) > 0)
        
    # Offset the free
    translate = np.empty(n, dtype=int)
    translate[new_vertices] = len(x) + np.arange(len(new_vertices))
    # Those at master positions take slave values
    translate[master_vertices] = slave_vertices

    # Verices of the glued tiles
    x = np.vstack([x, y[new_vertices] + shift])

    # Cells of the glued tiles
    new_cells = np.empty_like(y_cells)
    new_cells.ravel()[:] = translate[y_cells.flat]

    x_cells = np.vstack([x_cells, new_cells])

    return x, x_cells, translate


def TileMesh(tile, shape, mesh_data={}, TOL=1E-9):
    '''
    [tile tile;
     tile tile;
     tile tile;
     tile tile]

    The shape is an ntuple describing the number of pieces put next 
    to each other in the i-th axis. mesh_data : (tdim, tag) -> [entities] 
    is the way to encode mesh data of the tile.
    '''
    # Sanity for glueing
    gdim = tile.geometry().dim()
    assert len(shape) <= gdim, (shape, gdim)
    # While evolve is general mesh writing is limited to simplices only (FIXME)
    # so we bail out early
    assert str(tile.ufl_cell()) in ('interval', 'triangle', 'tetrahedron')

    t = Timer('evolve')
    # Do nothing
    if all(v == 1 for v in shape):
        return tile, mesh_data

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
    cells = np.empty((tile.num_cells(), tile.ufl_cell().num_vertices()), dtype='uintp')
    cells.ravel()[:] = tile.cells().flat
        
    # Evolve the mesh by tiling along the last direction in shape
    while shape:
        x, cells, vertex_mappings, shape = evolve(x,
                                                  cells,
                                                  vertex_mappings,
                                                  shape, shifts_x, mesh_data=mesh_data)
    info('\tEvolve took %g s ' % t.stop())

    # Mesh data is evolved, (x cells) -> to mesh
    mesh = make_mesh(x, cells, tdim=tile.topology().dim(), gdim=gdim)

    return mesh, mesh_data

        
def evolve(x, cells, vertex_mappings, shape, shifts_x, mesh_data={}):
    '''Evolve tile along the last axis'''
    axis, gdim = len(shape) - 1, x.shape[1]
    assert gdim > axis >= 0

    # We're done evolving if only one tile is to be plae in the axis dir
    if shape[axis] == 1:
        vertex_mappings.pop()  # No longer needed
        shifts_x.pop()  # Do not touch x and cells
        return x, cells, vertex_mappings, shape[:-1]

    # We have x, cells, ... representing a single cell. To get the tile
    # repeated shape[::-1] times I don't simply stack next to each other.
    # Consider 1+1+1+1+1+1+1+1 has 7 additions
    # But      (4+4) where 4 = 2+2 where 2 = 1+1 is three! So providided
    # that `+` scales reasonably with input size this is more efficient.
    # NOTE that here + is really unary. For input other than power of 2
    # it is necessary to extend it to binary, e.g. 5 = (2+2)+1
    refine = shape[axis]
    # Get the sizes of tiles (the real sizes are powers of two of these)
    tile_sizes = powers2(refine)
    # Iff refine was a power of two then there is no need to merge as
    # at the end of tile evolution is the final mesh data. Otherwise
    # we develop the tiles in their LOCAL numbering => they each start
    # from their copy
    vertex_mapping, shift_x = vertex_mappings.pop(), shifts_x.pop()
    master_vertices = vertex_mapping.values()
    slave_vertices = vertex_mapping.keys()

    tiles = []
    # Are we even or odd (to keep the initial tile)
    if 0 in tile_sizes:
        tiles = [make_tile(x, cells,
                           master_vertices, slave_vertices, vertex_mappings,
                           mesh_data)]
        # Done with this size
        tile_sizes.pop()

    # We need to evolve tiles of sizes which are power of 2. The idea
    # is to get them while evolving to the largest one
    size = 1
    max_size = max(tile_sizes)

    target_size = tile_sizes.pop()
    # The body of the loop corresponds to unary + 
    while size <= max_size:
        x, cells, translate = merge(x, cells, shift_x, master_vertices, slave_vertices)

        # For the directions that do not evolve we add the periodic pairs
        for vm in vertex_mappings:
            keys, values = np.array(vm.items()).T
            vm.update(dict(izip(translate[keys], translate[values])))

        # Update the periodicty mapping - slaves are new
        slave_vertices = translate[slave_vertices]
        
        # Add the entities defined in terms of the vertices
        if mesh_data:
            mesh_data = evolve_data(mesh_data, translate)

        # We might be at a needed size
        if size == target_size:
            tiles.append(make_tile(x, cells,
                                   master_vertices, slave_vertices, vertex_mappings,
                                   mesh_data))
            if tile_sizes:
                target_size = tile_sizes.pop()
      
        # Iterate
        size += 1
        shift_x *= 2
    assert not tile_sizes

    tile0 = tiles.pop()
    if not tiles:
        return tile0.coords, tile0.cells, tile0.mappings, shape[:-1]

    x, cells, vertex_mappings = tile0.coords, tile0.cells, tile0.mappings
    slave_vertices = tile0.slave_vertices

    # Now the binary +
    dx = shift_x / 2**max_size
    shift_x = 0

    shifts = (2**p for p in powers2(refine))
    while tiles:
        next_tile = tiles.pop()

        shift_x = shift_x + dx*next(shifts)        
        x, cells, translate = merge(x, cells, shift_x,
                                    next_tile.master_vertices, slave_vertices,
                                    next_tile.coords, next_tile.cells)

        # Data
        for vm, next_vm in zip(vertex_mappings, next_tile.mappings):
            keys, values = np.array(next_vm.items()).T
            vm.update(dict(izip(translate[keys], translate[values])))

        # Data evolve
        if mesh_data:
            mesh_data = evolve_data(mesh_data, translate, next_tile.data)
        
        # For next round
        slave_vertices = translate[next_tile.slave_vertices]
    # Discard data not needed in next evolution
    return x, cells, vertex_mappings, shape[:-1]


def evolve_data(data, mapping, other=None):
    '''
    If mapping holds (tdim, tag) -> [tuple of indices]) where indices are 
    w.r.t of old numbering and mapping is old to new we simply add the mapped 
    entities.
    '''
    if other is None: other = data
    
    for key in data.keys():
        old = data[key]
            
        new = np.empty_like(other[key])
        new.ravel()[:] = mapping[other[key].flat]
        data[key] = np.vstack([old, new])
    return data


def make_mesh(coordinates, cells, tdim, gdim):
    '''Mesh by MeshEditor from vertices and cells'''
    mesh = Mesh()
    assert mesh.mpi_comm().tompi4py().size == 1

    fill_mesh(coordinates.flatten(), cells.flatten(), tdim, gdim, mesh)
    
    return mesh


def mf_from_data(mesh, data):
    '''Build tdim -> mesh function from the data of TileMesh'''
    return _mx_from_data(mesh, data,
                         fill=fill_mesh_function,
                         init_container=lambda m, t: MeshFunction('size_t', m, t, 0))


def mvc_from_data(mesh, data):
    '''Build tdim -> mesh value collection from data of TileMesh'''
    return _mx_from_data(mesh, data,
                         fill=fill_mesh_valuecollection,
                         init_container=lambda m, t: MeshValueCollection('size_t', m, t))


def groupby(pairs, index):
    '''Organize pairs by pairs[index]'''
    groups = defaultdict(list)
    for pair in pairs: groups[pair[index]].append(pair)

    for item in groups.iteritems():
        yield item

 
def _mx_from_data(mesh, data, fill, init_container):
    '''Fill the contained over mesh by data'''
    assert mesh.mpi_comm().tompi4py().size == 1

    containers = {}
    # We have define entities in terms of vertex numbering
    # Order keys such by tdim (the first key)
    for tdim, keys in groupby(data.keys(), 0):
        # So we'll be getting the entity index by lookup
        mesh.init(tdim)
        mesh.init(0, tdim)
        # Build the meshfunction from data
        f = init_container(mesh, tdim)
        for key in keys:
            indices = data[key]
            # These entity indices get the 'color'
            fill(mesh, indices.flatten(), tdim, key[1], f)
        containers[tdim] = f

    return containers


def as_meshf(mvc, init_value=0):
    '''Make a mesh function out of mesh value collection'''
    if isinstance(mvc, (tuple, list)):
        return [as_meshf(x, init_value) for x in mvc]

    if isinstance(mvc, dict):
        return dict(zip(mvc.keys(), as_meshf(mvc.values())))

    # Base case
    mesh_f = MeshFunction('size_t', mvc.mesh(), mvc.dim(), init_value)
    fill_mf_from_mvc(mvc, mesh_f)

    return mesh_f


def load_data(mesh, h5_file, data_set, dim, data):
    '''
    Fill the data dictionary with data_set representing mesh function with 
    dim over mesh read from h5_file according to key spec expected by tiling 
    algorithm.
    '''
    mf = MeshFunction('size_t', mesh, dim, 0)
    h5_file.read(mf, data_set)
            
    # Data to evolve
    mesh.init(dim, 0)
    e2v = tile.topology()(dim, 0)

    tags = set(mf.array())
    # Don't evolve zero - we initialize to it
    if 0 in tags: tags.remove(0)
    info('%s evolves tags %r' % (data_set, tags))

    for tag in tags:
        data[(dim, tag)] = np.array([e2v(e.index()) for e in SubsetIterator(mf, tag)],
                                    dtype='uintp')
    return data

# ------------------------------------------------------------------------------

if __name__ == '__main__':
    from dolfin import mpi_comm_world, HDF5File, Timer, File
    import argparse, os

    parser = argparse.ArgumentParser(description='Put n tiles in x axis, m in y axis.')
    parser.add_argument('tile', type=str, help='H5 file that is the file')
    parser.add_argument('-n', type=int, default=1)
    parser.add_argument('-m', type=int, default=1)
    parser.add_argument('-facet_tags', type=str, default='surfaces',
                        help='name under which H5 stores facet tags')
    parser.add_argument('-cell_tags', type=str, default='volumes',
                        help='name under which H5 stores volume tags')

    save_pvd_parser = parser.add_mutually_exclusive_group(required=False)
    save_pvd_parser.add_argument('--save_pvd', dest='save_pvd', action='store_true')
    save_pvd_parser.add_argument('--no_save_pvd', dest='save_pvd', action='store_false')
    parser.set_defaults(save_pvd=False)

    args = parser.parse_args()

    # Some sanity
    root, ext = os.path.splitext(args.tile)
    assert ext == '.h5'

    shape = (args.n, args.m)
    assert all((((v & (v - 1)) == 0) and v > 0) for v in shape)
    
    # Load the tile mesh
    comm = mpi_comm_world()
    h5 = HDF5File(comm, args.tile, 'r')
    tile = Mesh()
    h5.read(tile, 'mesh', False)

    data = {}
    cell_dim = tile.topology().dim()
    facet_dim = cell_dim - 1

    if args.facet_tags: 
        data = load_data(tile, h5, args.facet_tags, facet_dim, data)
    
    if args.cell_tags: 
        data = load_data(tile, h5, args.cell_tags, cell_dim, data)

    t = Timer('tile')
    mesh, mesh_data = TileMesh(tile, shape, mesh_data=data)
    info('\nTiling took %g s; nvertices %d, ncells %d' % (t.stop(),
                                                          mesh.num_vertices(),
                                                          mesh.num_cells()))

    # Saving
    t = Timer('save')
    h5_file = '%s_%d_%d.h5' % (root, shape[0], shape[1])
        
    out = HDF5File(mesh.mpi_comm(), h5_file, 'w')
    out.write(mesh, 'mesh')
    
    tt = Timer('data')
    # To mesh functions
    if mesh_data:
        mfs = mf_from_data(mesh, mesh_data)

        for dim, name in (zip((facet_dim, cell_dim), (args.facet_tags, args.cell_tags))):
            if name:
                out.write(mfs[dim], name)
            
                if args.save_pvd:
                    File('%s_%d_%d_%s.pvd' % (root, shape[0], shape[1], name)) << mfs[dim]
                
    info('\t\tGetting data as MeshFoo took %g s' % tt.stop())
    
    info('\tSaving took %g' % t.stop())
