from tiling_cpp import (fill_mesh, fill_mesh_function, fill_mesh_valuecollection,
                        fill_mf_from_mvc)
from test_periodic import compute_vertex_periodicity

from dolfin import CompiledSubDomain, Mesh, MeshEditor, MeshFunction, MeshValueCollection
from collections import defaultdict
from itertools import izip
import numpy as np
import operator

# FIXME: 
#
# mesh connectivity to identify entities takes the most time
#
# if make_mesh stored MeshValueCollections and there was way to 
# make mvc to MeshFunctions then a lot of space could be saved


def TileMesh(tile, shape, mesh_data=None, TOL=1E-9):
    '''
    [tile tile;
     tile tile;
     tile tile;
     tile tile]

    The shape is an ntuple describing the number of pieces put next 
    to each other in the i-th axis. mesh_data : (tdim, tag) -> [entities] 
    is the way to encode mesh data of the tile.
    '''
    # All the axis shapes needs to be power of two
    assert all((((v & (v - 1)) == 0) and v > 0) for v in shape)
    # Sanity for glueing
    gdim = tile.geometry().dim()
    assert len(shape) <= gdim
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
    cells = np.fromiter(tile.cells().flat, dtype='uintp').reshape(tile.cells().shape)
        
    # Evolve
    while shape:
        # Evolve is a bang method on vertex_mappings, shifts_x
        x, cells, shape = evolve(x, cells, vertex_mappings, shifts_x, shape, mesh_data=mesh_data)
    info('\tEvolve took %g s ' % t.stop())

    # Mesh data is evolved, (x cells) -> to mesh
    mesh = make_mesh(x, cells, tdim=tile.topology().dim(), gdim=gdim)

    return mesh, mesh_data

        
def evolve(x, cells, vertex_mappings, shifts_x, shape, mesh_data=None):
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
        new_vertices = np.fromiter(sorted(set(range(n)) - set(master_vertices)),
                                   dtype=int, count=n-len(master_vertices))
        # Verices of the glued tiles
        x = np.vstack([x, x[new_vertices] + shift_x])

        # NOTE: using getitem and arrays seems to be on par in efficiency
        # with dicts. So then I keep translate as array because efficiency
        translate = np.arange(n)
        # Offset the free
        translate[new_vertices] = n + np.arange(len(new_vertices))
        # Those at master positions take slave values
        translate[master_vertices] = slave_vertices

        # Cells of the glued tiles
        new_cells = np.zeros_like(cells)
        new_cells.ravel()[:] = translate[cells.flatten()]

        cells = np.vstack([cells, new_cells])
        # Update the periodicty mapping - slaves are new
        slave_vertices = translate[slave_vertices]
        # For the directions that do not evolve we add the periodic pairs
        for vm in vertex_mappings:
            vm.update(dict(izip(translate[vm.keys()], translate[vm.values()])))
        # Add the entities defined in terms of the vertices
        if mesh_data is not None:
            evolve_data(mesh_data, translate)
            
        # Iterate
        refine /= 2
        shift_x *= 2
    # Discard data not needed in next evolution
    return x, cells, shape[:-1]


def evolve_data(data, mapping):
    '''
    If mapping holds (tdim, tag) -> [tuple of indices]) where indices are 
    w.r.t of old numbering and mapping is old to new we simply add the mapped 
    entities.
    '''
    for key in data.keys():
        old = data[key]
            
        new = np.zeros_like(old)
        new.ravel()[:] = mapping[old.flatten()]
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
            indices = mesh_data[key]
            # These entity indices get the 'color'
            fill(mesh, indices.flatten(), tdim, key[1], f)
        containers[tdim] = f

    return containers


def as_meshf(mvc, init_value=0):
    '''Make a mesh function out of mesh value collection'''
    if isinstance(mvc, (tuple, list)):
        return [as_meshf(x, init_value) for x in mvc]

    if isinstance(mvc, dict):
        print as_meshf(mvc.values())
        return dict(zip(mvc.keys(), as_meshf(mvc.values())))

    # Base case
    mesh_f = MeshFunction('size_t', mvc.mesh(), mvc.dim(), init_value)
    fill_mf_from_mvc(mvc, mesh_f)

    return mesh_f

# --------------------------------------------------------------------

if __name__ == '__main__':
    from dolfin import (File, HDF5File, mpi_comm_world, Mesh, UnitSquareMesh,
                        Timer, info, SubsetIterator)

    if False:
        tile = UnitSquareMesh(1, 1)
    else:
        mesh_file = 'tile_2x2.h5'  #'cell_grid_2d.h5'
    
        comm = mpi_comm_world()
        h5 = HDF5File(comm, mesh_file, 'r')
        tile = Mesh()
        h5.read(tile, 'mesh', False)

    for n in (2, ):
        cell_dim = tile.topology().dim()
        facet_dim = cell_dim - 1

        surfaces = MeshFunction('size_t', tile, facet_dim, 0)
        h5.read(surfaces, 'facet')
        # Encode the data for evolve_data
        tile.init(facet_dim, 0)
        f2v = tile.topology()(facet_dim, 0)
        # Only want to evolve tag 1 (interfaces) for the facets. 
        facet_data = np.array([f2v(f.index()) for f in SubsetIterator(surfaces, 1)],
                              dtype='uintp')

        volumes = MeshFunction('size_t', tile, cell_dim, 0)
        h5.read(volumes, 'physical')
        # Encode the data for evolve_data
        tile.init(cell_dim, 0)
        c2v = tile.topology()(cell_dim, 0)
        # Only want to evolve tag 1
        cell_data = np.array([c2v(f.index()) for f in SubsetIterator(volumes, 1)],
                             dtype='uintp')

        data = {(facet_dim, 1): facet_data,
                (cell_dim, 1): cell_data}

        t = Timer('x')
        mesh, mesh_data = TileMesh(tile, (n, n), mesh_data=data)
        info('\tTiling took %g s\n' % t.stop())

        mvcs = mvc_from_data(mesh, mesh_data)
        foos = as_meshf(mvcs)
        
    # File('test.pvd') << mesh
    File('test_facet_marker.pvd') << foos[facet_dim]
    # File('test_cell_marker.pvd') << foos[cell_dim]


    