from dolfin import CompiledSubDomain, Mesh, MeshEditor, MeshFunction
from test_periodic import compute_vertex_periodicity
import numpy as np
import operator


def TileMesh(tile, shape, mesh_data=None, TOL=1E-9):
    '''
    [tile tile;
     tile tile;
     tile tile;
     tile tile]

    The shape is an ntuple describing the number of pieces put next 
    to each other in the i-th axis. mesh_data : tdim -> (tag -> [entities]) 
    is the way to encode mesh data of the tile.
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
        x, cells, shape = evolve(x, cells, vertex_mappings, shifts_x, shape, mesh_data=mesh_data)

    # Done evolving, we can write data
    tdim = tile.topology().dim()
    ctype = str(tile.ufl_cell())

    return make_mesh(x, cells, ctype, tdim, gdim, mesh_data=mesh_data)

        
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
        new_cells.ravel()[:] = map(translate.__getitem__, cells.flat)  # FIXME: numpy access

        cells = np.r_[cells, new_cells]
        # Update the periodicty mapping - slaves are new
        slave_vertices = map(translate.__getitem__, slave_vertices)  # FIXME: numpy access
        # For the directions that do not evolve we add the periodic pairs
        for vm in vertex_mappings:
            vm.update(dict(zip(translate[vm.keys()], translate[vm.values()])))
        # Add the entities defined in terms of the vertices
        if mesh_data is not None: evolve_data(mesh_data, translate)
            
        # Iterate
        refine /= 2
        shift_x *= 2
    # Discard data not needed in next evolution
    return x, cells, shape[:-1]


def evolve_data(data, mapping):
    '''
    If mapping holds tdim -> (tag -> [tuple of indices]) where indices are 
    w.r.t of old numbering and mapping is old to new we simply add the mapped 
    entities.
    '''
    for tdim in data:
        data_tdim = data[tdim]

        for tag in data_tdim:
            old = data_tdim[tag]
            
            new = np.zeros_like(old)
            new.ravel()[:] = map(mapping.__getitem__, old.flat)
            data_tdim[tag] = np.r_[old, new]


def make_mesh(vertices, cells, ctype, tdim, gdim, mesh_data=None):
    '''Mesh by MeshEditor from vertices and cells'''
    mesh = Mesh()
    editor = MeshEditor()
    editor.open(mesh, ctype, tdim, gdim)

    editor.init_vertices(len(vertices))
    editor.init_cells(len(cells))

    for vi, x in enumerate(vertices): editor.add_vertex(vi, x)
    
    for ci, c in enumerate(cells): editor.add_cell(ci, *c)

    editor.close()
    # For now data we're done
    if mesh_data is None: return mesh

    # FIXME: do vertex lookup once here!
    # FIXME: should work mesh_data copy?

    mesh_functions = {}
    # We have define entities in terms of vertex numbering
    for tdim in mesh_data:
        # So we'll be getting the entity index by lookup
        mesh.init(tdim)
        mesh.init(0, tdim)
        v2entity = mesh.topology()(0, tdim)
        # And color corresponding mesh function
        f = MeshFunction('size_t', mesh, tdim, 0)
        f_values = f.array()
        for tag, entities in mesh_data[tdim].iteritems():
            # Our entity should be the single intersection of those connected
            # to its vertices
            entity_indices = [reduce(operator.and_, (set(v2entity(v)) for v in entity)).pop()
                              for entity in entities]
            f_values[entity_indices] = tag
        # Add
        mesh_functions[tdim] = f
        
    return mesh, mesh_functions

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

    for n in (2, 4, 8, 16):
        facet_dim = tile.topology().dim()-1
        surfaces = MeshFunction('size_t', tile, facet_dim, 0)
        h5.read(surfaces, 'facet')
        # Encode the data for evolve_data
        tile.init(facet_dim, 0)
        f2v = tile.topology()(facet_dim, 0)
        # Only want to evolve tag 1 (interfaces) for the facets. 
        facet_data = {1: np.array([f2v(f.index()) for f in SubsetIterator(surfaces, 1)])}
        data = {facet_dim: facet_data}

        t = Timer('x')
        mesh, foos = TileMesh(tile, (n, n), mesh_data=data)
        info('\nTiling took %g s' % t.stop())
        
    File('test.pvd') << mesh
    File('test_marker.pvd') << foos[facet_dim]


    
