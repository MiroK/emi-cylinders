from evtk.hl import unstructuredGridToVTK
from evtk.vtk import VtkTriangle, VtkGroup


pvtu_code = '''<?xml version="1.0"?>
<VTKFile type="PUnstructuredGrid" version="0.1">
<PUnstructuredGrid GhostLevel="0">
<PPoints>
<PDataArray type="Float64" NumberOfComponents="3" />
</PPoints>
<PCellData>
<PDataArray type="UInt32" Name="connectivity" />
<PDataArray type="UInt32" Name="offsets" />
<PDataArray type="UInt8" Name="types" />
</PCellData>
<PCellData Scalars="%(f)s">
<PDataArray type="Float64" Name="%(f)s" NumberOfComponents="0" />
</PCellData>
%(pieces)s
</PUnstructuredGrid>
</VTKFile>
'''


def facet_to_DLTdofs(V, facets):
    '''Facet -> degree of freedom of DLT space'''
    assert V.ufl_element().family() == 'HDiv Trace', V.ufl_element().family()
    assert V.ufl_element().degree() == 0
    assert V.ufl_element().value_shape() == ()

    dm = V.dofmap()

    mesh = V.mesh()
    fdim, cdim = mesh.topology().dim()-1, mesh.topology().dim()
    mesh.init(fdim, cdim)
    mesh.init(cdim, fdim)
    f2c = mesh.topology()(fdim, cdim)
    c2f = mesh.topology()(cdim, fdim)
        
    mapping = []
    for lid, fid in enumerate(facets):
        # Get dof as the one at which the cells agree
        cs = f2c(fid)
        dofs = list(map(dm.cell_dofs, cs))
        # A boundary facet
        if len(cs) == 1:
            c0, = cs
            dofs, = dofs
            dof = dofs[list(c2f(c0)).index(fid)]
        else:
            c0, c1 = cs
            dof,  = (set(dofs[0]) & set(dofs[1]))
        mapping.append(dof)

    mapping = np.array(mapping)
    # Chop it to local
    first, last = dm.ownership_range()
    offset = mapping + first

    mask,  = np.where(np.logical_and(offset >= first, offset < last))

    return mask, mapping[mask]


def triangulation(mesh, facets):
    '''Triangulation of facets of the mesh
    
    idx   = such that mesh.coordinates[idx] are vertices of the triangulation
    cells = cells of the triangulation defined in terms of idx (so local)
    '''
    # NOTE: these are per CPU calculations so I don't care about parallel
    fdim = mesh.topology().dim()-1
    mesh.init(fdim, 0)
    f2v = mesh.topology()(fdim, 0)

    vertex_mapping, cells = {}, []
    local_id = 0
    # We build as global to local
    for facet in facets:
        as_vertices = f2v(facet)
        cell_local = []
        for v in as_vertices:
            if v not in vertex_mapping:
                vertex_mapping[v] = local_id
                cell_local.append(local_id)
                local_id += 1
            else:
                cell_local.append(vertex_mapping[v])
        # Collect the cells in local
        cells.append(cell_local)
    cells = np.array(cells)
                    
    # Since keys in vertex_mapping are basically a range
    global_v, local_v = map(np.array, zip(*vertex_mapping.items()))
    # Now local are a range
    vertex_mapping = global_v[np.argsort(local_v)]
    # So (facet mesh local vertex index) -> (cpu local full mesh vertex index)
    return vertex_mapping, cells
    

class DltWriter(object):
    '''
    VTK writer for scalar order 0 discontinuous Lagrange trace functions.
    '''
    def __init__(self, path, u, active_facets=None):
        '''Consider only on active_facets'''
        V = u.function_space()
        mesh = V.mesh()
        fdim = mesh.topology().dim() - 1

        # All of them
        if active_facets is None: active_facets = range(mesh.num_entities(fdim))

        # Now how woud we do the update
        mask, mapping = facet_to_DLTdofs(V, active_facets)

        values = np.zeros(len(mapping), dtype=float)
        cell_data = {u.name(): values}  # We will update these

        
        active_facets = active_facets[mask]        
        # Let's compute things for the grid
        vertices, cells = triangulation(mesh, facets=active_facets)

        x, y, z = map(np.array, mesh.coordinates()[vertices].T)
        ncells, nvertices_cell = cells.shape
        assert nvertices_cell == 3  # FIXME: only triangles
        
        connectivity = cells.flatten()
        offsets = np.arange(nvertices_cell,
                            len(connectivity) + nvertices_cell,
                            nvertices_cell)

        cell_types = VtkTriangle.tid*np.ones(ncells)
        # The grid is made of x, y, z, connectivity, offsets, cell_types
        # And cell_data dictionary: alloc

        self.write_vtu_piece = lambda path, rank, counter: (
            np.copyto(values, u.vector().get_local()[mapping]),
            #
            unstructuredGridToVTK("%s_p%d_%06d" % (path, rank, counter),
                                  x, y, z,
                                  connectivity,
                                  offsets,
                                  cell_types,
                                  cellData=cell_data),
            #
            "%s_p%d_%06d.vtu" % (path, rank, counter)
        )[-1]

        comm = mesh.mpi_comm()
        # Only other piece to remember is for parallel writeing on root
        self.counter = 0  # Of times write_vtu_piece was called
        self.path = path
        self.world_rank = comm.rank
        self.world_size = comm.size
        self.u_name = u.name()

        self.world_rank == 0 and setattr(self, 'group', VtkGroup(path))

    def __enter__(self):
        return self
    
    def write(self, t):
        '''Write a new piece'''
        # Each process writes the vtu file
        group_path = self.write_vtu_piece(self.path, self.world_rank, self.counter)

        # Root write one pvtu file
        if self.world_size > 1:
            group_path = '%s%06d.pvtu' % (self.path, self.counter)

            if self.world_rank == 0:
                with open(group_path, 'w') as group_file:
                    group_file.write(pvtu_code %
                        {'f': self.u_name,
                         'pieces': '\n'.join(['<Piece Source="%s_p%d_%06d.vtu" />' % (self.path, rank, self.counter)
                                              for rank in range(self.world_size)])
}
)

        self.counter += 1  # Of times write_vtu_piece was called
                                     
        # Update group file
        self.world_rank == 0 and self.group.addFile(filepath=group_path, sim_time=t)
            
    def __exit__(self, exc_type, exc_value, exc_traceback):
        '''Complete the pvd file'''
        self.world_rank == 0 and self.group.save()
                                     
# --------------------------------------------------------------------


if __name__ == '__main__':
    from dolfin import *
    import numpy as np
    
    mesh = UnitCubeMesh(16, 16, 16)
    gdim = mesh.geometry().dim()
    fdim = mesh.topology().dim() - 1
        
    V = FunctionSpace(mesh, 'Discontinuous Lagrange Trace', 0)
    g = interpolate(Expression('x[0]+x[1]', degree=1), V)

    facet_f = MeshFunction('size_t', mesh, fdim, 0)
    CompiledSubDomain('near(x[0]*(1-x[0]), 0) || near(x[1]*(1-x[1]), 0) || near(x[0], x[1])').mark(facet_f, 1)
    
    active_facets = np.array([f.index() for f in SubsetIterator(facet_f, 1)])
    
    dlt_w = DltWriter(path='test_DLT', u=g, active_facets=active_facets)
    with dlt_w as f:
        f.write(0)

        g.assign(interpolate(Constant(2), V))
        f.write(0.1)

    mesh_file = 'tile_1_hein_GMSH307_20_4.h5'
    
    # Test with mesh for EMI
    comm = MPI.comm_world
    h5 = HDF5File(comm, mesh_file, 'r')
    mesh = Mesh()
    h5.read(mesh, 'mesh', False)  

    surfaces = MeshFunction('size_t', mesh, mesh.topology().dim()-1, 0)
    h5.read(surfaces, 'surfaces')

    V = FunctionSpace(mesh, 'Discontinuous Lagrange Trace', 0)
    g = interpolate(Expression('x[0]+x[1]', degree=1), V)

    active_facets = np.array([f.index() for f in SubsetIterator(surfaces, 1)])
    print(V.dim(), comm.allreduce(len(active_facets)))
    dlt_w = DltWriter(path='test_DLT', u=g, active_facets=active_facets)
    with dlt_w as f:
        f.write(0)
