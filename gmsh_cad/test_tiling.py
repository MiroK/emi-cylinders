from dolfin import *
from tiling import TileMesh, mf_from_data, load_data
from mpi4py import MPI
import numpy as np



def test(path, type='mf'):
    '''Evolve the tile in (n, n) pattern checking volume/surface properties'''

    comm = MPI.COMM_WORLD
    h5 = HDF5File(comm, path, 'r')
    tile = Mesh()
    h5.read(tile, 'mesh', False)

    init_container = lambda type, dim: (MeshFunction('size_t', tile, dim, 0)
                                        if type == 'mf' else
                                        MeshValueCollection('size_t', tile, dim))
        
    for n in (2, 4):
        data = {}
        checks = {}
        for dim, name in zip((2, 3), ('surfaces', 'volumes')):
            # Get the collection
            collection = init_container(type, dim)
            h5.read(collection, name)
                
            if type == 'mvc': collection = as_meshf(collection)
            
            load_data(tile, h5, name, dim, data)

            if dim == 2:
                # Interface area area
                check = (lambda m, f: assemble(FacetArea(m)*ds(domain=m, subdomain_data=f, subdomain_id=1)+
                                               avg(FacetArea(m))*dS(domain=m, subdomain_data=f, subdomain_id=1)), )
            else:
                check = (lambda m, f: assemble(CellVolume(m)*dx(domain=m, subdomain_data=f, subdomain_id=0)),
                         lambda m, f: assemble(CellVolume(m)*dx(domain=m, subdomain_data=f, subdomain_id=1)))

            checks[dim] = [lambda m, f, t=tile, c=collection, n=n, check=c: (abs(check(m, f)-n**2*check(t, c))/(n**2*check(t, c)), check(t, c))
                           for c in check]
                
        t = Timer('x')
        mesh, mesh_data = TileMesh(tile, (n, n), mesh_data=data)
        info('\tTiling took %g s. Ncells %d, nvertices %d, \n' % (t.stop(), mesh.num_vertices(), mesh.num_cells()))
            
        foos = mf_from_data(mesh, mesh_data)
        # Mesh Functions
        from_mf = np.array([check(mesh, foos[dim]) for dim in (2, 3) for check in checks[dim]])
        
        # Ignote this as mvc support dropped
        # mvcs = mvc_from_data(mesh, mesh_data)
        # foos = as_meshf(mvcs)
        # # Mesh ValueCollections
        # from_mvc = np.array([checks[dim](mesh, foos[dim]) for dim in (2, 3)])

        # assert np.linalg.norm(from_mf - from_mvc) < 1E-13
        # I ignore shared facets so there is bound to be some error in facets
        # Volume should match well
        print(from_mf)
    
# -----------------------------------------------------------------------------

if __name__ == '__main__':
    # test('tile_2x2.h5')

    test('tile_1_narrow.h5')
