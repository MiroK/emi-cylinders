from tiling import TileMesh, mf_from_data
from dolfin import *
import numpy as np


def load_data(mf, dim, data):
    '''
    Fill the data dictionary with data_set representing mesh function with 
    dim over mesh read from h5_file according to key spec expected by tiling 
    algorithm.
    '''
    mesh = mf.mesh()
    # Data to evolve
    mesh.init(dim, 0)
    e2v = tile.topology()(dim, 0)

    tags = set(mf.array())
    # Don't evolve zero - we initialize to it
    if 0 in tags: tags.remove(0)

    for tag in tags:
        data[(dim, tag)] = np.array([e2v(e.index()) for e in SubsetIterator(mf, tag)],
                                    dtype='uintp')
    return data


tile = UnitSquareMesh(2, 2)

mf = MeshFunction('size_t', tile, tile.topology().dim()-1, 0)
CompiledSubDomain('near(x[0], 0.5) || near(x[1], 0.5)').mark(mf, 1)

mesh_data = {}
mesh_data = load_data(mf, 1, mesh_data)

mesh, mesh_data = TileMesh(tile, shape=(23, 13), mesh_data=mesh_data)
fs = mf_from_data(mesh, mesh_data)

File('mesh_f.pvd') << fs[1]

tile = UnitSquareMesh(1, 1)
mesh, mesh_data = TileMesh(tile, shape=(13, 8), mesh_data={})

if True:
    ns = [128, 256, 1024, 2048, 4096]
    dts = []
    for n in ns:
        shape = (n+1, n-1)
    
        t = Timer('s')
        mesh, mesh_data = TileMesh(tile, shape, mesh_data={})
        dts.append(t.stop())
        print mesh.num_cells()
        
    import matplotlib.pyplot as plt
    import numpy as np

    a, b = np.polyfit(np.log(ns), np.log(dts), 1)
    print a, b

    ns, dts = map(np.array, (ns, dts))
    
    plt.figure()
    plt.loglog(ns, dts, basex=2., basey=2., marker='x')
    plt.loglog(ns, np.exp(b)*ns**a, basex=2., basey=2., linestyle='dashed')

    plt.show()
