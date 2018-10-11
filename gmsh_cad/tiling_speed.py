from tiling import TileMesh, mf_from_data
import numpy as np
from dolfin import *


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
# shape = (2, 1)

# mf = MeshFunction('size_t', tile, tile.topology().dim()-1, 0)
# CompiledSubDomain('near(x[0], 0.5) || near(x[1], 0.5)').mark(mf, 1)

# mesh_data = {}
# mesh_data = load_data(mf, 1, mesh_data)
# fprint mesh_data

# mesh, mesh_data = TileMesh(tile, shape, mesh_data)
# # #print mesh_data

# fs = mf_from_data(mesh, mesh_data)
# # print mesh_data
# # print fs
# # print mesh_data.keys()
# # print fs[mesh_data.keys()[0]]

# File('x.pvd') << fs[1]
# # # # ee132e
ns = [128, 256, 512, 1024, 2048]#, 4096]
dts = []
for n in ns:
    shape = (n, n)
    
    t = Timer('s')
    mesh, mesh_data = TileMesh(tile, shape, mesh_data={})
    dts.append(t.stop())
    print mesh.num_cells()

File('test.pvd') << mesh
import matplotlib.pyplot as plt
import numpy as np

print np.polyfit(np.log(ns), np.log(dts), 1)
plt.figure()
plt.loglog(ns, dts, basex=2., basey=2., marker='x')
plt.show()
