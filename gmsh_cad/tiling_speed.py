from tiling import TileMesh
from dolfin import UnitSquareMesh, Timer


tile = UnitSquareMesh(1, 1)

ns = [128, 256, 1024, 2048, 4096]
dts = []
for n in ns:
    shape = (n, )*2
    
    t = Timer('s')
    mesh, mesh_data = TileMesh(tile, shape, mesh_data={})
    dts.append(t.stop())
    print mesh.num_cells()
    
import matplotlib.pyplot as plt
import numpy as np

print np.polyfit(np.log(ns), np.log(dts), 1)
plt.figure()
plt.loglog(ns, dts, basex=2., basey=2.)
plt.show()
