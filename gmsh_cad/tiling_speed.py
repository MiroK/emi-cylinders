from tiling import TileMesh
from dolfin import UnitSquareMesh, Timer, File


tile = UnitSquareMesh(1, 1)
mesh, mesh_data = TileMesh(tile, shape=(13, 7), mesh_data={})
File('foo.pvd') << mesh

if False:
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

    a, b = np.polyfit(np.log(ns), np.log(dts), 1)
    print a, b

    ns, dts = map(np.array, (ns, dts))
    
    plt.figure()
    plt.loglog(ns, dts, basex=2., basey=2., marker='x')
    plt.loglog(ns, np.exp(b)*ns**a, basex=2., basey=2., linestyle='dashed')

    plt.show()
