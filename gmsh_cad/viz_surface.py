from dolfin import SubsetIterator


def vtk_surface(surfaces, tag, output, value):
    assert surfaces.dim() == 2
    
    mesh = surfaces.mesh()
    mesh.init(2, 0)
    f2v = mesh.topology()(2, 0)
    
    ncells = sum(1 for _ in SubsetIterator(surfaces, tag))
    nvertices = 3*ncells

    x = mesh.coordinates()

    with open(output, 'w') as f:

        f.write("<?xml version=\"1.0\"?>\n")
        f.write("<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n")
        f.write("  <UnstructuredGrid>\n")

        f.write("    <Piece NumberOfPoints=\"%i\" NumberOfCells=\"%i\">\n" % (nvertices, ncells))

        f.write("      <Points>\n")
        f.write("        <DataArray type=\"Float32\" Name=\"vertices\" NumberOfComponents=\"3\" format=\"ascii\">\n")

        for e in SubsetIterator(surfaces, tag):
            coords = x[f2v(e.index())]
            for v in coords:
                f.write("%g %g %g " % tuple(v))

        f.write("        </DataArray>\n")
        f.write("      </Points>\n")
        f.write("      <Cells>\n")


        f.write("        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n")
        f.write(" ".join(map(str, range(nvertices))))
        f.write("        </DataArray>\n")

        
        f.write("        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n")
        f.write(' '.join(map(str, range(0, 3*ncells, 3))))
        f.write("        </DataArray>\n")

        f.write("        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n")
        f.write('5 '*ncells)
        f.write('        </DataArray>\n')

        f.write("      </Cells>\n")

        value = '%d ' % value
        print value
        f.write('      <CellData Scalars="tags">\n')
        f.write("        <DataArray type=\"Float32\" Name=\"tags\" format=\"ascii\">\n")
        f.write(value*ncells)
        f.write("        </DataArray>\n")
        f.write("      </CellData>\n")

        f.write("    </Piece>\n")
        f.write("  </UnstructuredGrid>\n")
        f.write("</VTKFile>\n")

# -----------------------------------------------------------------------------

if __name__ == '__main__':
    from dolfin import *
    import os

    n = 128
    mesh_file = 'tile_1_narrow_%d_%d.h5' % (n, n)

    comm = mpi_comm_world()
    h5 = HDF5File(comm, mesh_file, 'r')
    mesh = Mesh()
    h5.read(mesh, 'mesh', False)

    surfaces = MeshFunction('size_t', mesh, mesh.topology().dim()-1, 0)
    h5.read(surfaces, 'surfaces')

    rank = MPI.rank(comm)
    dir = './results%d' % n
    if rank == 0 and not os.path.exists(dir):
        os.mkdir(dir)
    MPI.barrier(comm)

    output = '%s/surf_%s_p%d.vtu' % (dir, mesh_file, rank)
    tag = 1

    vtk_surface(surfaces, tag, output, tag+rank)
