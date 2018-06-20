from dolfin import Mesh, MeshFunction, HDF5File, MeshValueCollection, SubsetIterator
from msh_convert_cpp import fill_mvc_from_mf
import subprocess, os


def convert(msh_file, h5_file, save_mvc=False):
    '''Temporary version of convertin from msh to h5'''
    root, _ = os.path.splitext(msh_file)
    assert os.path.splitext(msh_file)[1] == '.msh'
    assert os.path.splitext(h5_file)[1] == '.h5'

    # Get the xml mesh
    xml_file = '.'.join([root, 'xml'])
    subprocess.call(['dolfin-convert %s %s' % (msh_file, xml_file)], shell=True)
    # Success?
    assert os.path.exists(xml_file)

    mesh = Mesh(xml_file)
    out = HDF5File(mesh.mpi_comm(), h5_file, 'w')
    out.write(mesh, 'mesh')

    print('Mesh has %d cells' % mesh.num_cells())
    print('Mesh size %g %g' % (mesh.hmin(), mesh.hmax()))
    
    # Save ALL data as facet_functions
    names = ('surfaces', 'volumes')
    if not save_mvc:
        for name, region in zip(names, ('facet_region.xml', 'physical_region.xml')):
            r_xml_file = '_'.join([root, region])

            f = MeshFunction('size_t', mesh, r_xml_file)
            print('%d %s with 1' % (sum(1 for _ in SubsetIterator(f, 1)), name))
            out.write(f, name)

        return True

    for name, region in zip(names, ('facet_region.xml', 'physical_region.xml')):
        r_xml_file = '_'.join([root, region])

        f = MeshFunction('size_t', mesh, r_xml_file)
        # With mesh value collection we only store nonzero tags
        mvc = MeshValueCollection('size_t', mesh, f.dim())
        # Fill
        fill_mvc_from_mf(f, mvc)
        # And save
        out.write(mvc, name)
                    
    return True
    

def cleanup(files=None, exts=()):
    '''Get rid of xml'''
    if files is not None:
        return map(os.remove, files)
    else:
        files = filter(lambda f: any(map(f.endswith, exts)), os.listdir('.'))
        print('Removing', files)
        return cleanup(files)
                    
# --------------------------------------------------------------------

if __name__ == '__main__':
    from dolfin import File, mpi_comm_world
    import argparse


    parser = argparse.ArgumentParser(description='Convert msh file to h5')
    parser.add_argument('io', type=str, nargs='+', help='input [output]')
    parser.add_argument('--cleanup', type=str, nargs='+',
                        help='extensions to delete', default=('.xml'))

    # Save the nonzero markers as mesh value collections
    save_mvc_parser = parser.add_mutually_exclusive_group(required=False)
    save_mvc_parser.add_argument('--save_mvc', dest='save_mvc', action='store_true')
    save_mvc_parser.add_argument('--no_save_mvc', dest='save_mvc', action='store_false')
    parser.set_defaults(save_mvc=False)

    # Save the mesh markers for visualizing
    save_pvd_parser = parser.add_mutually_exclusive_group(required=False)
    save_pvd_parser.add_argument('--save_pvd', dest='save_pvd', action='store_true')
    save_pvd_parser.add_argument('--no_save_pvd', dest='save_pvd', action='store_false')
    parser.set_defaults(save_pvd=False)

    args = parser.parse_args()

    # Protecting self
    assert not(set(('geo', '.geo')) & set(args.cleanup))
    
    try:
        msh_file, h5_file = args.io[:2]
    except ValueError:
        msh_file = args.io[0]

        root, ext = os.path.splitext(msh_file)
        h5_file = '.'.join([root, 'h5'])

    assert convert(msh_file, h5_file, args.save_mvc)

    if args.save_pvd:
        h5 = HDF5File(mpi_comm_world(), h5_file, 'r')
        mesh = Mesh()
        h5.read(mesh, 'mesh', False)

        surfaces = MeshFunction('size_t', mesh, mesh.topology().dim()-1, 0)
        volumes = MeshFunction('size_t', mesh, mesh.topology().dim(), 0)

        if not args.save_mvc:
            h5.read(surfaces, 'surfaces')
            h5.read(volumes, 'volumes')
        # The data is mesh value collections
        else:
            from tiling_cpp import fill_mf_from_mvc

            surfaces_mvc = MeshValueCollection('size_t', mesh, mesh.topology().dim()-1)
            h5.read(surfaces_mvc, 'surfaces')
            fill_mf_from_mvc(surfaces_mvc, surfaces)

            volumes_mvc = MeshValueCollection('size_t', mesh, mesh.topology().dim())
            h5.read(volumes_mvc, 'volumes')
            fill_mf_from_mvc(volumes_mvc, volumes)

        File('results/%s_surf.pvd' % root) << surfaces
        File('results/%s_vols.pvd' % root) << volumes    

    cleanup(exts=args.cleanup)
