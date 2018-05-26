from dolfin import Mesh, MeshFunction, HDF5File, MeshValueCollection
from msh_convert_cpp import fill_mvc_from_mf
import subprocess, os


def convert(msh_file, h5_file, save_mvc=True):
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

    # Save ALL data as facet_functions
    if not save_mvc:
        for region in ('facet_region.xml', 'physical_region.xml'):
            name, _ = region.split('_')
            r_xml_file = '_'.join([root, region])

            f = MeshFunction('size_t', mesh, r_xml_file)
            print name, f.array()
            out.write(f, name)

    for region in ('facet_region.xml', 'physical_region.xml'):
        name, _ = region.split('_')
        r_xml_file = '_'.join([root, region])

        f = MeshFunction('size_t', mesh, r_xml_file)
        # With mesh value collection we only store nonzero tags
        mvc = MeshValueCollection('size_t', mesh, f.dim())
        # Fill
        fill_mvc_from_mf(f, mvc)
        # And save
        out.write(f, name)
                    
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

    parser.add_argument('--save', type=int, help='save as pvd', default=0)
  
    parser.add_argument('--cleanup', type=str, nargs='+',
                        help='extensions to delete', default=('.xml'))
    args = parser.parse_args()

    # Protecting self
    assert not(set(('geo', '.geo')) & set(args.cleanup))
    
    try:
        msh_file, h5_file = args.io[:2]
    except ValueError:
        msh_file = args.io[0]

        root, ext = os.path.splitext(msh_file)
        h5_file = '.'.join([root, 'h5'])

    assert convert(msh_file, h5_file)

    h5 = HDF5File(mpi_comm_world(), h5_file, 'r')
    mesh = Mesh()
    h5.read(mesh, 'mesh', False)

    surfaces = MeshFunction('size_t', mesh, mesh.topology().dim()-1, 0)
    h5.read(surfaces, 'facet')

    volumes = MeshFunction('size_t', mesh, mesh.topology().dim(), 0)
    h5.read(volumes, 'physical')

    if args.save:
        File('results/%s_surf.pvd' % root) << surfaces
        File('results/%s_vols.pvd' % root) << volumes    

    cleanup(exts=args.cleanup)
