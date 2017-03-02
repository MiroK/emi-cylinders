import os, subprocess
from dolfin import Mesh, FacetFunction, CellFunction, File, HDF5File

        
def geo_from_template(template_file, specs, out_file):
    '''
    Given template geo file we update its parameters from specs and produce
    an out_file which can be used to generate meshes with gmsh.
    '''
    assert os.path.exists(template_file)
    root, ext = os.path.splitext(out_file)
    assert ext == '.geo'

    # Figure out where to cut
    with open(template_file) as f:
        content = f.readlines()
        for i, line in enumerate(content):
            if line.startswith('//! Cut'):
                break
            cut_line = i
    header, template = content[:cut_line], content[cut_line:]

    specs_keys = set(specs.keys())
    found_keys = set()
    # Now peform substitutions 
    new_header = []
    for line in header:
        if '=' in line:
            key, value = line.split('=')
            key = key.strip()
            found_keys.add(key)

            new_value = specs[key]
            new_line = '%s = %g;\n' % (key, new_value)
        else:
            new_line = line
        new_header.append(new_line)
    assert specs_keys == found_keys, '%s !+ %s' % (specs_keys, found_keys)

    content = new_header + template
    with open(out_file, 'w') as f:
        for line in content: f.write(line)

    return root
    


def generate_mesh(root, scale=1):
    '''
    Generate mesh from root.geo file using gmsh and dolfin-convert. Optionally
    we can scale the sizes with scale.
    '''
    assert os.path.splitext(root)[1] == ''

    geo_file = '%s.geo' % root
    msh_file = '%s.msh' % root

    # Make sure we have file
    assert os.path.exists(geo_file), 'Missing geo file %s .' % geo_file

    # Make the gmsh file for current size
    subprocess.call(['gmsh -clscale %g -3 -optimize %s' % (scale, geo_file)], shell=True)
    assert os.path.exists(msh_file)

    # Convert to xdmf
    xml_file = '%s.xml' % root              
    xml_facets = '%s_facet_region.xml' % root
    xml_volumes = '%s_physical_region.xml' % root
    
    subprocess.call(['dolfin-convert %s %s' % (msh_file, xml_file)], shell=True)
    # All 3 xml files should exist
    assert all(os.path.exists(f) for f in (xml_file, xml_facets, xml_volumes))
    
    # Convert
    h5_file = '%s.h5' % root
    cmd = r'''python -c"from dolfin import Mesh, HDF5File, MeshFunction;\
mesh=Mesh('%s');\
facet_f=MeshFunction('size_t', mesh, '%s');\
cell_f=MeshFunction('size_t', mesh, '%s');
out=HDF5File(mesh.mpi_comm(), '%s', 'w');\
out.write(mesh, '/mesh');\
out.write(facet_f, '/facet_markers');\
out.write(cell_f, '/cell_markers');\
"''' % (xml_file, xml_facets, xml_volumes, h5_file)
        
    subprocess.call([cmd], shell=True)
    # Success?
    assert os.path.exists(h5_file)

    # Cleanup
    [os.remove(f) for f in (xml_file, xml_facets, xml_volumes)]
    os.remove(msh_file)
    
    return h5_file


def as_pvd(h5_file):
    '''Store facet and cell function for pvd'''
    root, ext = os.path.splitext(h5_file)

    mesh = Mesh()
    hdf = HDF5File(mesh.mpi_comm(), h5_file, 'r')
    hdf.read(mesh, '/mesh', False)
    
    facet_markers = FacetFunction('size_t', mesh)
    hdf.read(facet_markers, '/facet_markers')

    cell_markers = CellFunction('size_t', mesh)
    hdf.read(cell_markers, '/cell_markers')

    File(root + 'facets' + '.pvd') << facet_markers
    File(root + 'volumes' + '.pvd') << cell_markers

    return True
