from driver import geo_from_template, generate_mesh, as_pvd
from dolfin import Mesh, HDF5File

specs = {'n_cylinders': 4,
         # The first cylinder has center at x0, y0, z0
         'x0': 0,
         'y0': 0,
         'z0': 0,
         # The height of each cylinder is H with char mesh size size_c
         'R': 0.5,
         'H': 1.5,
         'size_R': 0.2,
         # The height of the joint regions is h/2 with char mesh size size_r
         'r': 0.25,
         'h': 0.5,
         'size_r': 0.2,
         # The cylinders are enclosed in a bbox which leaves dx, dy gaps around the
         # square that bounds the cylinder crossection. In z direction the gap is dz.
         # For bbox the char size is size_b
         'dx': 0.2,
         'dy': 0.2,
         'dz': 0.2,
          'size_b': 0.3
         # Each cylinder volume gets number 1, ..., n_cylinders. These markers are also
         # inherited by surfaces that bound the cylinder except the surface that is
         # shared by 2 cylinders. The label for the cylinder is (volume label +
         # n_cylinders). The outer volume, i.e. bounding box - cylinder union is tagged
         }

assert specs['n_cylinders'] > 1
assert all(specs[key] > 0 for key in specs if key not in ('n_cylinders', 'x0', 'y0', 'z0'))
assert specs['R'] > specs['r']
assert specs['H'] > specs['h']

root = geo_from_template('array.geo', specs, 'array_test.geo')
mesh_file = generate_mesh(root, scale=1)
# Dump to paraview for inspection
assert as_pvd(mesh_file)
