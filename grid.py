from driver import geo_from_template, generate_mesh, as_pvd
from dolfin import Mesh, HDF5File

# NOTE: here the cylinders look like AAA batteries which are coupled in z and y
# directions to create a sheate. The coupling is respectively thought top and
# side surface which leads to a natural anisotropy in the system. For a sheat of
# 'watch batteries' where the coupling is via sides/curved surfaces thus
# resulting in an isotropic sheat see grid-iso

specs = {'z_cylinders': 3,
         'y_cylinders': 4,
         # The height of each cylinder is H, radius is R with char mesh size size_R
         'R': 0.8,
         'H': 1.5,
         'size_R': 0.2,
         # The height of the joint regions in z-direction is h with char mesh size size_r
         'r': 0.5,
         'h': 0.5,
         'size_r': 0.2,
         # The length of the joint regions in y-direction is l with char mesh size size_l
         't': 0.25,
         'l': 0.5,
         'size_l': 0.2,
         # Note that to make the transition of the joint to the cylindrical surface
         # there is a weld which results in effective radius of q > t. Below q is the
         # scaling factor so that q[radius] := q*t
         'q': 1.1,
         # The cylinders are enclosed in a bbox which leaves dx, dy gaps around the
         # square that bounds the cylinder crossection. In z direction the gap is dz.
         # For bbox the char size is size_b
         'dx': 0.2,
         'dy': 0.2,
         'dz': 0.2,
         'size_b': 0.3
         # The cylinders are numbered (the volume) starting from on first growing in 
         # y-direction then in z-direction. The volume id determines also id of all the
         # surfaces that are not shared. For the shared surfaces the y-shared surface is
         # volume_id + (z_cylinders*y_cylinders). In z direction this is 
         # volume_id + 2*(z_cylinders*y_cylinders). Note that in this policy it is the first
         # cylinder which hits the surface that gets to name it.
         # The outer volume, i.e. bounding box - cylinder union is tagged as 0
         }

# NOTE: these are some sanity checks. It is not guaranteed that the geometry is
# well-defined even if all pases
assert specs['z_cylinders'] > 1 and specs['y_cylinders'] > 1
assert all(specs[key] > 0 for key in specs)
q, t, r, R = specs['q'], specs['t'], specs['r'], specs['R']
assert q > 1
q = t*q
assert R > r and R > q

root = geo_from_template('grid.geo', specs, 'grid_test.geo')
mesh_file = generate_mesh(root)
# Dump to paraview for inspection
assert as_pvd(mesh_file)
