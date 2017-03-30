<p align="center">
  <img src="https://github.com/MiroK/emi-cylinders/blob/master/doc/grid-iso.png">
</p>

Scripts for generating geometries (and meshes) for sheets of cells. Sheet here is a
a grid or array of connected cells. The main workhorse is `GMSH`. The Python scripts are
here simply to make adjusting the domain parameters more convenient and also to convert
the mesh with `dolfin-convert` to `FEniCS` readable format. The first and the following 
picture are created using `grid-iso.py`. As you can see the cells are now connected *only*
over the curved surfaces and the resulting mesh can model some isotropic sheet.

<p align="center">
  <img src="https://github.com/MiroK/emi-cylinders/blob/master/doc/grid-iso-large.png">
</p>

On the other hand, `grid.py` gives you cells that are connected over the flat cylinder tops
*and* the curved wall. The mesh is therefore more suitable for some anisotropic sheets.

<p align="center">
  <img src="https://github.com/MiroK/emi-cylinders/blob/master/doc/grid.png">
</p>
