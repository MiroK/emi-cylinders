from __future__ import print_function
from dolfin import Cell, Point
from mpi4py import MPI as pyMPI
import numpy as np


class PointProbe(object):
    '''Perform efficient evaluation of function u at fixed points'''
    def __init__(self, u, locations):
        # The idea here is that u(x) means: search for cell containing x,
        # evaluate the basis functions of that element at x, restrict
        # the coef vector of u to the cell. Of these 3 steps the first
        # two don't change. So we cache them

        # Locate each point
        mesh = u.function_space().mesh()
        limit = mesh.num_entities(mesh.topology().dim())
        bbox_tree = mesh.bounding_box_tree()
        # In parallel x might not be on process, the cell is None then
        cells_for_x = [None]*len(locations)
        for i, x in enumerate(locations):
            cell = bbox_tree.compute_first_entity_collision(Point(*x))
            from dolfin import info
            if -1 < cell < limit:
                cells_for_x[i] = cell

        V = u.function_space()
        element = V.dolfin_element()

        size = V.ufl_element().value_size()
        # Build the sampling matrix
        evals = []
        dm = V.dofmap()
        for x, cell in zip(locations, cells_for_x):
            # If we own the cell we alloc stuff and precompute basis matrix
            if cell is not None:
                basis_matrix = np.zeros(size*element.space_dimension())
                coefficients = np.zeros(element.space_dimension())
                # NOTE: avoid using DOLFIN's restric; instead reach into
                # function's vector
                cell_dofs = dm.cell_dofs(cell) + dm.ownership_range()[0]

                cell = Cell(mesh, cell)
                vertex_coords, orientation = cell.get_vertex_coordinates(), cell.orientation()
                # Eval the basis once
                element.evaluate_basis_all(basis_matrix, x, vertex_coords, orientation)

                basis_matrix = basis_matrix.reshape((element.space_dimension(), size)).T
                # Make sure foo is bound to right objections
                def foo(u_vec, c=coefficients, A=basis_matrix, dofs=cell_dofs):
                    # Restrict for each call using the bound cell, vc ...
                    c[:] = u_vec.getValues(dofs)
                    # A here is bound to the right basis_matrix
                    return np.dot(A, c)
            # Otherwise we use the value which plays nicely with MIN reduction
            else:
                foo = lambda u, size=size: (np.finfo(float).max)*np.ones(size)

            evals.append(foo)

        self.probes = evals
        # To get the correct sampling on all cpus we reduce the local samples across
        # cpus
        self.comm = pyMPI.COMM_WORLD
        self.readings = np.zeros(size*len(locations), dtype=float)
        self.readings_local = np.zeros_like(self.readings)
        # Return the value in the shape of vector/matrix
        self.nprobes = len(locations)

    def sample(self, u):
        '''Evaluate the probes listing the time as t'''
        u_vec = as_backend_type(u.vector()).vec()  # This is PETSc
        self.readings_local[:] = np.hstack([f(u_vec) for f in self.probes])    # Get local
        self.comm.Reduce(self.readings_local, self.readings, op=pyMPI.MIN)  # Sync

        return self.readings.reshape((self.nprobes, -1))

# --------------------------------------------------------------------

if __name__ == '__main__':
    from dolfin import *

    mesh = UnitSquareMesh(256, 256)

    theta = np.linspace(0, 2*np.pi, 10)
    pts = np.c_[0.5 + 0.25*np.sin(theta), 0.5 + 0.25*np.cos(theta)]

    V = FunctionSpace(mesh, 'CG', 1)
    f = Expression('t*(x[0] + 2*x[1])', degree=1, t=1)
    
    v = interpolate(f, V)
    probes = PointProbe(v, pts)

    table = pts
    for i in range(10):
        f.t += i*0.1
        v.assign(interpolate(f, V))
        
        probe_values = probes.sample(v)
        if probes.comm.rank == 0:
            table = np.c_[table, probe_values.flatten()]

    if probes.comm.rank == 0:
        x, y = pts.T

        t = 1
        error = 0
        for i in range(10):
            t += i*0.1
            v = t*(x + 2*y)
            error = max(error, np.linalg.norm(v - table[:, 2+i]))
        print(error)
