from __future__ import print_function
from dolfin import Cell, Point, as_backend_type
from mpi4py import MPI as pyMPI
import numpy as np


def probe_cell_at(m, n, level, point):
    '''Sampling points for [m, n] Hein cell.

    Level specifies Z coordinate: the points are then as
    
    l=2 or -2
     ________________
    |_|            |_|
    |_|-1         1|_|
    |_|____________|_|
    
    l=1 or -1
      ________________
      |_|            |_|
    -2|_|-1         1|_|2
      |_|____________|_|

    l=0
      ________________
    -4|_|-3         4|_|3
    -2|_|-1         1|_|2
      |_|____________|_|
    '''
    # In base cell that touches [0, 0, 0]
    z = {0: 0,
         1: 8,
         -1: -8,
         2: 11.5,
         -2: -11.5}[level]

    (x, y) = {0: {1: (98, 8), 2: (100, 8), 3: (100, 15), 4: (98, 15),
                  -1: (2, 8), -2: (0, 8), -3: (2, 15), -4: (0, 15)},
              #
              1: {1: (98, 11.5), -1: (2, 11.5), 2: (100, 11.5), -2: (0, 11.5)},
              -1: {1: (98, 11.5), -1: (2, 11.5), 2: (100, 11.5), -2: (0, 11.5)},
              #
              2: {1: (98, 11.5), -1: (2, 11.5)},
              -2: {1: (98, 11.5), -1: (2, 11.5)},
    }[level][point]

    point = np.array([x, y, z])
    # Shift
    point += np.array([m*100, n*23, 0])
    # In mm
    point *= 1E-3

    return point


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
                basis_matrix.ravel()[:] = element.evaluate_basis_all(x, vertex_coords, orientation)
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

    
class ApproxPointProbe(object):
    '''Eval u at locations by samling the nearest dof.'''
    # So this makes sence for DOFs which are point evaluations
    def __init__(self, u, locations):
        assert u.ufl_shape == ()
        # Let's do the snap to nearest search
        V = u.function_space()
        mesh = V.mesh()
        dofs_x = V.tabulate_dof_coordinates().reshape((-1, mesh.geometry().dim()))

        # The idea is to get dofs closest to the candidate points
        nearest, distances = np.zeros(len(locations), dtype=int), np.zeros(len(locations))
        for i, y in enumerate(locations):
            dist = np.linalg.norm(dofs_x - y, 2, axis=1)
            nearest_dof_id = np.argmin(dist)

            nearest[i] = nearest_dof_id
            distances[i] = dist[nearest_dof_id]

        # We plan to ask each local vector for those values
        comm = mesh.mpi_comm()
        point_ranks = np.zeros_like(nearest)
        for i, d in enumerate(distances):
            point_ranks[i] = np.argmin(comm.allgather(d))
            
        # But when putting things together we should only use values from
        # those dofs that were closest
        mask = point_ranks == comm.rank 
        
        values = np.zeros(len(nearest), dtype=float)
        # So we ask: values <-- nearest from local, and filter: zero-ing
        # those that were far
        evals = lambda u, dof=nearest, values=values: (np.copyto(values, u.array_r[nearest]),
                                                       np.copyto(values, values*mask))

        self.probes = evals

        self.comm = pyMPI.COMM_WORLD
        self.values = values
        self.readings = np.zeros_like(self.values)
        # Return the value in the shape of vector/matrix
        self.nprobes = len(locations)

    def sample(self, u):
        '''Evaluate the probes listing the time as t'''
        u_vec = as_backend_type(u.vector())
        u_vec = u_vec.vec()  # This is PETSc
        
        self.probes(u_vec)
        # Summing Zero and the Right value
        self.comm.Reduce(self.values, self.readings, op=pyMPI.SUM)  # Sync

        return self.readings.reshape((self.nprobes, -1))

# --------------------------------------------------------------------

if __name__ == '__main__':
    from dolfin import *

    
    parameters['ghost_mode'] = 'shared_facet'

    
    mesh_file = 'tile_1_hein_GMSH307_10_1.h5'
    # Get mesh to setup the ode solver
    comm = MPI.comm_world  # FIXME!
    h5 = HDF5File(comm, mesh_file, 'r')
    mesh = Mesh()
    h5.read(mesh, 'mesh', False)
    
    # Sample a tidge
    pts = np.array([probe_cell_at(m=m, n=0, level=2, point=p)
                    for m in range(2) for p in (-1, 1)])

    V = FunctionSpace(mesh, 'CG', 1)
    d = Function(V)

    f = Expression('t*(x[0] + 2*x[1])', degree=1, t=1)

    v = interpolate(f, V)

    probes = ApproxPointProbe(v, pts)

    stop = 1
    
    table = pts
    gdim = table.shape[1]
    for i in range(stop):
        f.t += i*0.1
        v.assign(interpolate(f, V))
        
        probe_values = probes.sample(v)
        if probes.comm.rank == 0:
            print(probe_values.flatten())
            table = np.c_[table, probe_values.flatten()]

    if probes.comm.rank == 0:
        x, y, z = pts.T

        t = 1
        error = 0
        for i in range(stop):
            t += i*0.1
            v = t*(x + 2*y)
            error = max(error, np.linalg.norm(v - table[:, gdim+i]))
        print('>>>', error)
