try:
    from dolfin import compile_extension_module as compile_cpp

    code="""
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/CellType.h>
#include <dolfin/mesh/MeshTopology.h>
#include <dolfin/mesh/MeshConnectivity.h>
#include <dolfin/mesh/MeshValueCollection.h>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <map>

namespace dolfin {
  // Fills a SIMPLICIAL mesh
  void fill_mesh(const Array<double>& coordinates,
                 const Array<std::size_t>& cells, 
                 const int tdim, 
                 const int gdim, 
                 std::shared_ptr<Mesh> mesh)
  {
     int nvertices = coordinates.size()/gdim;     

     int nvertices_per_cell = tdim + 1;
     int ncells = cells.size()/nvertices_per_cell;   

     MeshEditor editor;
     if (tdim == 1){
         editor.open(*mesh, CellType::Type::interval, tdim, gdim);
     }
     else if (tdim == 2){
         editor.open(*mesh, CellType::Type::triangle, tdim, gdim);
     }
     else{
         editor.open(*mesh, CellType::Type::tetrahedron, tdim, gdim);
     }

     editor.init_vertices(nvertices);
     editor.init_cells(ncells);

     std::vector<double> vertex(gdim);
     for(std::size_t index = 0; index < nvertices; index++){
         for(std::size_t i = 0; i < gdim; i++){
             vertex[i] = coordinates[gdim*index  + i];
         }
         editor.add_vertex(index, vertex);
     }

     std::vector<std::size_t> cell(nvertices_per_cell);
     for(std::size_t index = 0; index < ncells; index++){
         for(std::size_t i = 0; i < nvertices_per_cell; i++){
             cell[i] = cells[nvertices_per_cell*index  + i];
         }
         editor.add_cell(index, cell);
     }

     editor.close();
  }

  // Assign tag to entities of mesh function over simplicial mesh
  void fill_mesh_function(const std::shared_ptr<Mesh> mesh, 
                          const Array<std::size_t>& indices,
                          const int tdim,
                          const std::size_t tag,
                          std::shared_ptr<MeshFunction<std::size_t>> values)
  {
    const MeshConnectivity& v2entity = (mesh->topology())(0, tdim);

    const int nindices_per_entity = tdim + 1;
    const int n_entities = indices.size()/nindices_per_entity;

    std::vector<unsigned int> entities0, entities1, intersect;
    std::vector<unsigned int>::iterator last;
    std::size_t v, count0, count1, entity;

    for(std::size_t e=0; e < n_entities; e++){
        // There's always at least one entity
        v = indices[e*nindices_per_entity];
            
        count0 = v2entity.size(v);
        entities0.assign(v2entity(v), v2entity(v) + count0);
        // The tagged entity is the one which is connected to all 
        // the vertices
        for(std::size_t ei=1; ei < nindices_per_entity; ei++){
            v = indices[e*nindices_per_entity + ei];
            
            count1 = v2entity.size(v);
            entities1.assign(v2entity(v), v2entity(v) + count1);
            // set_intersection needs sorted containers
            std::sort(entities0.begin(), entities0.end());
            std::sort(entities1.begin(), entities1.end());
            // Alloc
            intersect.resize(std::max(count0, count1));

            // intersection of the two
            last = std::set_intersection(entities0.begin(), entities0.end(),
                                         entities1.begin(), entities1.end(),
                                         intersect.begin());
            // Next time around we compare with the intersection
            entities0.assign(intersect.begin(), last);
            count0 = last - intersect.begin();
        }
        entity = *entities0.begin();
        // Tag it
        (*values)[entity] = tag;
    }
  }
};
"""

except ImportError:
    from dolfin import compile_cpp_code as compile_cpp

    code="""
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshEditor.h>
#include <dolfin/mesh/CellType.h>
#include <dolfin/mesh/MeshTopology.h>
#include <dolfin/mesh/MeshConnectivity.h>
#include <dolfin/mesh/MeshValueCollection.h>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <map>

#include <Eigen/Core>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

using IntVecIn = Eigen::Ref<const Eigen::VectorXi>;
using DoubleVecIn = Eigen::Ref<const Eigen::VectorXd>;

namespace dolfin {
  // Fills a SIMPLICIAL mesh
  void fill_mesh(const DoubleVecIn coordinates,
                 const IntVecIn cells, 
                 const int tdim, 
                 const int gdim, 
                 std::shared_ptr<Mesh> mesh)
  {
     int nvertices = coordinates.size()/gdim;     

     int nvertices_per_cell = tdim + 1;
     int ncells = cells.size()/nvertices_per_cell;   

     MeshEditor editor;
     if (tdim == 1){
         editor.open(*mesh, CellType::Type::interval, tdim, gdim);
     }
     else if (tdim == 2){
         editor.open(*mesh, CellType::Type::triangle, tdim, gdim);
     }
     else{
         editor.open(*mesh, CellType::Type::tetrahedron, tdim, gdim);
     }

     editor.init_vertices(nvertices);
     editor.init_cells(ncells);

     std::vector<double> vertex(gdim);
     for(std::size_t index = 0; index < nvertices; index++){
         for(std::size_t i = 0; i < gdim; i++){
             vertex[i] = coordinates[gdim*index  + i];
         }
         editor.add_vertex(index, vertex);
     }

     std::vector<std::size_t> cell(nvertices_per_cell);
     for(std::size_t index = 0; index < ncells; index++){
         for(std::size_t i = 0; i < nvertices_per_cell; i++){
             cell[i] = cells[nvertices_per_cell*index  + i];
         }
         editor.add_cell(index, cell);
     }

     editor.close();
  }

  // Assign tag to entities of mesh function over simplicial mesh
  void fill_mesh_function(const std::shared_ptr<Mesh> mesh, 
                          const IntVecIn indices,
                          const int tdim,
                          const std::size_t tag,
                          std::shared_ptr<MeshFunction<std::size_t>> values)
  {
    const MeshConnectivity& v2entity = (mesh->topology())(0, tdim);

    const int nindices_per_entity = tdim + 1;
    const int n_entities = indices.size()/nindices_per_entity;

    std::vector<unsigned int> entities0, entities1, intersect;
    std::vector<unsigned int>::iterator last;
    std::size_t v, count0, count1, entity;

    for(std::size_t e=0; e < n_entities; e++){
        // There's always at least one entity
        v = indices[e*nindices_per_entity];
            
        count0 = v2entity.size(v);
        entities0.assign(v2entity(v), v2entity(v) + count0);
        // The tagged entity is the one which is connected to all 
        // the vertices
        for(std::size_t ei=1; ei < nindices_per_entity; ei++){
            v = indices[e*nindices_per_entity + ei];
            
            count1 = v2entity.size(v);
            entities1.assign(v2entity(v), v2entity(v) + count1);
            // set_intersection needs sorted containers
            std::sort(entities0.begin(), entities0.end());
            std::sort(entities1.begin(), entities1.end());
            // Alloc
            intersect.resize(std::max(count0, count1));

            // intersection of the two
            last = std::set_intersection(entities0.begin(), entities0.end(),
                                         entities1.begin(), entities1.end(),
                                         intersect.begin());
            // Next time around we compare with the intersection
            entities0.assign(intersect.begin(), last);
            count0 = last - intersect.begin();
        }
        entity = *entities0.begin();
        // Tag it
        (*values)[entity] = tag;
    }
  }
};

PYBIND11_MODULE(SIGNATURE, m)
{
    m.def("fill_mesh", &dolfin::fill_mesh);
    m.def("fill_mesh_function", &dolfin::fill_mesh_function);
}

"""


# -----------

module = compile_cpp(code)

# Exports
fill_mesh = module.fill_mesh
fill_mesh_function = module.fill_mesh_function
