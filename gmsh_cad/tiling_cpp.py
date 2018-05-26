from dolfin import compile_extension_module

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
         editor.open(*mesh, CellType::interval, tdim, gdim);
     }
     else if (tdim == 2){
         editor.open(*mesh, CellType::triangle, tdim, gdim);
     }
     else{
         editor.open(*mesh, CellType::tetrahedron, tdim, gdim);
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

  // Assign tag to entities of mesh value collection over simplicial mesh
  void fill_mesh_valuecollection(const std::shared_ptr<Mesh> mesh, 
                                 const Array<std::size_t>& indices,
                                 const int tdim,
                                 const std::size_t tag,
                                 std::shared_ptr<MeshValueCollection<std::size_t>> values)
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
        values->set_value(entity, tag);
    }
  }

  void fill_mf_from_mvc(const std::shared_ptr<MeshValueCollection<std::size_t>> mvc,
                        std::shared_ptr<MeshFunction<std::size_t>> mesh_f)
  {
    const std::shared_ptr<const Mesh> _mesh = mesh_f->mesh();
    const std::size_t _dim = mesh_f->dim();
    const std::size_t _size = _mesh->num_entities(_dim);
    _mesh->init(_dim);

    // Get mesh connectivity D --> d
    const std::size_t d = _dim;
    const std::size_t D = _mesh->topology().dim();
    dolfin_assert(d <= D);

    // Generate connectivity if it does not exist
    _mesh->init(D, d);
    const MeshConnectivity& connectivity = _mesh->topology()(D, d);
    dolfin_assert(!connectivity.empty());

    // Iterate over all values
    std::unordered_set<std::size_t> entities_values_set;
    std::map<std::pair<std::size_t, std::size_t>, std::size_t>::const_iterator it;
    const std::map<std::pair<std::size_t, std::size_t>, std::size_t>& values = mvc->values();
    for (it = values.begin(); it != values.end(); ++it)
    {
      // Get value collection entry data
      const std::size_t cell_index = it->first.first;
      const std::size_t local_entity = it->first.second;
      const std::size_t value = it->second;

      std::size_t entity_index = 0;
      if (d != D)
      {
        // Get global (local to to process) entity index
        dolfin_assert(cell_index < _mesh->num_cells());
        entity_index = connectivity(cell_index)[local_entity];
      }
      else
      {
        entity_index = cell_index;
        dolfin_assert(local_entity == 0);
      }

      dolfin_assert(entity_index < _size);
      mesh_f->set_value(entity_index, value);

      // Add entity index to set (used to check that all values are set)
      entities_values_set.insert(entity_index);
    }

    // Check that all values have been set, if not issue a debug message
    if (entities_values_set.size() != _size)
      dolfin_debug("Mesh value collection does not contain all values for all entities");

  }
};
"""

module = compile_extension_module(code)

# Exports
fill_mesh = module.fill_mesh
fill_mesh_function = module.fill_mesh_function
fill_mesh_valuecollection = module.fill_mesh_valuecollection
fill_mf_from_mvc = module.fill_mf_from_mvc
