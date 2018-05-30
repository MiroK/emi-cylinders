from dolfin import compile_extension_module

code="""
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/MeshValueCollection.h>

namespace dolfin {

  void fill_mvc_from_mf(const std::shared_ptr<MeshFunction<std::size_t>> mesh_f,
                        std::shared_ptr<MeshValueCollection<std::size_t>> mvc)
  {
    std::size_t value;
    std::size_t* values = mesh_f->values();
    for(std::size_t i = 0; i < mesh_f->size(); i++){
        value = *values;
        if(value != 0){
            mvc->set_value(i, value);
        }
        values++;
    }
  }
};
"""

module = compile_extension_module(code)

# Exports
fill_mvc_from_mf = module.fill_mvc_from_mf

