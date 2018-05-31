# Keeping fenics 2017 and 2018 compatibility
try:
    from dolfin import compile_extension_module as compile_cpp
    is_2017 = True
except ImportError:
    from dolfin import compile_cpp_code as compile_cpp
    is_2017 = False

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

if not is_2017:
    includes = '''
    #include <pybind11/pybind11.h>
    '''
    
    pybind = '''
    PYBIND11_MODULE(SIGNATURE, m)
    {
        m.def("fill_mvc_from_mf", &dolfin::fill_mvc_from_mf);
    }
    '''

    code = '\n'.join([includes, code, pybind])

module = compile_cpp(code)

# Exports
fill_mvc_from_mf = module.fill_mvc_from_mf

