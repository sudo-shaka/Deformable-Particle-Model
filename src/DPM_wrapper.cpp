#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_Cell(py::module &);
void init_monolayer(py::module &);

namespace dpm {

    PYBIND11_MODULE(DPM, m){
        m.doc() = "Deformable Particle Model";
        init_Cell(m);
        init_monolayer(m);
    }
}
