#include "../include/Cell/monolayer.hpp"
#include <pybind11/stl.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_monolayer(py::module &m){
    py::class_<DPM::monolayer>(m, "monolayer")
    .def(py::init<std::vector<DPM::Cell>, double>(),
    py::arg("cells"),
    py::arg("phi0")
    )
    .def_readonly("phi0",&DPM::monolayer::phi0)
    .def_readonly("NCELLS",&DPM::monolayer::NCELLS)
    .def_readonly("vertDOF",&DPM::monolayer::VertDOF)
    .def_readonly("PotentialEnergy",&DPM::monolayer::U)
    .def_readwrite("Cells",&DPM::monolayer::Cells)
    .def_readwrite("BoxLength",&DPM::monolayer::L)
    .def_readwrite("Kc",&DPM::monolayer::Kc)
    .def_readwrite("Ftol",&DPM::monolayer::Ftol)
    .def_readwrite("dt0",&DPM::monolayer::dt0)
    .def_readonly("overlaps",&DPM::monolayer::overlaps)
    .def("disperse",&DPM::monolayer::disperse)
    .def("VertexFIRE",&DPM::monolayer::VertexFIRE)
    .def("UpdateEuler",py::overload_cast<int, double>(&DPM::monolayer::UpdateEuler),py::arg("nsteps"),py::arg("dt"))
    .def("UpdateVV",py::overload_cast<int, double>(&DPM::monolayer::UpdateVV),py::arg("nsteps"),py::arg("dt"))
    .def("InteractingForceUpdate",&DPM::monolayer::InteractingForceUpdate)
    .def("FindOverlaps",&DPM::monolayer::FindOverlaps)
    .def("GetPackingFraction",&::DPM::monolayer::GetPackingFraction)
    .def("ResetForces",&DPM::monolayer::ResetForces)
    .def("CellDivide",py::overload_cast<int>(&DPM::monolayer::CellDivision),py::arg("CellIdx"))
    ;
}
