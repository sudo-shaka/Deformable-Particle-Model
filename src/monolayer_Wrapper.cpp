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
    .def("disperse",&DPM::monolayer::disperse)
    .def("VertexFIRE",py::overload_cast<double, double, int, double>(&DPM::monolayer::VertexFIRE),py::arg("Alpha"),py::arg("dt"),py::arg("Max Iters"),py::arg("Force tolerance"))
    .def("UpdateEuler",py::overload_cast<int, double>(&DPM::monolayer::UpdateEuler),py::arg("nsteps"),py::arg("dt"))
    .def("UpdateEuler",py::overload_cast<double>(&DPM::monolayer::UpdateEuler),py::arg("dt"))
    .def("UpdateVV",py::overload_cast<int, double>(&DPM::monolayer::UpdateVV),py::arg("nsteps"),py::arg("dt"))
    .def("GetPackingFraction",&::DPM::monolayer::GetPackingFraction)
    .def("ResetForces",&DPM::monolayer::ResetForces)
    .def("CellDivide",py::overload_cast<int>(&DPM::monolayer::CellDivision),py::arg("CellIdx"))
    .def("RepulsiveForceUpdate",&DPM::monolayer::RepulsiveForces)
    .def("AttractiveForceUpdate",&DPM::monolayer::AttractiveForces)
    .def("RetractingForceUpdate",&DPM::monolayer::RetractingForceUpdate)
    .def("DrivingForceUpdate",&DPM::monolayer::DrivingForceUpdate)
    .def("InteractingForceUpdate",&DPM::monolayer::InteractingForceUpdate)
    .def("UpdateDirectorDiffusion",py::overload_cast<double>(&DPM::monolayer::UpdateDirectorDiffusion),py::arg("dt"))
    .def("UpdateDirectorVicsek",py::overload_cast<double, double>(&DPM::monolayer::UpdateDirectorVicsek),py::arg("eta"),py::arg("InteractingRadius"))
    .def("UpdateNearestNeighbors",&DPM::monolayer::NearestNeighborUpdate)
    .def("PinnedForceUpdate",&DPM::monolayer::PinnedForces)
    .def("UpdateShapeForces",&DPM::monolayer::ShapeForceUpdate)
    .def("SetInteractingMethod",&DPM::monolayer::SetInteractingMethod)
    ;
}
