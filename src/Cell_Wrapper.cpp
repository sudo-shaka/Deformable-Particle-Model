#include "../include/Cell/Cell.hpp"
#include <pybind11/stl.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_Cell(py::module &m){
    py::class_<DPM::Cell>(m, "Cell")
    .def(py::init<double, double, double, int, double, double,double,double,double,double,double,double>(),
    py::arg("x1"),
    py::arg("y1"),
    py::arg("calA0"),
    py::arg("NV"),
    py::arg("Kl"),
    py::arg("Kb"),
    py::arg("Ka"),
    py::arg("v0"),
    py::arg("Dr"),
    py::arg("Ds"),
    py::arg("a0"),
    py::arg("psi")
    )
    .def(py::init<int>(),py::arg("NV"))
    .def_readwrite("v0",&DPM::Cell::v0)
    .def_readwrite("Ds",&DPM::Cell::Dr)
    .def_readwrite("Dr",&DPM::Cell::Ds)
    .def_readwrite("Kl",&DPM::Cell::Kl)
    .def_readwrite("Kb",&DPM::Cell::Kb)
    .def_readwrite("Ka",&DPM::Cell::Ka)
    .def_readwrite("a0",&DPM::Cell::a0)
    .def_readwrite("l0",&DPM::Cell::l0)
    .def_readwrite("r0",&DPM::Cell::r0)
    .def_readwrite("X",&DPM::Cell::X)
    .def_readwrite("Y",&DPM::Cell::Y)
    .def_readwrite("Fx",&DPM::Cell::Fx)
    .def_readwrite("Fy",&DPM::Cell::Fy)
    .def_readwrite("l1",&DPM::Cell::l1)
    .def_readwrite("l2",&DPM::Cell::l2)
    .def_readwrite("psi",&DPM::Cell::psi)
    .def_readonly("radii",&DPM::Cell::radii)
    .def("ResetForces",&DPM::Cell::ResetForces)
    .def("PerimeterForceUpdate",&DPM::Cell::PerimeterForceUpdate)
    .def("AreaForceUpdate", &DPM::Cell::AreaForceUpdate)
    .def("BendingForceUpdate",&DPM::Cell::BendingForceUpdate)
    .def("DrivingForceUpdate",py::overload_cast<double>(&DPM::Cell::DrivingForceUpdate),py::arg("dt"))
    .def("UpdateShapeForces",&DPM::Cell::UpdateShapeForces)
    .def("UpdateEuler",py::overload_cast<int,double>(&DPM::Cell::UpdateEuler),py::arg("nsteps"),py::arg("dt"))
    .def("UpdateVV",py::overload_cast<int,double>(&DPM::Cell::UpdateVV),py::arg("nsteps"),py::arg("dt"))
    .def("UpdateDirectorDiffusion",py::overload_cast<double>(&DPM::Cell::UpdateDirectorDiffusion),py::arg("dt"))
    .def("SetCalA0",py::overload_cast<double>(&DPM::Cell::SetCalA0),py::arg("calA0"))
    .def("GetPerim", &DPM::Cell::GetPerim)
    .def("GetArea",&DPM::Cell::GetArea)
    .def("GetCenterX",&DPM::Cell::GetCenterX)
    .def("GetCenterY",&DPM::Cell::GetCenterY)
    .def("GetShapeParameter", &DPM::Cell::GetShapeParam)
    .def("isConvex",&DPM::Cell::isConvex)
;
}
