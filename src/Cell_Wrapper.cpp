#include "../include/Cell/Cell.hpp"
#include <pybind11/stl.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_Cell(py::module &m){
    py::class_<DPM::Cell>(m, "Cell")
    .def(py::init<int, double, double, double, int, double, double,double,double,double,double,double,double>(),
    py::arg("idx"),
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
    .def_readwrite("v0",&DPM::Cell::v0)
    .def_readwrite("Kl",&DPM::Cell::Kl)
    .def_readwrite("Kb",&DPM::Cell::Kb)
    .def_readwrite("Ka",&DPM::Cell::Ka)
    .def_readwrite("a0",&DPM::Cell::a0)
    .def_readwrite("r0",&DPM::Cell::r0)
    .def_readwrite("X",&DPM::Cell::X)
    .def_readwrite("Y",&DPM::Cell::Y)
    .def("PerimeterForceUpdate",&DPM::Cell::PerimeterForceUpdate)
    .def("AreaForceUpdate", &DPM::Cell::AreaForceUpdate)
    .def("BendingForceUpdate",&DPM::Cell::BendingForceUpdate)
    .def("DrivingForceUpdate",py::overload_cast<double>(&DPM::Cell::DrivingForceUpdate),py::arg("dt"))
    .def("UpdateShapeForces",&DPM::Cell::UpdateShapeForces)
    .def("UpdateEuler",py::overload_cast<double>(&DPM::Cell::UpdateEuler),py::arg("dt"))
    .def("UpdateDirectorDiffusion",py::overload_cast<double>(&DPM::Cell::UpdateDirectorDiffusion),py::arg("dt"))
    .def("UpdateDirectorPsi",py::overload_cast<double>(&DPM::Cell::UpdateDirectorPsi),py::arg("psi"))
    .def("UpdateStickyness",py::overload_cast<int, double, double>(&DPM::Cell::UpdateStickness),py::arg("vertexindex"),py::arg("l1"),py::arg("l2"))
    .def("GetPerim", &DPM::Cell::GetPerim)
    .def("GetArea",&DPM::Cell::GetArea)
    .def("GetCenterX",&DPM::Cell::GetCenterX)
    .def("GetCenterY",&DPM::Cell::GetCenterY)
    .def("GetShapeParameter", &DPM::Cell::GetShapeParam)
    .def("isConvex",&DPM::Cell::isConvex)
;
}
