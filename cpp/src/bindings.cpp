#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "include.hpp"
#include "pauliword.hpp"
#include "hamiltonian.hpp"
#include "simulation.hpp"

namespace py = pybind11;

static py::array_t<Complex>
evolve_py(const Hamiltonian &ham,
          py::array_t<Complex, py::array::c_style | py::array::forcecast> psi,
          Complex coeff,
          const std::string &method)
{
    py::buffer_info info = psi.request();

    const Size dim = (info.shape[0]);
    const auto *ptr = static_cast<const Complex *>(info.ptr);
    Complex *out_ptr = evolve(ham, ptr, coeff, method);
    py::capsule owner(out_ptr, [](void *p)
                      { delete[] static_cast<Complex *>(p); });
    return py::array_t<Complex>(static_cast<py::ssize_t>(dim), out_ptr, owner);
}

PYBIND11_MODULE(quantlop_cpp, module_py)
{
    module_py.doc() = "Quantlop C++ core bindings";

    py::class_<PauliWord>(module_py, "PauliWord")
        .def(py::init<Complex, String>(), py::arg("coeff"), py::arg("string"));

    py::class_<Hamiltonian>(module_py, "Hamiltonian")
        .def(py::init<std::vector<PauliWord>>(), py::arg("pauli_words"));

    module_py.def("evolve", &evolve_py, py::arg("ham"), py::arg("psi"), py::arg("coeff"), py::arg("method"));
}
