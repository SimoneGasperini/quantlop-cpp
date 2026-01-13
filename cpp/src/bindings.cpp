#include <algorithm>
#include <complex>
#include <cstddef>
#include <string>
#include <vector>

#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "hamiltonian.hpp"
#include "pauliword.hpp"
#include "simulation.hpp"
#include "statevector.hpp"
#include "utils.hpp"

namespace py = pybind11;

static py::array_t<Complex>
evolve_py(const Hamiltonian &ham,
          py::array_t<Complex, py::array::c_style | py::array::forcecast> psi,
          Complex coeff)
{
    py::buffer_info info = psi.request();

    const int dim = static_cast<int>(info.shape[0]);
    const int nq = static_cast<int>(std::log2(dim));
    const auto *ptr = static_cast<const Complex *>(info.ptr);
    std::vector<Complex> data(ptr, ptr + dim);
    StateVector state(nq, std::move(data));
    StateVector out = evolve(ham, state, coeff);

    py::array_t<Complex> result(dim);
    py::buffer_info out_info = result.request();
    auto *out_ptr = static_cast<Complex *>(out_info.ptr);
    const auto &out_data = out.data();
    std::copy(out_data.begin(), out_data.end(), out_ptr);
    return result;
}

static py::array_t<Complex>
trace_evolve_py(const Hamiltonian &ham,
          py::array_t<Complex, py::array::c_style | py::array::forcecast> psi,
          Complex coeff)
{
    py::buffer_info info = psi.request();

    const int dim = static_cast<int>(info.shape[0]);
    const int nq = static_cast<int>(std::log2(dim));
    const auto *ptr = static_cast<const Complex *>(info.ptr);
    std::vector<Complex> data(ptr, ptr + dim);
    StateVector state(nq, std::move(data));
    std::vector<StateVector> out = trace_evolve(ham, state, coeff);

    py::array_t<Complex> result({static_cast<py::ssize_t>(out.size()),
                             static_cast<py::ssize_t>(dim)});
    auto r = result.mutable_unchecked<2>();

    //py::buffer_info out_info = result.request();
    //auto **out_ptr = static_cast<Complex **>(out_info.ptr);
    //std::vector<StateVector>::iterator it;
    //for (it = out.begin(); it != out.end(); ++it) {
    //    const auto &it_data = it->data();
    //    std::copy(it_data.begin(), it_data.end(), out_ptr[std::distance(out.begin(), it)]);
    //}

    for (py::ssize_t i = 0; i < (py::ssize_t)out.size(); ++i) {
       const auto& v = out[i].data();              // assuming vector<Complex>
       std::copy(v.begin(), v.end(), &r(i, 0));
    }
    //const auto &out_data = out.data();
    //std::copy(out_data.begin(), out_data.end(), out_ptr);
    return result;
}


PYBIND11_MODULE(_core, m)
{
    m.doc() = "Quantlop C++ core bindings";

    py::class_<PauliWord>(m, "PauliWord")
        .def(py::init<Complex, std::string>(), py::arg("coeff"), py::arg("string"));

    py::class_<Hamiltonian>(m, "Hamiltonian")
        .def(py::init<std::vector<PauliWord>>(), py::arg("pauli_words"));

    m.def("evolve", &evolve_py, py::arg("ham"), py::arg("psi"),
          py::arg("coeff") = Complex(1.0, 0.0));
    m.def("trace_evolve", &trace_evolve_py, py::arg("ham"), py::arg("psi"),
          py::arg("coeff") = Complex(1.0, 0.0));

}
