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
#include "utils.hpp"

namespace py = pybind11;

static py::array_t<Complex>
evolve_py(const Hamiltonian &ham,
          py::array_t<Complex, py::array::c_style | py::array::forcecast> psi,
          Complex coeff)
{
    py::buffer_info info = psi.request();

    const std::size_t dim = static_cast<std::size_t>(info.shape[0]);
    const auto *ptr = static_cast<const Complex *>(info.ptr);
    py::array_t<Complex> result(static_cast<py::ssize_t>(dim));
    py::buffer_info out_info = result.request();
    auto *out_ptr = static_cast<Complex *>(out_info.ptr);
    evolve_into(ham, ptr, dim, out_ptr, coeff);
    return result;
}

static py::array_t<Complex>
trace_evolve_py(const Hamiltonian &ham,
                py::array_t<Complex, py::array::c_style | py::array::forcecast> psi,
                Complex coeff)
{
    py::buffer_info info = psi.request();

    const std::size_t dim = static_cast<std::size_t>(info.shape[0]);
    const auto *ptr = static_cast<const Complex *>(info.ptr);
    const std::size_t max_steps = trace_evolve_max_steps(ham);
    std::vector<Complex> out(max_steps * dim);
    const std::size_t steps = trace_evolve_into(ham, ptr, dim, out.data(), coeff);

    py::array_t<Complex> result({static_cast<py::ssize_t>(steps),
                                 static_cast<py::ssize_t>(dim)});
    auto r = result.mutable_unchecked<2>();

    // py::buffer_info out_info = result.request();
    // auto **out_ptr = static_cast<Complex **>(out_info.ptr);
    // std::vector<StateVector>::iterator it;
    // for (it = out.begin(); it != out.end(); ++it) {
    //     const auto &it_data = it->data();
    //     std::copy(it_data.begin(), it_data.end(), out_ptr[std::distance(out.begin(), it)]);
    // }

    for (py::ssize_t i = 0; i < static_cast<py::ssize_t>(steps); ++i)
    {
        const Complex *row = out.data() + (static_cast<std::size_t>(i) * dim);
        std::copy(row, row + dim, &r(i, 0));
    }
    // const auto &out_data = out.data();
    // std::copy(out_data.begin(), out_data.end(), out_ptr);
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
