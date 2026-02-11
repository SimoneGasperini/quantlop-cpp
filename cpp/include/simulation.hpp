#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include "include.hpp"
#include "hamiltonian.hpp"

Complex *expm_multiply_higham(const Hamiltonian &ham, const Complex *psi);
Complex *expm_multiply_krylov(const Hamiltonian &ham, const Complex *psi);

inline Complex *evolve(const Hamiltonian &ham, const Complex *psi, Complex coeff, const std::string &method)
{
    const Complex i(0.0, 1.0);
    const Hamiltonian expm = (-i * coeff) * ham;
    if (method == "higham")
    {
        return expm_multiply_higham(expm, psi);
    }
    if (method == "krylov")
    {
        return expm_multiply_krylov(expm, psi);
    }
    throw std::invalid_argument("Invalid method: " + method);
}
