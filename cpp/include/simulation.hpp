#pragma once
#include "include.hpp"

Complex *expm_multiply_higham(const Hamiltonian &ham, const Complex *psi);
Complex *expm_multiply_krylov(const Hamiltonian &ham, const Complex *psi, Complex coeff);
Complex *evolve(const Hamiltonian &ham, const Complex *psi, Complex coeff, const std::string &method);
