#pragma once
#include "include.hpp"
#include "hamiltonian.hpp"

Complex *expm_multiply_krylov(const Hamiltonian &ham, const Complex *psi, Complex coeff);
Complex *evolve(const Hamiltonian &ham, const Complex *psi, Complex coeff);
