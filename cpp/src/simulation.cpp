#include "hamiltonian.hpp"
#include "simulation.hpp"

Complex *evolve(const Hamiltonian &ham, const Complex *psi, Complex coeff, const std::string &method)
{
    const Complex i(0.0, 1.0);
    const Hamiltonian expm = (-i * coeff) * ham;

    if (method == "higham") return expm_multiply_higham(expm, psi);
    if (method == "krylov") return expm_multiply_krylov(expm, psi);
}
