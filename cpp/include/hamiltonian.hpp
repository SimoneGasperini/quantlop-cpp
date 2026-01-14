#pragma once

#include <cstddef>
#include <complex>
#include <string>
#include <vector>

#include "pauliword.hpp"

class Hamiltonian
{
public:
    Hamiltonian(std::vector<PauliWord> pwords);

    const std::vector<PauliWord> &pwords() const;
    std::vector<std::complex<double>> coeffs() const;
    std::vector<std::string> strings() const;
    void matvec_into(const std::complex<double> *in, std::complex<double> *out) const;

    Hamiltonian operator*(std::complex<double> c) const;
    friend Hamiltonian operator*(std::complex<double> c, const Hamiltonian &ham);

    double lcu_norm() const;

private:
    std::vector<PauliWord> _pwords;
    std::size_t _dim;
};
