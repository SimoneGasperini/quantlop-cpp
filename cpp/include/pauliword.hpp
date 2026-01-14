#pragma once

#include <cstddef>
#include <complex>
#include <string>

#include "matvec.hpp"

class PauliWord
{
public:
    PauliWord(std::complex<double> coeff, std::string string);

    std::complex<double> coeff() const;
    const std::string &string() const;
    std::size_t num_qubits() const;

    PauliWord operator*(std::complex<double> c) const;
    friend PauliWord operator*(std::complex<double> c, const PauliWord &pw);

private:
    friend class Hamiltonian;
    std::complex<double> _coeff;
    std::string _string;
    MatVec _matvec;
};
