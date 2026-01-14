#pragma once

#include <cstddef>
#include <string>

#include "utils.hpp"

class PauliWord
{
public:
    PauliWord(Complex coeff, std::string string);

    Complex coeff() const;
    const std::string &string() const;
    std::size_t num_qubits() const;
    void matvec_into(const Complex *in, Complex *out) const;

    PauliWord operator*(Complex c) const;
    friend PauliWord operator*(Complex c, const PauliWord &pw);

private:
    friend class Hamiltonian;
    static MatVecFn make_matvec(Complex coeff, std::string string);

    Complex _coeff;
    std::string _string;
    MatVecFn _matvec;
};
