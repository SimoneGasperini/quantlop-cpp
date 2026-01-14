#pragma once

#include <cstddef>
#include "utils.hpp"

class PauliWord
{
public:
    PauliWord(Complex coeff, std::string string);

    Complex coeff() const;
    const std::string &string() const;
    std::size_t num_qubits() const;
    Vector matvec(const Vector &vec) const;

    PauliWord operator*(Complex c) const;
    friend PauliWord operator*(Complex c, const PauliWord &pw);

private:
    friend class Hamiltonian;
    static MatVecFn make_matvec(Complex coeff, std::string string);

    Complex _coeff;
    std::string _string;
    MatVecFn _matvec;
};
