#pragma once
#include "include.hpp"
#include "matvec.hpp"

class PauliWord
{
public:
    PauliWord(Complex c, String str);
    Size num_qubits() const;
    PauliWord operator*(Complex c) const;
    friend PauliWord operator*(Complex c, const PauliWord &pw);

private:
    friend class Hamiltonian;

    Complex coeff;
    String string;
    MatVec matvec;
};
