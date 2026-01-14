#pragma once
#include "include.hpp"
#include "pauliword.hpp"

class Hamiltonian
{
public:
    Hamiltonian(std::vector<PauliWord> pws);
    Size num_qubits() const;
    void matvec_into(const Complex *in, Complex *out) const;
    Hamiltonian operator*(Complex c) const;
    friend Hamiltonian operator*(Complex c, const Hamiltonian &ham);
    double lcu_norm() const;

private:
    std::vector<PauliWord> pwords;
};
