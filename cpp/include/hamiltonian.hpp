#pragma once

#include "pauliword.hpp"
#include "utils.hpp"

class Hamiltonian
{
public:
    Hamiltonian(std::vector<PauliWord> pwords);

    const std::vector<PauliWord> &pwords() const;
    std::vector<Complex> coeffs() const;
    std::vector<std::string> strings() const;
    Vector matvec(const Vector &vec) const;

    Hamiltonian operator*(Complex c) const;
    friend Hamiltonian operator*(Complex c, const Hamiltonian &ham);

    double lcu_norm() const;

private:
    static MatVecFn make_matvec(const std::vector<PauliWord> &pwords);

    std::vector<PauliWord> _pwords;
    MatVecFn _matvec;
};
