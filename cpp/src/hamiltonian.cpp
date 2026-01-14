#include "hamiltonian.hpp"

Hamiltonian::Hamiltonian(std::vector<PauliWord> pws)
    : pwords(std::move(pws)) {}

Size Hamiltonian::num_qubits() const { return pwords.front().num_qubits(); }

void Hamiltonian::matvec_into(const Complex *in, Complex *out) const
{
    Size dim = Size(1) << num_qubits();
    std::fill(out, out + dim, 0.0);
    for (const PauliWord &pw : pwords)
    {
        pw.matvec(in, out);
    }
}

Hamiltonian Hamiltonian::operator*(Complex c) const
{
    std::vector<PauliWord> pws;
    pws.reserve(pwords.size());
    for (const PauliWord &pw : pwords)
    {
        pws.push_back(pw * c);
    }
    return Hamiltonian(std::move(pws));
}

Hamiltonian operator*(Complex c, const Hamiltonian &ham) { return ham * c; }

double Hamiltonian::lcu_norm() const
{
    double norm = 0.0;
    for (const PauliWord &pw : pwords)
    {
        norm += std::abs(pw.coeff);
    }
    return norm;
}
