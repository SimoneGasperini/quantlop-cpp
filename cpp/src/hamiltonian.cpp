#include "hamiltonian.hpp"

#include <utility>

Hamiltonian::Hamiltonian(std::vector<PauliWord> pwords)
    : _pwords(std::move(pwords)),
      _dim(std::size_t{1} << _pwords.front().string().size()) {}

const std::vector<PauliWord> &Hamiltonian::pwords() const { return _pwords; }

std::vector<std::complex<double>> Hamiltonian::coeffs() const
{
    std::vector<std::complex<double>> out;
    out.reserve(_pwords.size());
    for (const PauliWord &pw : _pwords)
    {
        out.push_back(pw.coeff());
    }
    return out;
}

std::vector<std::string> Hamiltonian::strings() const
{
    std::vector<std::string> out;
    out.reserve(_pwords.size());
    for (const PauliWord &pw : _pwords)
    {
        out.push_back(pw.string());
    }
    return out;
}

void Hamiltonian::matvec_into(const std::complex<double> *in, std::complex<double> *out) const
{
    for (std::size_t i = 0; i < _dim; ++i)
    {
        out[i] = std::complex<double>();
    }
    for (const auto &pw : _pwords)
    {
        pw._matvec(in, out);
    }
}

Hamiltonian Hamiltonian::operator*(std::complex<double> c) const
{
    std::vector<PauliWord> pwords;
    pwords.reserve(_pwords.size());
    for (const PauliWord &pw : _pwords)
    {
        pwords.push_back(pw * c);
    }
    return Hamiltonian(std::move(pwords));
}

Hamiltonian operator*(std::complex<double> c, const Hamiltonian &ham) { return ham * c; }

double Hamiltonian::lcu_norm() const
{
    double norm = 0.0;
    for (const PauliWord &pw : _pwords)
    {
        norm += std::abs(pw.coeff());
    }
    return norm;
}
