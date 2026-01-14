#include "hamiltonian.hpp"

Hamiltonian::Hamiltonian(std::vector<PauliWord> pwords)
    : _pwords(std::move(pwords)),
      _matvec(make_matvec(_pwords)) {}

const std::vector<PauliWord> &Hamiltonian::pwords() const { return _pwords; }

std::vector<Complex> Hamiltonian::coeffs() const
{
    std::vector<Complex> out;
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

void Hamiltonian::matvec_into(const Complex *in, Complex *out) const { _matvec(in, out); }

Hamiltonian Hamiltonian::operator*(Complex c) const
{
    std::vector<PauliWord> pwords;
    pwords.reserve(_pwords.size());
    for (const PauliWord &pw : _pwords)
    {
        pwords.push_back(pw * c);
    }
    return Hamiltonian(std::move(pwords));
}

Hamiltonian operator*(Complex c, const Hamiltonian &ham) { return ham * c; }

double Hamiltonian::lcu_norm() const
{
    double norm = 0.0;
    for (const PauliWord &pw : _pwords)
    {
        norm += std::abs(pw.coeff());
    }
    return norm;
}

MatVecFn Hamiltonian::make_matvec(const std::vector<PauliWord> &pwords)
{
    std::vector<MatVecFn> matvecs;
    matvecs.reserve(pwords.size());
    for (const auto &pw : pwords)
    {
        matvecs.push_back(PauliWord::make_matvec(pw.coeff(), pw.string()));
    }

    const std::size_t dim = pwords.empty()
                                ? 0
                                : (std::size_t{1} << pwords.front().string().size());
    return [matvecs = std::move(matvecs), dim](const Complex *in, Complex *out)
    {
        for (std::size_t i = 0; i < dim; ++i)
        {
            out[i] = Complex();
        }
        for (const auto &matvec : matvecs)
        {
            matvec(in, out);
        }
    };
}
