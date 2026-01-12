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

Vector Hamiltonian::matvec(const Vector &vec) const { return _matvec(vec); }

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

    return [matvecs = std::move(matvecs)](const Vector &vec)
    {
        const int dim = static_cast<int>(vec.size());
        Vector out(dim, Complex());
        for (const auto &matvec : matvecs)
        {
            const Vector partial = matvec(vec);
            for (int i = 0; i < dim; ++i)
            {
                out[i] += partial[i];
            }
        }
        return out;
    };
}
