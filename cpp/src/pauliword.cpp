#include "pauliword.hpp"

PauliWord::PauliWord(Complex coeff, std::string string)
    : _coeff(coeff),
      _string(std::move(string)),
      _matvec(make_matvec(_coeff, _string)) {}

Complex PauliWord::coeff() const { return _coeff; }

const std::string &PauliWord::string() const { return _string; }

int PauliWord::num_qubits() const { return _string.size(); }

Vector PauliWord::matvec(const Vector &vec) const { return _matvec(vec); }

PauliWord PauliWord::operator*(Complex c) const { return PauliWord(_coeff * c, _string); }

PauliWord operator*(Complex c, const PauliWord &pw) { return pw * c; }

MatVecFn PauliWord::make_matvec(Complex coeff, std::string string)
{
    return [coeff, string](const Vector &vec)
    {
        const int nq = string.size();
        const int dim = 1 << nq;
        Vector out(dim, Complex());
        for (int i = 0; i < dim; ++i)
        {
            int j = i;
            Complex phase(1.0, 0.0);
            for (int q = 0; q < nq; ++q)
            {
                const int bit_pos = nq - 1 - q;
                const bool bit = ((i >> bit_pos) & 1) != 0;
                switch (string[q])
                {
                case 'I':
                    break;
                case 'X':
                    j ^= (1 << bit_pos);
                    break;
                case 'Y':
                    j ^= (1 << bit_pos);
                    phase *= bit ? Complex(0.0, -1.0) : Complex(0.0, 1.0);
                    break;
                case 'Z':
                    phase *= bit ? -1.0 : 1.0;
                    break;
                }
            }
            out[j] += coeff * phase * vec[i];
        }
        return out;
    };
}
