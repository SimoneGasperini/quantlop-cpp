#include "pauliword.hpp"

PauliWord::PauliWord(Complex coeff, std::string string)
    : _coeff(coeff),
      _string(std::move(string)),
      _matvec(make_matvec(_coeff, _string)) {}

Complex PauliWord::coeff() const { return _coeff; }

const std::string &PauliWord::string() const { return _string; }

std::size_t PauliWord::num_qubits() const { return _string.size(); }

void PauliWord::matvec_into(const Complex *in, Complex *out) const { _matvec(in, out); }

PauliWord PauliWord::operator*(Complex c) const { return PauliWord(_coeff * c, _string); }

PauliWord operator*(Complex c, const PauliWord &pw) { return pw * c; }

MatVecFn PauliWord::make_matvec(Complex coeff, std::string string)
{
    return [coeff, string](const Complex *in, Complex *out)
    {
        const std::size_t nq = string.size();
        const std::size_t dim = std::size_t{1} << nq;
        for (std::size_t i = 0; i < dim; ++i)
        {
            std::size_t j = i;
            Complex phase(1.0, 0.0);
            for (std::size_t q = 0; q < nq; ++q)
            {
                const std::size_t bit_pos = nq - 1 - q;
                const bool bit = ((i >> bit_pos) & std::size_t{1}) != 0;
                switch (string[q])
                {
                case 'I':
                    break;
                case 'X':
                    j ^= (std::size_t{1} << bit_pos);
                    break;
                case 'Y':
                    j ^= (std::size_t{1} << bit_pos);
                    phase *= bit ? Complex(0.0, -1.0) : Complex(0.0, 1.0);
                    break;
                case 'Z':
                    phase *= bit ? -1.0 : 1.0;
                    break;
                }
            }
            out[j] += coeff * phase * in[i];
        }
    };
}
