#include "matvec.hpp"

#include <cstddef>
#include <utility>

MatVec::MatVec(std::complex<double> coeff, std::string string)
    : _coeff(coeff),
      _string(std::move(string)) {}

void MatVec::operator()(const std::complex<double> *in, std::complex<double> *out) const
{
    const std::size_t nq = _string.size();
    const std::size_t dim = std::size_t{1} << nq;
    for (std::size_t i = 0; i < dim; ++i)
    {
        std::size_t j = i;
        std::complex<double> phase(1.0, 0.0);
        for (std::size_t q = 0; q < nq; ++q)
        {
            const std::size_t bit_pos = nq - 1 - q;
            const bool bit = ((i >> bit_pos) & std::size_t{1}) != 0;
            switch (_string[q])
            {
            case 'I':
                break;
            case 'X':
                j ^= (std::size_t{1} << bit_pos);
                break;
            case 'Y':
                j ^= (std::size_t{1} << bit_pos);
                phase *= bit ? std::complex<double>(0.0, -1.0) : std::complex<double>(0.0, 1.0);
                break;
            case 'Z':
                phase *= bit ? -1.0 : 1.0;
                break;
            }
        }
        out[j] += _coeff * phase * in[i];
    }
}
