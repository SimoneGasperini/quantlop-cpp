#include "matvec.hpp"

#include <cstddef>
#include <cstdint>
#include <utility>

MatVec::MatVec(std::complex<double> coeff, std::string string)
    : _coeff(coeff),
      _string(std::move(string)),
      _flip_mask(0),
      y_mask(0),
      z_mask(0)
{
    const std::size_t nq = _string.size();
    for (std::size_t q = 0; q < nq; ++q)
    {
        const std::size_t bit_pos = nq - 1 - q;
        const std::uint64_t bit = std::uint64_t{1} << bit_pos;
        switch (_string[q])
        {
        case 'X':
            _flip_mask |= bit;
            break;
        case 'Y':
            _flip_mask |= bit;
            y_mask |= bit;
            break;
        case 'Z':
            z_mask |= bit;
            break;
        default:
            break;
        }
    }
}

void MatVec::operator()(const std::complex<double> *in, std::complex<double> *out) const
{
    const std::size_t nq = _string.size();
    const std::size_t dim = std::size_t{1} << nq;
    const std::uint64_t flip_mask = _flip_mask;
    const std::uint64_t yz_mask = y_mask | z_mask;
    const int y_count = __builtin_popcountll(y_mask);
    const std::complex<double> base_phase = [&]()
    {
        switch (y_count & 3)
        {
        case 0:
            return std::complex<double>(1.0, 0.0);
        case 1:
            return std::complex<double>(0.0, 1.0);
        case 2:
            return std::complex<double>(-1.0, 0.0);
        default:
            return std::complex<double>(0.0, -1.0);
        }
    }();
    for (std::size_t i = 0; i < dim; ++i)
    {
        const std::uint64_t idx = static_cast<std::uint64_t>(i);
        const std::uint64_t j = idx ^ flip_mask;
        const int parity = __builtin_popcountll(yz_mask & idx) & 1;
        const std::complex<double> phase = parity ? -base_phase : base_phase;
        out[static_cast<std::size_t>(j)] += _coeff * phase * in[i];
    }
}
