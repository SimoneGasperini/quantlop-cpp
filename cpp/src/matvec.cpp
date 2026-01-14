#include "matvec.hpp"

MatVec::MatVec(Complex c, String str)
    : coeff(c),
      string(str),
      flip_mask(0),
      y_mask(0),
      z_mask(0)
{
    const Size nq = string.size();
    for (Index q = 0; q < nq; ++q)
    {
        const Mask bit = Mask(1) << (nq - 1 - q);
        switch (string[q])
        {
        case 'X':
            flip_mask |= bit;
            break;
        case 'Y':
            flip_mask |= bit;
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

void MatVec::operator()(const Complex *in, Complex *out) const
{
    const Complex I(0.0, 1.0);
    const Size dim = Size(1) << string.size();
    const int y_count = std::popcount(y_mask);
    for (Index i = 0; i < dim; ++i)
    {
        const Mask j = Mask(i) ^ flip_mask;
        const Mask yz_mask = y_mask | z_mask;
        const int parity = std::popcount(yz_mask & i) & 1;
        Complex phase;
        switch (y_count & 3)
        {
        case 0:
            phase = 1;
            break;
        case 1:
            phase = I;
            break;
        case 2:
            phase = -1;
            break;
        default:
            phase = -I;
            break;
        }
        phase = parity ? -phase : phase;
        out[j] += coeff * phase * in[i];
    }
}
