#pragma once
#include "include.hpp"

class MatVec
{
public:
    MatVec(Complex c, String str);
    void operator()(const Complex *in, Complex *out) const;

private:
    Complex coeff;
    String string;
    Mask flip_mask, y_mask, z_mask;
};
