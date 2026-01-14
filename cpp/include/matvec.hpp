#pragma once

#include <complex>
#include <string>

class MatVec
{
public:
    MatVec(std::complex<double> coeff, std::string string);

    void operator()(const std::complex<double> *in, std::complex<double> *out) const;

private:
    std::complex<double> _coeff;
    std::string _string;
};
