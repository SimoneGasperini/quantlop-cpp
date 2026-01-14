#pragma once

#include <complex>
#include <functional>

using Complex = std::complex<double>;
using MatVecFn = std::function<void(const Complex *, Complex *)>;
