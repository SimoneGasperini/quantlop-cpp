#pragma once

#include <complex>
#include <vector>
#include <functional>

using Complex = std::complex<double>;
using Vector = std::vector<Complex>;
using MatVecFn = std::function<Vector(const Vector &)>;
