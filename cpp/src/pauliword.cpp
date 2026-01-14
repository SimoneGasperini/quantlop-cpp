#include "pauliword.hpp"

#include <cstddef>
#include <utility>

PauliWord::PauliWord(std::complex<double> coeff, std::string string)
    : _coeff(coeff),
      _string(std::move(string)),
      _matvec(_coeff, _string) {}

std::complex<double> PauliWord::coeff() const { return _coeff; }

const std::string &PauliWord::string() const { return _string; }

std::size_t PauliWord::num_qubits() const { return _string.size(); }

PauliWord PauliWord::operator*(std::complex<double> c) const { return PauliWord(_coeff * c, _string); }

PauliWord operator*(std::complex<double> c, const PauliWord &pw) { return pw * c; }
