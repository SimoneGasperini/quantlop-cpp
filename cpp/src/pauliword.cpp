#include "pauliword.hpp"

PauliWord::PauliWord(Complex c, String str)
    : coeff(c),
      string(str),
      matvec(c, str) {}

Size PauliWord::num_qubits() const { return string.size(); }

PauliWord PauliWord::operator*(Complex c) const { return PauliWord(coeff * c, string); }

PauliWord operator*(Complex c, const PauliWord &pw) { return pw * c; }
