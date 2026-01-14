#include "statevector.hpp"

StateVector::StateVector(std::size_t num_qubits, Vector data)
    : _num_qubits(num_qubits),
      _data(std::move(data)) {}

StateVector::StateVector(const StateVector &other)
    : _num_qubits(other._num_qubits),
      _data(other._data) {}

std::size_t StateVector::num_qubits() const { return _num_qubits; }

Vector &StateVector::data() { return _data; }

const Vector &StateVector::data() const { return _data; }

Complex &StateVector::operator[](std::size_t idx) { return _data[idx]; }

const Complex &StateVector::operator[](std::size_t idx) const { return _data[idx]; }

double StateVector::inf_norm() const
{
    double norm = 0.0;
    for (const Complex &c : _data)
    {
        double val = std::abs(c);
        if (val > norm)
        {
            norm = val;
        }
    }
    return norm;
}
