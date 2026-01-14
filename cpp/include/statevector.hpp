#pragma once

#include <cstddef>
#include "utils.hpp"

class StateVector
{
public:
    StateVector(std::size_t num_qubits, Vector data);

    StateVector(const StateVector &other);

    std::size_t num_qubits() const;
    Vector &data();
    const Vector &data() const;

    Complex &operator[](std::size_t idx);
    const Complex &operator[](std::size_t idx) const;

    double inf_norm() const;

private:
    std::size_t _num_qubits;
    Vector _data;
};
