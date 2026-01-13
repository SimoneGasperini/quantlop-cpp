#pragma once

#include "utils.hpp"

class StateVector
{
public:
    StateVector(int num_qubits, Vector data);

    StateVector(const StateVector& other);

    int num_qubits() const;
    Vector &data();
    const Vector &data() const;

    Complex &operator[](int idx);
    const Complex &operator[](int idx) const;

    double inf_norm() const;

private:
    int _num_qubits;
    Vector _data;
};
