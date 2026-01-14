#pragma once

#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include "hamiltonian.hpp"
#include "statevector.hpp"
#include "utils.hpp"

inline const std::vector<std::pair<int, double>> &theta_table()
{
    static const std::vector<std::pair<int, double>> table = {
        {5, 1.3e-1},
        {10, 1.0e0},
        {15, 2.2e0},
        {20, 3.6e0},
        {25, 4.9e0},
        {30, 6.3e0},
        {35, 7.7e0},
        {40, 9.1e0},
        {45, 1.1e1},
        {50, 1.2e1},
        {55, 1.3e1},
    };
    return table;
}

inline std::pair<int, std::size_t> fragment_3_1(const Hamiltonian &A)
{
    int m_star = -1;
    std::size_t s_star = 0;
    double best_cost = 0.0;
    const double A_one_norm = A.lcu_norm();
    for (const auto &entry : theta_table())
    {
        const int m = entry.first;
        const double theta = entry.second;
        const double s_m = std::ceil(A_one_norm / theta);
        const double cost = m * s_m;
        if (m_star < 0 || cost < best_cost)
        {
            m_star = m;
            s_star = static_cast<std::size_t>(s_m);
            best_cost = cost;
        }
    }
    return std::make_pair(m_star, s_star);
}

inline StateVector expm_multiply(const Hamiltonian &A, const StateVector &psi)
{
    const double tol = std::ldexp(1.0, -24);
    const auto [m_star, s] = fragment_3_1(A);
    StateVector b = psi;
    StateVector f = psi;
    auto &f_data = f.data();
    for (std::size_t step = 0; step < s; ++step)
    {
        double c1 = b.inf_norm();
        for (int j = 1; j <= m_star; ++j)
        {
            Vector tmp = A.matvec(b.data());
            const double scale = 1.0 / static_cast<double>(s * j);
            for (auto &v : tmp)
            {
                v *= scale;
            }
            b = StateVector(psi.num_qubits(), std::move(tmp));
            double c2 = b.inf_norm();
            const auto &b_data = b.data();
            for (std::size_t i = 0; i < f_data.size(); ++i)
            {
                f_data[i] += b_data[i];
            }
            if (c1 + c2 <= tol * f.inf_norm())
            {
                break;
            }
            c1 = c2;
        }
        b = f;
    }
    return f;
}

inline std::vector<StateVector> expm_multiply_trace(const Hamiltonian &A, const StateVector &psi)
{
    std::vector<StateVector> trace;
    const double tol = std::ldexp(1.0, -16);
    const auto [m_star, s] = fragment_3_1(A);
    StateVector b = psi;
    StateVector f = psi;
    trace.push_back(StateVector(f));
    auto &f_data = f.data();
    for (std::size_t step = 0; step < s; ++step)
    {
        double c1 = b.inf_norm();
        for (int j = 1; j <= m_star; ++j)
        {
            Vector tmp = A.matvec(b.data());
            const double scale = 1.0 / static_cast<double>(s * j);
            for (auto &v : tmp)
            {
                v *= scale;
            }
            b = StateVector(psi.num_qubits(), std::move(tmp));
            double c2 = b.inf_norm();
            const auto &b_data = b.data();
            for (std::size_t i = 0; i < f_data.size(); ++i)
            {
                f_data[i] += b_data[i];
            }
            if (c1 + c2 <= tol * f.inf_norm())
            {
                break;
            }
            c1 = c2;

            trace.push_back(StateVector(f));
        }
        b = f;
    }
    return trace;
}

inline StateVector evolve(const Hamiltonian &ham, const StateVector &psi, Complex coeff = Complex(1.0, 0.0))
{
    const Complex i(0.0, 1.0);
    const Hamiltonian expm = (-i * coeff) * ham;
    return expm_multiply(expm, psi);
}

inline std::vector<StateVector> trace_evolve(const Hamiltonian &ham, const StateVector &psi, Complex coeff = Complex(1.0, 0.0))
{
    const Complex i(0.0, 1.0);
    const Hamiltonian expm = (-i * coeff) * ham;
    return expm_multiply_trace(expm, psi);
}
