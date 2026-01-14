#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include "hamiltonian.hpp"
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

inline std::pair<int, std::size_t> fragment_3_1(const Hamiltonian &ham)
{
    int m_star = -1;
    std::size_t s_star = 0;
    double best_cost = 0.0;
    const double one_norm = ham.lcu_norm();
    for (const auto &entry : theta_table())
    {
        const int m = entry.first;
        const double theta = entry.second;
        const double s_m = std::ceil(one_norm / theta);
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

inline double inf_norm(const Complex *psi, std::size_t dim)
{
    double norm = 0.0;
    for (std::size_t i = 0; i < dim; ++i)
    {
        const double val = std::abs(psi[i]);
        if (val > norm)
        {
            norm = val;
        }
    }
    return norm;
}

inline void expm_multiply_into(const Hamiltonian &ham, const Complex *psi, std::size_t dim,
                               Complex *out)
{
    const double tol = std::ldexp(1.0, -24);
    const auto [m_star, s] = fragment_3_1(ham);
    std::vector<Complex> b(psi, psi + dim);
    std::vector<Complex> tmp(dim);
    std::copy(psi, psi + dim, out);
    for (std::size_t step = 0; step < s; ++step)
    {
        double c1 = inf_norm(b.data(), dim);
        for (int j = 1; j <= m_star; ++j)
        {
            ham.matvec_into(b.data(), tmp.data());
            const double scale = 1.0 / static_cast<double>(s * j);
            for (auto &v : tmp)
            {
                v *= scale;
            }
            b.swap(tmp);
            double c2 = inf_norm(b.data(), dim);
            for (std::size_t i = 0; i < dim; ++i)
            {
                out[i] += b[i];
            }
            if (c1 + c2 <= tol * inf_norm(out, dim))
            {
                break;
            }
            c1 = c2;
        }
        std::copy(out, out + dim, b.begin());
    }
}

inline std::size_t expm_multiply_trace_into(const Hamiltonian &ham, const Complex *psi,
                                            std::size_t dim, Complex *out)
{
    const double tol = std::ldexp(1.0, -16);
    const auto [m_star, s] = fragment_3_1(ham);
    std::vector<Complex> b(psi, psi + dim);
    std::vector<Complex> f(psi, psi + dim);
    std::vector<Complex> tmp(dim);
    std::copy(f.begin(), f.end(), out);
    std::size_t trace_len = 1;
    for (std::size_t step = 0; step < s; ++step)
    {
        double c1 = inf_norm(b.data(), dim);
        for (int j = 1; j <= m_star; ++j)
        {
            ham.matvec_into(b.data(), tmp.data());
            const double scale = 1.0 / static_cast<double>(s * j);
            for (auto &v : tmp)
            {
                v *= scale;
            }
            b.swap(tmp);
            double c2 = inf_norm(b.data(), dim);
            for (std::size_t i = 0; i < dim; ++i)
            {
                f[i] += b[i];
            }
            if (c1 + c2 <= tol * inf_norm(f.data(), dim))
            {
                break;
            }
            c1 = c2;

            std::copy(f.begin(), f.end(), out + trace_len * dim);
            ++trace_len;
        }
        std::copy(f.begin(), f.end(), b.begin());
    }
    return trace_len;
}

inline std::size_t trace_evolve_max_steps(const Hamiltonian &ham)
{
    const auto [m_star, s] = fragment_3_1(ham);
    if (m_star < 0)
    {
        return 0;
    }
    return 1 + s * static_cast<std::size_t>(m_star);
}

inline void evolve_into(const Hamiltonian &ham, const Complex *psi, std::size_t dim,
                        Complex *out, Complex coeff = Complex(1.0, 0.0))
{
    const Complex i(0.0, 1.0);
    const Hamiltonian expm = (-i * coeff) * ham;
    expm_multiply_into(expm, psi, dim, out);
}

inline std::size_t trace_evolve_into(const Hamiltonian &ham, const Complex *psi, std::size_t dim,
                                     Complex *out, Complex coeff = Complex(1.0, 0.0))
{
    const Complex i(0.0, 1.0);
    const Hamiltonian expm = (-i * coeff) * ham;
    return expm_multiply_trace_into(expm, psi, dim, out);
}
