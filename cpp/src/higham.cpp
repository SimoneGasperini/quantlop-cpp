#include "simulation.hpp"

namespace
{
    const std::vector<std::pair<int, double>> &theta_table()
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

    std::pair<int, int> fragment_3_1(const Hamiltonian &ham)
    {
        int m_star = -1;
        int s_star = 0;
        double best_cost = 0.0;
        const double one_norm = ham.lcu_norm();
        for (const auto &pair : theta_table())
        {
            const int m = pair.first;
            const double theta = pair.second;
            const double s_m = std::ceil(one_norm / theta);
            const double cost = m * s_m;
            if (m_star < 0 || cost < best_cost)
            {
                m_star = m;
                s_star = static_cast<int>(s_m);
                best_cost = cost;
            }
        }
        return std::make_pair(m_star, s_star);
    }

    double inf_norm(const Complex *psi, Size dim)
    {
        double norm = 0.0;
        for (Index i = 0; i < dim; ++i)
        {
            const double val = std::abs(psi[i]);
            if (val > norm)
            {
                norm = val;
            }
        }
        return norm;
    }
}

Complex *expm_multiply_higham(const Hamiltonian &ham, const Complex *psi)
{
    const double tol = std::ldexp(1.0, -24);
    const auto [m_star, s] = fragment_3_1(ham);
    const Size dim = Size(1) << ham.num_qubits();

    Complex *b = new Complex[dim];
    Complex *tmp = new Complex[dim];
    Complex *out = new Complex[dim];

    std::copy(psi, psi + dim, b);
    std::copy(psi, psi + dim, out);

    const double inv_s = 1.0 / static_cast<double>(s);
    for (Index step = 0; step < static_cast<Size>(s); ++step)
    {
        double c1 = inf_norm(b, dim);
        for (int j = 1; j <= m_star; ++j)
        {
            ham.matvec_into(b, tmp);
            const double scale = inv_s / static_cast<double>(j);
            for (Index k = 0; k < dim; ++k)
            {
                tmp[k] *= scale;
            }
            std::swap(b, tmp);
            const double c2 = inf_norm(b, dim);
            for (Index i = 0; i < dim; ++i)
            {
                out[i] += b[i];
            }
            if (c1 + c2 <= tol * inf_norm(out, dim))
            {
                break;
            }
            c1 = c2;
        }
        std::copy(out, out + dim, b);
    }

    delete[] tmp;
    delete[] b;
    return out;
}
