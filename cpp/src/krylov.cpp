#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

#include "hamiltonian.hpp"
#include "simulation.hpp"

namespace
{
    double l2_norm(const Complex *arr, Size dim)
    {
        double sum = 0.0;
        for (Index i = 0; i < dim; ++i)
        {
            sum += std::norm(arr[i]);
        }
        return std::sqrt(sum);
    }

    Complex dot_product(const Complex *lhs, const Complex *rhs, Size dim)
    {
        Complex out = 0.0;
        for (Index i = 0; i < dim; ++i)
        {
            out += std::conj(lhs[i]) * rhs[i];
        }
        return out;
    }

    double one_norm_dense(const std::vector<Complex> &arr, Size dim)
    {
        double best = 0.0;
        for (Index col = 0; col < dim; ++col)
        {
            double col_sum = 0.0;
            for (Index row = 0; row < dim; ++row)
            {
                col_sum += std::abs(arr[row * dim + col]);
            }
            best = std::max(best, col_sum);
        }
        return best;
    }

    std::vector<Complex> identity_dense(Size n)
    {
        std::vector<Complex> out(n * n, Complex(0.0, 0.0));
        for (Index i = 0; i < n; ++i)
        {
            out[i * n + i] = Complex(1.0, 0.0);
        }
        return out;
    }

    std::vector<Complex> matmul_dense(const std::vector<Complex> &a, const std::vector<Complex> &b, Size n)
    {
        std::vector<Complex> c(n * n, Complex(0.0, 0.0));
        for (Index i = 0; i < n; ++i)
        {
            for (Index k = 0; k < n; ++k)
            {
                const Complex aik = a[i * n + k];
                if (aik == Complex(0.0, 0.0))
                {
                    continue;
                }
                for (Index j = 0; j < n; ++j)
                {
                    c[i * n + j] += aik * b[k * n + j];
                }
            }
        }
        return c;
    }

    std::vector<Complex> solve_dense_system(std::vector<Complex> a, std::vector<Complex> b, Size n)
    {
        for (Index col = 0; col < n; ++col)
        {
            Index pivot = col;
            double pivot_abs = std::abs(a[col * n + col]);
            for (Index row = col + 1; row < n; ++row)
            {
                const double cand = std::abs(a[row * n + col]);
                if (cand > pivot_abs)
                {
                    pivot = row;
                    pivot_abs = cand;
                }
            }

            if (pivot != col)
            {
                for (Index j = 0; j < n; ++j)
                {
                    std::swap(a[col * n + j], a[pivot * n + j]);
                    std::swap(b[col * n + j], b[pivot * n + j]);
                }
            }

            const Complex diag = a[col * n + col];
            for (Index j = 0; j < n; ++j)
            {
                a[col * n + j] /= diag;
                b[col * n + j] /= diag;
            }

            for (Index row = 0; row < n; ++row)
            {
                if (row == col)
                {
                    continue;
                }
                const Complex factor = a[row * n + col];
                if (factor == Complex(0.0, 0.0))
                {
                    continue;
                }
                for (Index j = 0; j < n; ++j)
                {
                    a[row * n + j] -= factor * a[col * n + j];
                    b[row * n + j] -= factor * b[col * n + j];
                }
            }
        }
        return b;
    }

    std::vector<Complex> expm_dense(std::vector<Complex> a, Size n)
    {
        if (n == 1)
        {
            return {std::exp(a[0])};
        }

        constexpr std::array<double, 14> b = {
            64764752532480000.0,
            32382376266240000.0,
            7771770303897600.0,
            1187353796428800.0,
            129060195264000.0,
            10559470521600.0,
            670442572800.0,
            33522128640.0,
            1323241920.0,
            40840800.0,
            960960.0,
            16380.0,
            182.0,
            1.0};

        constexpr double theta_13 = 5.371920351148152;
        const double norm_a = one_norm_dense(a, n);
        int s = 0;
        if (norm_a > theta_13)
        {
            s = static_cast<int>(std::ceil(std::log2(norm_a / theta_13)));
        }
        const double scale = std::ldexp(1.0, s);
        for (Complex &z : a)
        {
            z /= scale;
        }

        const std::vector<Complex> I = identity_dense(n);
        const std::vector<Complex> A2 = matmul_dense(a, a, n);
        const std::vector<Complex> A4 = matmul_dense(A2, A2, n);
        const std::vector<Complex> A6 = matmul_dense(A4, A2, n);

        std::vector<Complex> tmp1(n * n, Complex(0.0, 0.0));
        std::vector<Complex> tmp2(n * n, Complex(0.0, 0.0));
        std::vector<Complex> U_inner(n * n, Complex(0.0, 0.0));
        std::vector<Complex> V(n * n, Complex(0.0, 0.0));

        for (Index i = 0; i < n * n; ++i)
        {
            tmp1[i] = b[13] * A6[i] + b[11] * A4[i] + b[9] * A2[i];
        }
        tmp2 = matmul_dense(A6, tmp1, n);
        for (Index i = 0; i < n * n; ++i)
        {
            U_inner[i] = tmp2[i] + b[7] * A6[i] + b[5] * A4[i] + b[3] * A2[i] + b[1] * I[i];
        }
        const std::vector<Complex> U = matmul_dense(a, U_inner, n);

        for (Index i = 0; i < n * n; ++i)
        {
            tmp1[i] = b[12] * A6[i] + b[10] * A4[i] + b[8] * A2[i];
        }
        tmp2 = matmul_dense(A6, tmp1, n);
        for (Index i = 0; i < n * n; ++i)
        {
            V[i] = tmp2[i] + b[6] * A6[i] + b[4] * A4[i] + b[2] * A2[i] + b[0] * I[i];
        }

        std::vector<Complex> P(n * n), Q(n * n);
        for (Index i = 0; i < n * n; ++i)
        {
            P[i] = V[i] - U[i];
            Q[i] = V[i] + U[i];
        }
        std::vector<Complex> R = solve_dense_system(std::move(P), std::move(Q), n);

        for (int k = 0; k < s; ++k)
        {
            R = matmul_dense(R, R, n);
        }
        return R;
    }

    std::vector<Complex> extract_scaled_dense(
        const std::vector<Complex> &h,
        Size ld,
        Size n,
        Complex scale)
    {
        std::vector<Complex> out(n * n, Complex(0.0, 0.0));
        for (Index row = 0; row < n; ++row)
        {
            for (Index col = 0; col < n; ++col)
            {
                out[row * n + col] = scale * h[row * ld + col];
            }
        }
        return out;
    }

    Size build_lanczos_tridiagonal(
        const Hamiltonian &ham,
        const Complex *psi,
        double bnorm,
        std::vector<Complex> &T,
        Size m)
    {
        const Size dim = Size(1) << ham.num_qubits();
        const double norm_tol = std::numeric_limits<double>::epsilon() * 1e2;
        std::vector<Complex> prev(dim, Complex(0.0, 0.0));
        std::vector<Complex> curr(dim, Complex(0.0, 0.0));
        std::vector<Complex> w(dim, Complex(0.0, 0.0));
        double beta_prev = 0.0;

        for (Index row = 0; row < dim; ++row)
        {
            curr[row] = psi[row] / bnorm;
        }

        for (Index k = 0; k < m; ++k)
        {
            ham.matvec_into(curr.data(), w.data());

            if (k > 0)
            {
                for (Index row = 0; row < dim; ++row)
                {
                    w[row] -= beta_prev * prev[row];
                }
            }

            Complex alpha = dot_product(curr.data(), w.data(), dim);
            alpha = Complex(alpha.real(), 0.0);
            T[k * m + k] = alpha;

            for (Index row = 0; row < dim; ++row)
            {
                w[row] -= alpha * curr[row];
            }

            const double beta_next = l2_norm(w.data(), dim);
            if (k + 1 < m)
            {
                T[(k + 1) * m + k] = beta_next;
                T[k * m + (k + 1)] = beta_next;
            }

            if (beta_next < norm_tol || k + 1 == m)
            {
                return k + 1;
            }

            std::copy(curr.begin(), curr.end(), prev.begin());
            for (Index row = 0; row < dim; ++row)
            {
                curr[row] = w[row] / beta_next;
            }

            beta_prev = beta_next;
        }
        return m;
    }

    void reconstruct_lanczos_state(
        const Hamiltonian &ham,
        const Complex *psi,
        double bnorm,
        const std::vector<Complex> &T,
        Size basis_size,
        Complex coeff,
        Complex *out)
    {
        const Size dim = Size(1) << ham.num_qubits();
        std::vector<Complex> prev(dim, Complex(0.0, 0.0));
        std::vector<Complex> curr(dim, Complex(0.0, 0.0));
        std::vector<Complex> w(dim, Complex(0.0, 0.0));

        for (Index row = 0; row < dim; ++row)
        {
            curr[row] = psi[row] / bnorm;
        }

        const Complex scale = Complex(0.0, -1.0) * coeff;
        const std::vector<Complex> fT = expm_dense(extract_scaled_dense(T, basis_size, basis_size, scale), basis_size);

        for (Index row = 0; row < dim; ++row)
        {
            out[row] = bnorm * curr[row] * fT[0];
        }

        for (Index k = 1; k < basis_size; ++k)
        {
            const double beta_prev = T[(k - 1) * basis_size + k].real();
            const Complex alpha_prev = T[(k - 1) * basis_size + (k - 1)];
            ham.matvec_into(curr.data(), w.data());

            for (Index row = 0; row < dim; ++row)
            {
                w[row] -= alpha_prev * curr[row];
            }
            if (k > 1)
            {
                const double beta_prevprev = T[(k - 2) * basis_size + (k - 1)].real();
                for (Index row = 0; row < dim; ++row)
                {
                    w[row] -= beta_prevprev * prev[row];
                }
            }

            std::copy(curr.begin(), curr.end(), prev.begin());
            for (Index row = 0; row < dim; ++row)
            {
                curr[row] = w[row] / beta_prev;
                out[row] += bnorm * curr[row] * fT[k * basis_size + 0];
            }
        }
    }
}

Complex *expm_multiply_krylov(const Hamiltonian &ham, const Complex *psi, Complex coeff)
{
    const Size dim = Size(1) << ham.num_qubits();
    Complex *out = new Complex[dim];
    std::fill(out, out + dim, Complex(0.0, 0.0));

    const double bnorm = l2_norm(psi, dim);
    const Size m = std::min<Size>(30, dim);

    std::vector<Complex> T(m * m, Complex(0.0, 0.0));
    const Size basis_size = build_lanczos_tridiagonal(ham, psi, bnorm, T, m);
    reconstruct_lanczos_state(ham, psi, bnorm, T, basis_size, coeff, out);

    return out;
}
