#include "simulation.hpp"

namespace
{
    double l2_norm(const Complex *vec, Size dim)
    {
        double sum = 0.0;
        for (Index i = 0; i < dim; ++i)
        {
            sum += std::norm(vec[i]);
        }
        return std::sqrt(sum);
    }

    double l2_norm(const std::vector<Complex> &vec)
    {
        double sum = 0.0;
        for (const Complex &z : vec)
        {
            sum += std::norm(z);
        }
        return std::sqrt(sum);
    }

    double one_norm_dense(const std::vector<Complex> &a, Size n)
    {
        double best = 0.0;
        for (Index col = 0; col < n; ++col)
        {
            double col_sum = 0.0;
            for (Index row = 0; row < n; ++row)
            {
                col_sum += std::abs(a[row * n + col]);
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

            if (pivot_abs == 0.0)
            {
                throw std::runtime_error("Singular matrix encountered in Krylov matrix exponential");
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
        if (n == 0)
        {
            return {};
        }
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

    bool krylov_arnoldi_iteration(
        const Hamiltonian &ham,
        const Complex *b,
        double bnorm,
        std::vector<Complex> &V,
        std::vector<Complex> &H,
        Size h_ld,
        Size begin,
        Size m)
    {
        const Size dim = Size(1) << ham.num_qubits();
        const double norm_tol = std::numeric_limits<double>::epsilon() * 1e2;
        std::vector<Complex> w(dim, Complex(0.0, 0.0));

        for (Index i = 0; i < dim; ++i)
        {
            V[i] = b[i] / bnorm;
        }

        for (Index k = 0; k < m; ++k)
        {
            ham.matvec_into(V.data() + k * dim, w.data());
            std::copy(w.begin(), w.end(), V.begin() + (k + 1) * dim);

            for (Index i = 0; i <= k; ++i)
            {
                Complex hij = 0.0;
                const Complex *vi = V.data() + i * dim;
                Complex *vkp1 = V.data() + (k + 1) * dim;
                for (Index row = 0; row < dim; ++row)
                {
                    hij += std::conj(vi[row]) * vkp1[row];
                }
                H[(begin + i) * h_ld + (begin + k)] = hij;
                for (Index row = 0; row < dim; ++row)
                {
                    vkp1[row] -= hij * vi[row];
                }
            }

            Complex *vkp1 = V.data() + (k + 1) * dim;
            const double beta = l2_norm(vkp1, dim);
            H[(begin + k + 1) * h_ld + (begin + k)] = beta;
            if (beta < norm_tol)
            {
                return true;
            }
            for (Index row = 0; row < dim; ++row)
            {
                vkp1[row] /= beta;
            }
        }

        return false;
    }
}

Complex *expm_multiply_krylov(const Hamiltonian &ham, const Complex *psi)
{
    const Size dim = Size(1) << ham.num_qubits();
    Complex *y = new Complex[dim];
    std::fill(y, y + dim, Complex(0.0, 0.0));

    const Size m = std::min<Size>(10, dim);
    if (m == 0)
    {
        return y;
    }

    const int max_restarts_base = 3;
    const int max_restarts = std::min(max_restarts_base, static_cast<int>(dim / m) + 1);
    const Size mmax = m * static_cast<Size>(max_restarts);
    const double bnorm = l2_norm(psi, dim);
    if (bnorm == 0.0)
    {
        return y;
    }

    const double rtol = 1e-3;
    const double atol = rtol * bnorm;
    std::vector<Complex> V(dim * (m + 1), Complex(0.0, 0.0));
    std::vector<Complex> H((mmax + 1) * mmax, Complex(0.0, 0.0));

    const bool breakdown0 = krylov_arnoldi_iteration(ham, psi, bnorm, V, H, mmax, 0, m);
    Size j0 = m;
    if (breakdown0)
    {
        for (Index k = 0; k < m; ++k)
        {
            if (std::abs(H[(k + 1) * mmax + k]) < std::numeric_limits<double>::epsilon() * 1e2)
            {
                j0 = k + 1;
                break;
            }
        }
    }

    auto extract_h = [&H, mmax](Size k)
    {
        std::vector<Complex> hk(k * k, Complex(0.0, 0.0));
        for (Index row = 0; row < k; ++row)
        {
            for (Index col = 0; col < k; ++col)
            {
                hk[row * k + col] = H[row * mmax + col];
            }
        }
        return hk;
    };

    std::vector<Complex> fH = expm_dense(extract_h(j0), j0);
    for (Index row = 0; row < dim; ++row)
    {
        Complex accum = 0.0;
        for (Index col = 0; col < j0; ++col)
        {
            accum += V[col * dim + row] * fH[col * j0 + 0];
        }
        y[row] = bnorm * accum;
    }

    if (breakdown0)
    {
        return y;
    }

    double update_norm = 0.0;
    for (Index row = 0; row < j0; ++row)
    {
        update_norm += std::norm(bnorm * fH[row * j0 + 0]);
    }
    update_norm = std::sqrt(update_norm);

    int restart = 1;
    while (restart < max_restarts && update_norm > atol)
    {
        const Size begin = static_cast<Size>(restart) * m;
        const Size end = static_cast<Size>(restart + 1) * m;
        const bool breakdown = krylov_arnoldi_iteration(ham, V.data() + m * dim, 1.0, V, H, mmax, begin, m);

        Size j = m;
        if (breakdown)
        {
            for (Index k = 0; k < m; ++k)
            {
                if (std::abs(H[(begin + k + 1) * mmax + (begin + k)]) < std::numeric_limits<double>::epsilon() * 1e2)
                {
                    j = k + 1;
                    break;
                }
            }
        }

        const Size k = breakdown ? (begin + j) : end;
        fH = expm_dense(extract_h(k), k);

        const Size cols = breakdown ? j : m;
        for (Index row = 0; row < dim; ++row)
        {
            Complex accum = 0.0;
            for (Index col = 0; col < cols; ++col)
            {
                const Size global_row = begin + col;
                accum += V[col * dim + row] * fH[global_row * k + 0];
            }
            y[row] += bnorm * accum;
        }

        if (breakdown)
        {
            return y;
        }

        update_norm = 0.0;
        for (Index row = begin; row < end; ++row)
        {
            update_norm += std::norm(bnorm * fH[row * k + 0]);
        }
        update_norm = std::sqrt(update_norm);
        ++restart;
    }
    return y;
}
