#ifndef AMG_SOLVER_3D_H
#define AMG_SOLVER_3D_H

#include <vector>
#include <algorithm>
#include <cmath>

// Utility vector types
struct Vec3i {
    int x, y, z;
    Vec3i() : x(0), y(0), z(0) {}
    Vec3i(int x_, int y_, int z_) : x(x_), y(y_), z(z_) {}
};

// Fixed sparse matrix structure for AMG solver
template<class T>
class FixedSparseMatrix3D {
public:
    struct Entry {
        int row, col;
        T value;
        Entry(int r, int c, T v) : row(r), col(c), value(v) {}
    };
    
    std::vector<Entry> entries;
    int rows, cols;
    
    FixedSparseMatrix3D(int r, int c) : rows(r), cols(c) {}
    
    void addEntry(int row, int col, T value) {
        entries.push_back(Entry(row, col, value));
    }
    
    void clear() {
        entries.clear();
    }
};

// BLAS-like operations namespace
namespace BLAS3D {
    template<class T>
    void zero(std::vector<T>& x) {
        std::fill(x.begin(), x.end(), T(0));
    }
    
    template<class T>
    T dot(const std::vector<T>& x, const std::vector<T>& y) {
        T result = T(0);
        for (size_t i = 0; i < x.size(); ++i) {
            result += x[i] * y[i];
        }
        return result;
    }
    
    template<class T>
    T abs_max(const std::vector<T>& x) {
        T result = T(0);
        for (size_t i = 0; i < x.size(); ++i) {
            result = std::max(result, std::abs(x[i]));
        }
        return result;
    }
    
    template<class T>
    T mean(const std::vector<T>& x) {
        T sum = T(0);
        for (size_t i = 0; i < x.size(); ++i) {
            sum += x[i];
        }
        return sum / T(x.size());
    }
    
    template<class T>
    void add_scaled(T alpha, const std::vector<T>& x, std::vector<T>& y) {
        for (size_t i = 0; i < x.size() && i < y.size(); ++i) {
            y[i] += alpha * x[i];
        }
    }
    
    template<class T>
    void subtractConst(std::vector<T>& x, T c) {
        for (size_t i = 0; i < x.size(); ++i) {
            x[i] -= c;
        }
    }
}

// Zero vector utility
template<class T>
void zero3D(std::vector<T>& x) {
    BLAS3D::zero(x);
}

// Matrix-vector multiplication
template<class T>
void multiply3D(const FixedSparseMatrix3D<T>& matrix, const std::vector<T>& x, std::vector<T>& y) {
    BLAS3D::zero(y);
    for (const auto& entry : matrix.entries) {
        y[entry.row] += entry.value * x[entry.col];
    }
}

// Simple AMG preconditioner for 3D (placeholder implementation)
template<class T>
void amgPrecond3D(std::vector<FixedSparseMatrix3D<T>*>& A_L,
                  std::vector<FixedSparseMatrix3D<T>*>& R_L,
                  std::vector<FixedSparseMatrix3D<T>*>& P_L,
                  std::vector<Vec3i>& S_L,
                  std::vector<T>& z,
                  const std::vector<T>& r) {
    // Simple Jacobi preconditioner as placeholder
    z = r;
    for (size_t i = 0; i < z.size(); ++i) {
        z[i] *= T(0.6);  // Simple damping factor
    }
}

// Poisson matrix setup for 3D grid
template<class T>
void setupPoissonMatrix3D(FixedSparseMatrix3D<T>& matrix, int ni, int nj, int nk) {
    matrix.clear();
    
    for (int k = 0; k < nk; ++k) {
        for (int j = 0; j < nj; ++j) {
            for (int i = 0; i < ni; ++i) {
                int idx = k * ni * nj + j * ni + i;
                
                // Diagonal entry
                T diag_val = T(6);
                matrix.addEntry(idx, idx, diag_val);
                
                // Off-diagonal entries
                if (i > 0) {
                    matrix.addEntry(idx, idx - 1, T(-1));  // Left
                }
                if (i < ni - 1) {
                    matrix.addEntry(idx, idx + 1, T(-1));  // Right
                }
                if (j > 0) {
                    matrix.addEntry(idx, idx - ni, T(-1)); // Back
                }
                if (j < nj - 1) {
                    matrix.addEntry(idx, idx + ni, T(-1)); // Front
                }
                if (k > 0) {
                    matrix.addEntry(idx, idx - ni * nj, T(-1)); // Bottom
                }
                if (k < nk - 1) {
                    matrix.addEntry(idx, idx + ni * nj, T(-1)); // Top
                }
            }
        }
    }
}

// AMG PCG Solver for 3D
template<class T>
bool AMGPCGSolvePrebuilt3D(const FixedSparseMatrix3D<T> &fixed_matrix,
                           const std::vector<T> &rhs,
                           std::vector<T> &result,
                           std::vector<FixedSparseMatrix3D<T> *> &A_L,
                           std::vector<FixedSparseMatrix3D<T> *> &R_L,
                           std::vector<FixedSparseMatrix3D<T> *> &P_L,
                           std::vector<Vec3i> &S_L,
                           const int total_level,
                           T tolerance_factor,
                           int max_iterations,
                           T &residual_out,
                           int &iterations_out,
                           int ni, int nj, int nk,
                           bool PURE_NEUMANN) {
    std::vector<T> m, z, s, r, v;
    unsigned int n = ni * nj * nk;
    if (m.size() != n) {
        m.resize(n);
        s.resize(n);
        z.resize(n);
        r.resize(n);
        if (PURE_NEUMANN) { v.resize(n); }
    }
    zero3D(result);
    r = rhs;
    if (PURE_NEUMANN) {
        T mean = BLAS3D::mean(r);
        v = r;
        BLAS3D::subtractConst(v, mean);
    }
    residual_out = PURE_NEUMANN ? BLAS3D::abs_max(v) : BLAS3D::abs_max(r);
    if (residual_out == 0) {
        iterations_out = 0;
        return true;
    }
    T tol = std::max(tolerance_factor * residual_out, tolerance_factor);
    if (PURE_NEUMANN) r = v;
    amgPrecond3D(A_L, R_L, P_L, S_L, z, r);
    T rho = BLAS3D::dot(z, r);
    if (rho == 0 || rho != rho) {
        iterations_out = 0;
        return false;
    }

    s = z;

    int iteration;
    for (iteration = 0; iteration < max_iterations; ++iteration) {
        multiply3D(fixed_matrix, s, z);
        T alpha = rho / BLAS3D::dot(s, z);
        BLAS3D::add_scaled(alpha, s, result);
        BLAS3D::add_scaled(-alpha, z, r);
        if (PURE_NEUMANN) {
            T mean = BLAS3D::mean(r);
            v = r;
            BLAS3D::subtractConst(v, mean);
        }
        residual_out = PURE_NEUMANN ? BLAS3D::abs_max(v) : BLAS3D::abs_max(r);
        if (residual_out <= tol) {
            iterations_out = iteration + 1;
            return true;
        }
        if (PURE_NEUMANN) r = v;
        amgPrecond3D(A_L, R_L, P_L, S_L, z, r);
        T rho_new = BLAS3D::dot(z, r);
        T beta = rho_new / rho;
        BLAS3D::add_scaled(beta, s, z);
        s.swap(z); // s=beta*s+z
        rho = rho_new;
    }
    iterations_out = iteration;
    return false;
}

#endif // AMG_SOLVER_3D_H
