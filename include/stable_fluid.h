#ifndef STABLE_FLUID_3D_H
#define STABLE_FLUID_3D_H

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

/**
 * @brief Data aggregation of current and next frame
 * @tparam T the type of data stored in aggregation
 */
template <typename T> struct TexPair
{
  public:
    T &cur;
    T &nxt;
    TexPair(T &a, T &b) : cur(a), nxt(b) {}
    void Swap()
    {
        T tmp = cur;
        cur = nxt;
        nxt = tmp;
    }
};

template <typename T> inline T min(T a, T b)
{
    return a < b ? a : b;
}
template <typename T> inline T max(T a, T b)
{
    return a > b ? a : b;
}

/**
 * @brief Get index of 1d array out of 3d index (i,j,k)
 *
 * @param i index of x-axis
 * @param j index of y-axis
 * @param k index of z-axis
 * @param N width of 3d array
 * @param M height of 3d array
 */
inline int IXYZ(int i, int j, int k, int N, int M)
{
    return k * N * M + j * N + i;
}

template <typename T, typename SCALAR> inline T lerp(T l, T r, SCALAR t)
{
    return l + t * (r - l);
}

/**
 * @brief Trilinear interpolation for 3D scalar field
 *
 * @param N width
 * @param M height
 * @param L depth
 * @param field scalar field
 * @param x position x
 * @param y position y
 * @param z position z
 */
template <typename T, typename SCALAR>
T trilerp_scalar(const int N, const int M, const int L, const std::vector<T> &field, 
                 SCALAR x, SCALAR y, SCALAR z)
{
    // Clamp coordinates to grid bounds
    x = max(SCALAR(0.5), min(SCALAR(N - 0.5), x));
    y = max(SCALAR(0.5), min(SCALAR(M - 0.5), y));
    z = max(SCALAR(0.5), min(SCALAR(L - 0.5), z));
    
    int i = int(x - 0.5);
    int j = int(y - 0.5);
    int k = int(z - 0.5);
    
    SCALAR fx = x - 0.5 - i;
    SCALAR fy = y - 0.5 - j;
    SCALAR fz = z - 0.5 - k;
    
    // Ensure indices are within bounds
    i = max(0, min(N - 2, i));
    j = max(0, min(M - 2, j));
    k = max(0, min(L - 2, k));
    
    // Trilinear interpolation
    T c000 = field[IXYZ(i,     j,     k,     N, M)];
    T c001 = field[IXYZ(i,     j,     k + 1, N, M)];
    T c010 = field[IXYZ(i,     j + 1, k,     N, M)];
    T c011 = field[IXYZ(i,     j + 1, k + 1, N, M)];
    T c100 = field[IXYZ(i + 1, j,     k,     N, M)];
    T c101 = field[IXYZ(i + 1, j,     k + 1, N, M)];
    T c110 = field[IXYZ(i + 1, j + 1, k,     N, M)];
    T c111 = field[IXYZ(i + 1, j + 1, k + 1, N, M)];
    
    T c00 = lerp(c000, c001, fz);
    T c01 = lerp(c010, c011, fz);
    T c10 = lerp(c100, c101, fz);
    T c11 = lerp(c110, c111, fz);
    
    T c0 = lerp(c00, c01, fy);
    T c1 = lerp(c10, c11, fy);
    
    return lerp(c0, c1, fx);
}

/**
 * @brief Trilinear interpolation for u component in MAC grid
 */
template <typename SCALAR>
SCALAR interpolate_u(const int N, const int M, const int L, const std::vector<SCALAR> &u_vel, 
                     SCALAR x, SCALAR y, SCALAR z)
{
    // u component is stored at (i, j+0.5, k+0.5)
    // Clamp coordinates
    x = max(SCALAR(0), min(SCALAR(N), x));
    y = max(SCALAR(0.5), min(SCALAR(M - 0.5), y));
    z = max(SCALAR(0.5), min(SCALAR(L - 0.5), z));
    
    int i = int(x);
    int j = int(y - 0.5);
    int k = int(z - 0.5);
    
    SCALAR fx = x - i;
    SCALAR fy = y - 0.5 - j;
    SCALAR fz = z - 0.5 - k;
    
    // Ensure indices are within bounds
    i = max(0, min(N - 1, i));
    j = max(0, min(M - 2, j));
    k = max(0, min(L - 2, k));
    
    // Trilinear interpolation
    SCALAR c000 = u_vel[IXYZ(i,     j,     k,     N + 1, M)];
    SCALAR c001 = u_vel[IXYZ(i,     j,     k + 1, N + 1, M)];
    SCALAR c010 = u_vel[IXYZ(i,     j + 1, k,     N + 1, M)];
    SCALAR c011 = u_vel[IXYZ(i,     j + 1, k + 1, N + 1, M)];
    SCALAR c100 = u_vel[IXYZ(i + 1, j,     k,     N + 1, M)];
    SCALAR c101 = u_vel[IXYZ(i + 1, j,     k + 1, N + 1, M)];
    SCALAR c110 = u_vel[IXYZ(i + 1, j + 1, k,     N + 1, M)];
    SCALAR c111 = u_vel[IXYZ(i + 1, j + 1, k + 1, N + 1, M)];
    
    SCALAR c00 = lerp(c000, c001, fz);
    SCALAR c01 = lerp(c010, c011, fz);
    SCALAR c10 = lerp(c100, c101, fz);
    SCALAR c11 = lerp(c110, c111, fz);
    
    SCALAR c0 = lerp(c00, c01, fy);
    SCALAR c1 = lerp(c10, c11, fy);
    
    return lerp(c0, c1, fx);
}

/**
 * @brief Trilinear interpolation for v component in MAC grid
 */
template <typename SCALAR>
SCALAR interpolate_v(const int N, const int M, const int L, const std::vector<SCALAR> &v_vel, 
                     SCALAR x, SCALAR y, SCALAR z)
{
    // v component is stored at (i+0.5, j, k+0.5)
    x = max(SCALAR(0.5), min(SCALAR(N - 0.5), x));
    y = max(SCALAR(0), min(SCALAR(M), y));
    z = max(SCALAR(0.5), min(SCALAR(L - 0.5), z));
    
    int i = int(x - 0.5);
    int j = int(y);
    int k = int(z - 0.5);
    
    SCALAR fx = x - 0.5 - i;
    SCALAR fy = y - j;
    SCALAR fz = z - 0.5 - k;
    
    i = max(0, min(N - 2, i));
    j = max(0, min(M - 1, j));
    k = max(0, min(L - 2, k));
    
    SCALAR c000 = v_vel[IXYZ(i,     j,     k,     N, M + 1)];
    SCALAR c001 = v_vel[IXYZ(i,     j,     k + 1, N, M + 1)];
    SCALAR c010 = v_vel[IXYZ(i,     j + 1, k,     N, M + 1)];
    SCALAR c011 = v_vel[IXYZ(i,     j + 1, k + 1, N, M + 1)];
    SCALAR c100 = v_vel[IXYZ(i + 1, j,     k,     N, M + 1)];
    SCALAR c101 = v_vel[IXYZ(i + 1, j,     k + 1, N, M + 1)];
    SCALAR c110 = v_vel[IXYZ(i + 1, j + 1, k,     N, M + 1)];
    SCALAR c111 = v_vel[IXYZ(i + 1, j + 1, k + 1, N, M + 1)];
    
    SCALAR c00 = lerp(c000, c001, fz);
    SCALAR c01 = lerp(c010, c011, fz);
    SCALAR c10 = lerp(c100, c101, fz);
    SCALAR c11 = lerp(c110, c111, fz);
    
    SCALAR c0 = lerp(c00, c01, fy);
    SCALAR c1 = lerp(c10, c11, fy);
    
    return lerp(c0, c1, fx);
}

/**
 * @brief Trilinear interpolation for w component in MAC grid
 */
template <typename SCALAR>
SCALAR interpolate_w(const int N, const int M, const int L, const std::vector<SCALAR> &w_vel, 
                     SCALAR x, SCALAR y, SCALAR z)
{
    // w component is stored at (i+0.5, j+0.5, k)
    x = max(SCALAR(0.5), min(SCALAR(N - 0.5), x));
    y = max(SCALAR(0.5), min(SCALAR(M - 0.5), y));
    z = max(SCALAR(0), min(SCALAR(L), z));
    
    int i = int(x - 0.5);
    int j = int(y - 0.5);
    int k = int(z);
    
    SCALAR fx = x - 0.5 - i;
    SCALAR fy = y - 0.5 - j;
    SCALAR fz = z - k;
    
    i = max(0, min(N - 2, i));
    j = max(0, min(M - 2, j));
    k = max(0, min(L - 1, k));
    
    SCALAR c000 = w_vel[IXYZ(i,     j,     k,     N, M)];
    SCALAR c001 = w_vel[IXYZ(i,     j,     k + 1, N, M)];
    SCALAR c010 = w_vel[IXYZ(i,     j + 1, k,     N, M)];
    SCALAR c011 = w_vel[IXYZ(i,     j + 1, k + 1, N, M)];
    SCALAR c100 = w_vel[IXYZ(i + 1, j,     k,     N, M)];
    SCALAR c101 = w_vel[IXYZ(i + 1, j,     k + 1, N, M)];
    SCALAR c110 = w_vel[IXYZ(i + 1, j + 1, k,     N, M)];
    SCALAR c111 = w_vel[IXYZ(i + 1, j + 1, k + 1, N, M)];
    
    SCALAR c00 = lerp(c000, c001, fz);
    SCALAR c01 = lerp(c010, c011, fz);
    SCALAR c10 = lerp(c100, c101, fz);
    SCALAR c11 = lerp(c110, c111, fz);
    
    SCALAR c0 = lerp(c00, c01, fy);
    SCALAR c1 = lerp(c10, c11, fy);
    
    return lerp(c0, c1, fx);
}

/**
 * @brief Backtrace for 3D advection
 */
template <typename SCALAR>
void backtrace_3d(const int N, const int M, const int L, SCALAR x, SCALAR y, SCALAR z, SCALAR dt,
                  const std::vector<SCALAR> &u_vel, const std::vector<SCALAR> &v_vel, 
                  const std::vector<SCALAR> &w_vel, SCALAR &back_x, SCALAR &back_y, SCALAR &back_z)
{
    // RK2 integration for better accuracy
    SCALAR u_mid = interpolate_u(N, M, L, u_vel, x, y, z);
    SCALAR v_mid = interpolate_v(N, M, L, v_vel, x, y, z);
    SCALAR w_mid = interpolate_w(N, M, L, w_vel, x, y, z);
    
    SCALAR x_mid = x - 0.5f * dt * u_mid;
    SCALAR y_mid = y - 0.5f * dt * v_mid;
    SCALAR z_mid = z - 0.5f * dt * w_mid;
    
    SCALAR u_final = interpolate_u(N, M, L, u_vel, x_mid, y_mid, z_mid);
    SCALAR v_final = interpolate_v(N, M, L, v_vel, x_mid, y_mid, z_mid);
    SCALAR w_final = interpolate_w(N, M, L, w_vel, x_mid, y_mid, z_mid);
    
    back_x = x - dt * u_final;
    back_y = y - dt * v_final;
    back_z = z - dt * w_final;
}

/**
 * @brief Advection of scalar field using MAC grid velocity in 3D
 */
template <typename SCALAR>
void advection_dye_3d(const int N, const int M, const int L,
                      const std::vector<SCALAR> &u_vel, const std::vector<SCALAR> &v_vel, 
                      const std::vector<SCALAR> &w_vel, const std::vector<SCALAR> &old_field, 
                      std::vector<SCALAR> &new_field, SCALAR dt)
{
    for (int k = 0; k < L; ++k) {
        for (int j = 0; j < M; ++j) {
            for (int i = 0; i < N; ++i) {
                int idx = IXYZ(i, j, k, N, M);
                
                // Current cell center position
                SCALAR x = static_cast<SCALAR>(i) + 0.5f;
                SCALAR y = static_cast<SCALAR>(j) + 0.5f;
                SCALAR z = static_cast<SCALAR>(k) + 0.5f;
                
                // Backtrace to find source position
                SCALAR back_x, back_y, back_z;
                backtrace_3d(N, M, L, x, y, z, dt, u_vel, v_vel, w_vel, back_x, back_y, back_z);
                
                // Interpolate at source position
                new_field[idx] = trilerp_scalar(N, M, L, old_field, back_x, back_y, back_z);
            }
        }
    }
}

/**
 * @brief Velocity advection for MAC grid in 3D
 */
template <typename SCALAR>
void advection_velocity_3d(const int N, const int M, const int L,
                           const std::vector<SCALAR> &u_vel, const std::vector<SCALAR> &v_vel, 
                           const std::vector<SCALAR> &w_vel, std::vector<SCALAR> &new_u_vel, 
                           std::vector<SCALAR> &new_v_vel, std::vector<SCALAR> &new_w_vel, SCALAR dt)
{
    // Advect u component
    for (int k = 0; k < L; ++k) {
        for (int j = 0; j < M; ++j) {
            for (int i = 0; i <= N; ++i) {
                int idx = IXYZ(i, j, k, N + 1, M);
                
                // u component position
                SCALAR x = static_cast<SCALAR>(i);
                SCALAR y = static_cast<SCALAR>(j) + 0.5f;
                SCALAR z = static_cast<SCALAR>(k) + 0.5f;
                
                SCALAR back_x, back_y, back_z;
                backtrace_3d(N, M, L, x, y, z, dt, u_vel, v_vel, w_vel, back_x, back_y, back_z);
                
                new_u_vel[idx] = interpolate_u(N, M, L, u_vel, back_x, back_y, back_z);
            }
        }
    }
    
    // Advect v component
    for (int k = 0; k < L; ++k) {
        for (int j = 0; j <= M; ++j) {
            for (int i = 0; i < N; ++i) {
                int idx = IXYZ(i, j, k, N, M + 1);
                
                // v component position
                SCALAR x = static_cast<SCALAR>(i) + 0.5f;
                SCALAR y = static_cast<SCALAR>(j);
                SCALAR z = static_cast<SCALAR>(k) + 0.5f;
                
                SCALAR back_x, back_y, back_z;
                backtrace_3d(N, M, L, x, y, z, dt, u_vel, v_vel, w_vel, back_x, back_y, back_z);
                
                new_v_vel[idx] = interpolate_v(N, M, L, v_vel, back_x, back_y, back_z);
            }
        }
    }
    
    // Advect w component
    for (int k = 0; k <= L; ++k) {
        for (int j = 0; j < M; ++j) {
            for (int i = 0; i < N; ++i) {
                int idx = IXYZ(i, j, k, N, M);
                
                // w component position
                SCALAR x = static_cast<SCALAR>(i) + 0.5f;
                SCALAR y = static_cast<SCALAR>(j) + 0.5f;
                SCALAR z = static_cast<SCALAR>(k);
                
                SCALAR back_x, back_y, back_z;
                backtrace_3d(N, M, L, x, y, z, dt, u_vel, v_vel, w_vel, back_x, back_y, back_z);
                
                new_w_vel[idx] = interpolate_w(N, M, L, w_vel, back_x, back_y, back_z);
            }
        }
    }
}

/**
 * @brief Diffusion of velocity field for MAC grid in 3D
 * @param N width
 * @param M height  
 * @param L depth
 * @param u_vel current u component field (size: (N+1) x M x L)
 * @param v_vel current v component field (size: N x (M+1) x L)
 * @param w_vel current w component field (size: N x M x (L+1))
 * @param new_u_vel new u component field after diffusion
 * @param new_v_vel new v component field after diffusion
 * @param new_w_vel new w component field after diffusion
 * @param nu kinematic viscosity
 * @param dt time step
 */
template <typename SCALAR>
void diffusion_velocity_3d(const int N, const int M, const int L,
                           const std::vector<SCALAR> &u_vel, const std::vector<SCALAR> &v_vel, 
                           const std::vector<SCALAR> &w_vel,
                           std::vector<SCALAR> &new_u_vel, std::vector<SCALAR> &new_v_vel, 
                           std::vector<SCALAR> &new_w_vel,
                           SCALAR nu, SCALAR dt)
{
    // u^{t+1} = u^t + nu*dt*laplace(u^t)
    SCALAR dx = SCALAR(1.0);
    SCALAR coeff = nu * dt / (dx * dx);

    // Initialize output arrays
    new_u_vel = u_vel;
    new_v_vel = v_vel;
    new_w_vel = w_vel;

    // Diffusion of u component
    for (int k = 0; k < L; ++k) {
        for (int j = 0; j < M; ++j) {
            for (int i = 1; i < N; ++i) {  // Skip boundaries at i=0 and i=N
                int id = IXYZ(i, j, k, N + 1, M);
                SCALAR u_center = u_vel[id];
                
                // 6-point stencil for 3D Laplacian
                SCALAR u_left   = u_vel[IXYZ(i - 1, j, k, N + 1, M)];
                SCALAR u_right  = u_vel[IXYZ(i + 1, j, k, N + 1, M)];
                SCALAR u_down   = (j > 0)     ? u_vel[IXYZ(i, j - 1, k, N + 1, M)] : u_center;
                SCALAR u_up     = (j < M - 1) ? u_vel[IXYZ(i, j + 1, k, N + 1, M)] : u_center;
                SCALAR u_back   = (k > 0)     ? u_vel[IXYZ(i, j, k - 1, N + 1, M)] : u_center;
                SCALAR u_front  = (k < L - 1) ? u_vel[IXYZ(i, j, k + 1, N + 1, M)] : u_center;
                
                SCALAR lap = u_left + u_right + u_down + u_up + u_back + u_front - SCALAR(6.0) * u_center;
                new_u_vel[id] = u_center + coeff * lap;
            }
        }
    }

    // Diffusion of v component  
    for (int k = 0; k < L; ++k) {
        for (int j = 1; j < M; ++j) {  // Skip boundaries at j=0 and j=M
            for (int i = 0; i < N; ++i) {
                int id = IXYZ(i, j, k, N, M + 1);
                SCALAR v_center = v_vel[id];
                
                // 6-point stencil for 3D Laplacian
                SCALAR v_left   = (i > 0)     ? v_vel[IXYZ(i - 1, j, k, N, M + 1)] : v_center;
                SCALAR v_right  = (i < N - 1) ? v_vel[IXYZ(i + 1, j, k, N, M + 1)] : v_center;
                SCALAR v_down   = v_vel[IXYZ(i, j - 1, k, N, M + 1)];
                SCALAR v_up     = v_vel[IXYZ(i, j + 1, k, N, M + 1)];
                SCALAR v_back   = (k > 0)     ? v_vel[IXYZ(i, j, k - 1, N, M + 1)] : v_center;
                SCALAR v_front  = (k < L - 1) ? v_vel[IXYZ(i, j, k + 1, N, M + 1)] : v_center;
                
                SCALAR lap = v_left + v_right + v_down + v_up + v_back + v_front - SCALAR(6.0) * v_center;
                new_v_vel[id] = v_center + coeff * lap;
            }
        }
    }

    // Diffusion of w component
    for (int k = 1; k < L; ++k) {  // Skip boundaries at k=0 and k=L
        for (int j = 0; j < M; ++j) {
            for (int i = 0; i < N; ++i) {
                int id = IXYZ(i, j, k, N, M);
                SCALAR w_center = w_vel[id];
                
                // 6-point stencil for 3D Laplacian
                SCALAR w_left   = (i > 0)     ? w_vel[IXYZ(i - 1, j, k, N, M)] : w_center;
                SCALAR w_right  = (i < N - 1) ? w_vel[IXYZ(i + 1, j, k, N, M)] : w_center;
                SCALAR w_down   = (j > 0)     ? w_vel[IXYZ(i, j - 1, k, N, M)] : w_center;
                SCALAR w_up     = (j < M - 1) ? w_vel[IXYZ(i, j + 1, k, N, M)] : w_center;
                SCALAR w_back   = w_vel[IXYZ(i, j, k - 1, N, M)];
                SCALAR w_front  = w_vel[IXYZ(i, j, k + 1, N, M)];
                
                SCALAR lap = w_left + w_right + w_down + w_up + w_back + w_front - SCALAR(6.0) * w_center;
                new_w_vel[id] = w_center + coeff * lap;
            }
        }
    }
}

/**
 * @brief Apply global forces in 3D
 */
template <typename SCALAR>
void apply_force_3d(const int N, const int M, const int L, 
                    std::vector<SCALAR> &u_vel, std::vector<SCALAR> &v_vel, std::vector<SCALAR> &w_vel,
                    const SCALAR g = -9.8f, const SCALAR dt = 0.03f)
{
    // Apply gravity to v component (downward)
    int v_total = N * (M + 1) * L;
    for(int i = 0; i < v_total; i++)
    {
        v_vel[i] += g * 1.0 * dt;
    }
}

/**
 * @brief Get divergence of velocity field for MAC grid in 3D
 */
template <typename SCALAR>
void get_divergence_3d(const int N, const int M, const int L, 
                       const std::vector<SCALAR> &u_vel, const std::vector<SCALAR> &v_vel, 
                       const std::vector<SCALAR> &w_vel, std::vector<SCALAR> &divergence)
{
    for (int k = 0; k < L; ++k) {
        for (int j = 0; j < M; ++j) {
            for (int i = 0; i < N; ++i) {
                int idx = IXYZ(i, j, k, N, M);
                
                SCALAR du_dx = u_vel[IXYZ(i + 1, j, k, N + 1, M)] - u_vel[IXYZ(i, j, k, N + 1, M)];
                SCALAR dv_dy = v_vel[IXYZ(i, j + 1, k, N, M + 1)] - v_vel[IXYZ(i, j, k, N, M + 1)];
                SCALAR dw_dz = w_vel[IXYZ(i, j, k + 1, N, M)] - w_vel[IXYZ(i, j, k, N, M)];
                
                divergence[idx] = du_dx + dv_dy + dw_dz;
            }
        }
    }
}

/**
 * @brief Single iteration of Gauss-Seidel method for pressure in 3D
 */
template <typename SCALAR>
void pressure_gauss_seidel_3d(const int N, const int M, const int L, 
                              const std::vector<SCALAR> &divergence,
                              const std::vector<SCALAR> &pressure, std::vector<SCALAR> &new_pressure)
{
    for (int k = 0; k < L; ++k) {
        for (int j = 0; j < M; ++j) {
            for (int i = 0; i < N; ++i) {
                int idx = IXYZ(i, j, k, N, M);
                
                SCALAR sum = 0.0f;
                int count = 0;
                
                // Neighboring cells
                if (i > 0) {
                    sum += new_pressure[IXYZ(i - 1, j, k, N, M)];
                    count++;
                }
                if (i < N - 1) {
                    sum += new_pressure[IXYZ(i + 1, j, k, N, M)];
                    count++;
                }
                if (j > 0) {
                    sum += new_pressure[IXYZ(i, j - 1, k, N, M)];
                    count++;
                }
                if (j < M - 1) {
                    sum += new_pressure[IXYZ(i, j + 1, k, N, M)];
                    count++;
                }
                if (k > 0) {
                    sum += new_pressure[IXYZ(i, j, k - 1, N, M)];
                    count++;
                }
                if (k < L - 1) {
                    sum += new_pressure[IXYZ(i, j, k + 1, N, M)];
                    count++;
                }
                
                new_pressure[idx] = (sum - divergence[idx]) / SCALAR(count);
            }
        }
    }
}

/**
 * @brief Apply pressure gradient to velocity field for MAC grid in 3D
 */
template <typename SCALAR>
void subtract_gradient_3d(const int N, const int M, const int L, 
                          std::vector<SCALAR> &u_vel, std::vector<SCALAR> &v_vel, 
                          std::vector<SCALAR> &w_vel, const std::vector<SCALAR> &pressure)
{
    // Update u component
    for (int k = 0; k < L; ++k) {
        for (int j = 0; j < M; ++j) {
            for (int i = 1; i < N; ++i) {
                int idx = IXYZ(i, j, k, N + 1, M);
                u_vel[idx] -= (pressure[IXYZ(i, j, k, N, M)] - pressure[IXYZ(i - 1, j, k, N, M)]);
            }
        }
    }
    
    // Update v component
    for (int k = 0; k < L; ++k) {
        for (int j = 1; j < M; ++j) {
            for (int i = 0; i < N; ++i) {
                int idx = IXYZ(i, j, k, N, M + 1);
                v_vel[idx] -= (pressure[IXYZ(i, j, k, N, M)] - pressure[IXYZ(i, j - 1, k, N, M)]);
            }
        }
    }
    
    // Update w component
    for (int k = 1; k < L; ++k) {
        for (int j = 0; j < M; ++j) {
            for (int i = 0; i < N; ++i) {
                int idx = IXYZ(i, j, k, N, M);
                w_vel[idx] -= (pressure[IXYZ(i, j, k, N, M)] - pressure[IXYZ(i, j, k - 1, N, M)]);
            }
        }
    }
}

/**
 * @brief Apply dissipation to scalar field in 3D
 */
template <typename SCALAR>
void dissipate_3d(const int N, const int M, const int L, std::vector<SCALAR> &field, 
                  const SCALAR dissipation_rate, const SCALAR dt)
{
    SCALAR factor = SCALAR(1) / (SCALAR(1) + dt * dissipation_rate);
    for (int k = 0; k < L; ++k) {
        for (int j = 0; j < M; ++j) {
            for (int i = 0; i < N; ++i) {
                int idx = IXYZ(i, j, k, N, M);
                field[idx] *= factor;
            }
        }
    }
}

template <typename SCALAR> SCALAR smooth_step(SCALAR a, SCALAR x)
{
    if (x < 0) return 0;
    if (x > a) return 1;
    SCALAR t = x / a;
    return t * t * (3 - 2 * t);
}

#endif // STABLE_FLUID_3D_H
