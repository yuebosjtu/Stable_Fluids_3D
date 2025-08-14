#ifndef FLUID_SIMULATOR_3D_H
#define FLUID_SIMULATOR_3D_H

#include "stable_fluid.h"
#include "amg_solver.h"
#include <Eigen/Dense>
#include <string>
#include <fstream>
#include <memory>
#include <vector>

/**
 * @brief 3D Fluid simulation class implementing Jos Stam's stable fluid method
 */
class FluidSimulator3D
{
public:
    enum class AdvectMethod {
            SemiLagrange = 0,
            MacCormack = 1
        };

    typedef Eigen::Vector3f Vector3f;
    typedef Eigen::VectorXf VectorXf;
    
    struct SimulationParams
    {
        int width;                    // Grid width (X direction)
        int height;                   // Grid height (Y direction)
        int depth;                    // Grid depth (Z direction)
        float dt;                     // Time step
        float viscosity;              // Fluid viscosity
        float diffusion;              // Dye diffusion constant
        float dissipation;            // Dye dissipation rate
        int pressure_iterations;      // Gauss-Seidel iterations for pressure
        float total_time;             // Total simulation time
        bool enable_gravity;          // Enable gravity
        
        // Spherical source parameters
        float source_center_x;        // Sphere center x coordinate
        float source_center_y;        // Sphere center y coordinate
        float source_center_z;        // Sphere center z coordinate
        float source_radius;          // Sphere radius
        float source_velocity_y;      // Upward velocity magnitude
        float source_velocity_x;      // X direction velocity
        float source_velocity_z;      // Z direction velocity

        // Advect method to use
        AdvectMethod advect_method_ = AdvectMethod::SemiLagrange;
        
        // Constructor with default values
        SimulationParams() : 
            width(32), height(32), depth(32), dt(0.016f), viscosity(0.0001f), diffusion(0.0001f),
            dissipation(0.001f), pressure_iterations(80), total_time(10.0f), enable_gravity(true),
            source_center_x(16.0f), source_center_y(8.0f), source_center_z(16.0f), 
            source_radius(4.0f), source_velocity_y(5.0f), source_velocity_x(0.0f), source_velocity_z(0.0f),
            advect_method_(AdvectMethod::SemiLagrange) {}
    };
    

private:
    SimulationParams params_;
    float current_time_;
    int save_frame_;
    int current_frame_;
    
    // Velocity components fields (MAC grid)
    std::vector<float> u_field_;     // Size: (width+1) x height x depth
    std::vector<float> u_temp_;
    std::vector<float> v_field_;     // Size: width x (height+1) x depth
    std::vector<float> v_temp_;
    std::vector<float> w_field_;     // Size: width x height x (depth+1)
    std::vector<float> w_temp_;

    // Scalar fields stored at cell centers 
    std::vector<float> pressure_field_;    // Size: width x height x depth
    std::vector<float> pressure_temp_;
    std::vector<float> dye_field_;         // Size: width x height x depth
    std::vector<float> dye_temp_;
    std::vector<float> divergence_field_;  // Size: width x height x depth
    
    // TexPair wrappers for easy swapping
    TexPair<std::vector<float> >* u_pair_;
    TexPair<std::vector<float> >* v_pair_;
    TexPair<std::vector<float> >* w_pair_;
    TexPair<std::vector<float> >* pressure_pair_;
    TexPair<std::vector<float> >* dye_pair_;
    
    // AMG Solver data structures
    FixedSparseMatrix3D<float>* poisson_matrix_;
    std::vector<FixedSparseMatrix3D<float>*> A_L_;
    std::vector<FixedSparseMatrix3D<float>*> R_L_;
    std::vector<FixedSparseMatrix3D<float>*> P_L_;
    std::vector<Vec3i> S_L_;
    bool amg_initialized_;
    
    // Output files for saving data
    std::ofstream velocity_file_;
    std::ofstream pressure_file_;
    std::ofstream dye_file_;

public:
    explicit FluidSimulator3D(const SimulationParams& params = SimulationParams());
    
    ~FluidSimulator3D();
    
    /**
     * @brief Initialize the simulation
     * @param output_dir Directory to save output files
     * @param initial_u_field Optional initial u velocity field (size: (width+1) x height x depth)
     * @param initial_v_field Optional initial v velocity field (size: width x (height+1) x depth)
     * @param initial_w_field Optional initial w velocity field (size: width x height x (depth+1))
     * @param initial_dye_field Optional initial dye field (size: width x height x depth)
     */
    void Initialize(const std::string& output_dir = "output",
                    const std::string& image_dir = "image_file",
                   const std::vector<float>* initial_u_field = nullptr,
                   const std::vector<float>* initial_v_field = nullptr,
                   const std::vector<float>* initial_w_field = nullptr,
                   const std::vector<float>* initial_dye_field = nullptr);
    
    /**
     * @brief Run one simulation step
     */
    void Step(const std::string& img_dir = "image_file");
    
    /**
     * @brief Run the complete simulation
     */
    void RunSimulation(const std::string& img_dir = "image_file");
    
    /**
     * @brief Create initial velocity field with spherical source
     * @param u_field Output u velocity field (will be resized)
     * @param v_field Output v velocity field (will be resized)
     * @param w_field Output w velocity field (will be resized)
     */
    void CreateSphericalSourceField(std::vector<float>& u_field, std::vector<float>& v_field, 
                                   std::vector<float>& w_field) const;
    
    float GetCurrentTime() const { return current_time_; }
    
    int GetCurrentFrame() const { return current_frame_; }
    
    bool IsComplete() const { return current_time_ >= params_.total_time; }
    
    const SimulationParams& GetParams() const { return params_; }
    
    const std::vector<float>& GetVelocityField_u() const { return u_field_; }

    const std::vector<float>& GetVelocityField_v() const { return v_field_; }

    const std::vector<float>& GetVelocityField_w() const { return w_field_; }

    const std::vector<Vector3f>& GetVelocityField() const
    {
        static std::vector<Vector3f> velocity_field;
        velocity_field.resize(params_.width * params_.height * params_.depth);
        
        for (int k = 0; k < params_.depth; ++k)
        {
            for (int j = 0; j < params_.height; ++j)
            {
                for (int i = 0; i < params_.width; ++i)
                {
                    int idx = IXYZ(i, j, k, params_.width, params_.height);
                    
                    // Average u velocities at left and right faces
                    float u_left = u_field_[IXYZ(i, j, k, params_.width + 1, params_.height)];
                    float u_right = u_field_[IXYZ(i + 1, j, k, params_.width + 1, params_.height)];
                    float u_avg = 0.5f * (u_left + u_right);
                    
                    // Average v velocities at bottom and top faces
                    float v_bottom = v_field_[IXYZ(i, j, k, params_.width, params_.height + 1)];
                    float v_top = v_field_[IXYZ(i, j + 1, k, params_.width, params_.height + 1)];
                    float v_avg = 0.5f * (v_bottom + v_top);
                    
                    // Average w velocities at back and front faces
                    float w_back = w_field_[IXYZ(i, j, k, params_.width, params_.height)];
                    float w_front = w_field_[IXYZ(i, j, k + 1, params_.width, params_.height)];
                    float w_avg = 0.5f * (w_back + w_front);
                    
                    velocity_field[idx] = Vector3f(u_avg, v_avg, w_avg);
                }
            }
        }
        
        return velocity_field;
    }

    const std::vector<float>& GetPressureField() const { return pressure_field_; }

    const std::vector<float>& GetDyeField() const { return dye_field_; }

    /**
     * @brief Save current frame data to files
     */
    void SaveFrameData(const std::string& image_dir = "image_file");

private:
    /**
     * @brief Initialize all fields
     * @param initial_u_field Optional initial u velocity field (size: (width+1) x height x depth). If nullptr, zero initialized.
     * @param initial_v_field Optional initial v velocity field (size: width x (height+1) x depth). If nullptr, zero initialized.
     * @param initial_w_field Optional initial w velocity field (size: width x height x (depth+1)). If nullptr, zero initialized.
     * @param initial_dye_field Optional initial dye field (size: width x height x depth). If nullptr, zero initialized.
     */
    void InitializeFields(const std::vector<float>* initial_u_field = nullptr,
                         const std::vector<float>* initial_v_field = nullptr,
                         const std::vector<float>* initial_w_field = nullptr,
                         const std::vector<float>* initial_dye_field = nullptr);
    
    /**
     * @brief Apply external forces (gravity, user input, etc.)
     */
    void ApplyForces();
    
    /**
     * @brief Initialize AMG solver for pressure solve
     */
    void InitializeAMGSolver();
    
    /**
     * @brief Advect velocity field
     */
    void AdvectVelocity(AdvectMethod method);
    
    /**
     * @brief Solve for pressure to make velocity field divergence-free
     */
    void SolvePressure();

    /**
     * @brief Advect dye field
     */
    void AdvectDye();

    void DiffuseVelocity();
    
    /**
     * @brief Dissipate (decay) dye field
     */
    void DissipateDye();
    
    /**
     * @brief Apply boundary conditions (no-slip walls)
     */
    void ApplyBoundaryConditions();
    
    /**
     * @brief Apply inside source conditions (constant velocity in sphere)
     */
    void ApplyInsideSource();
    
    /**
     * @brief Save velocity field as VTK file for 3D visualization
     */
    void SaveVelocityFieldVTK(const std::string& image_dir, int current_frame);
};

#endif // FLUID_SIMULATOR_3D_H
