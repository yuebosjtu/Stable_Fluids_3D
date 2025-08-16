#include "fluid_simulator.h"
#include "../visualize/colorramp.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <direct.h>  // For _mkdir on Windows

FluidSimulator3D::FluidSimulator3D(const SimulationParams& params)
    : params_(params), current_time_(0.0f), save_frame_(0), current_frame_(0), 
      poisson_matrix_(nullptr), amg_initialized_(false)
{
    InitializeFields();
}

FluidSimulator3D::~FluidSimulator3D()
{
    if (velocity_file_.is_open())
        velocity_file_.close();
    if (pressure_file_.is_open())
        pressure_file_.close();
    if (dye_file_.is_open())
        dye_file_.close();
    
    delete u_pair_;
    delete v_pair_;
    delete w_pair_;
    delete pressure_pair_;
    delete dye_pair_;
    
    delete poisson_matrix_;
    for (auto* matrix : A_L_)
        delete matrix;
    for (auto* matrix : R_L_)
        delete matrix;
    for (auto* matrix : P_L_)
        delete matrix;
}

void FluidSimulator3D::InitializeFields(const std::vector<float>* initial_u_field,
                                       const std::vector<float>* initial_v_field,
                                       const std::vector<float>* initial_w_field,
                                       const std::vector<float>* initial_dye_field)
{
    int u_total_cells = (params_.width + 1) * params_.height * params_.depth;
    int v_total_cells = params_.width * (params_.height + 1) * params_.depth;
    int w_total_cells = params_.width * params_.height * (params_.depth + 1);
    int p_total_cells = params_.width * params_.height * params_.depth;
    
    // Initialize u velocity component field
    u_field_.resize(u_total_cells, 0.0f);
    u_temp_.resize(u_total_cells, 0.0f);
    
    // Set initial u velocity field if provided
    if (initial_u_field != nullptr)
    {
        if (initial_u_field->size() == u_total_cells)
        {
            u_field_ = *initial_u_field;
        }
        else
        {
            std::cerr << "Warning: Initial u field size mismatch. Expected " 
                      << u_total_cells << ", got " << initial_u_field->size() << std::endl;
        }
    }
    
    u_pair_ = new TexPair<std::vector<float> >(u_field_, u_temp_);
    
    // Initialize v velocity component field
    v_field_.resize(v_total_cells, 0.0f);
    v_temp_.resize(v_total_cells, 0.0f);
    
    // Set initial v velocity field if provided
    if (initial_v_field != nullptr)
    {
        if (initial_v_field->size() == v_total_cells)
        {
            v_field_ = *initial_v_field;
        }
        else
        {
            std::cerr << "Warning: Initial v field size mismatch. Expected " 
                      << v_total_cells << ", got " << initial_v_field->size() << std::endl;
        }
    }
    
    v_pair_ = new TexPair<std::vector<float> >(v_field_, v_temp_);
    
    // Initialize w velocity component field
    w_field_.resize(w_total_cells, 0.0f);
    w_temp_.resize(w_total_cells, 0.0f);
    
    // Set initial w velocity field if provided
    if (initial_w_field != nullptr)
    {
        if (initial_w_field->size() == w_total_cells)
        {
            w_field_ = *initial_w_field;
        }
        else
        {
            std::cerr << "Warning: Initial w field size mismatch. Expected " 
                      << w_total_cells << ", got " << initial_w_field->size() << std::endl;
        }
    }
    
    w_pair_ = new TexPair<std::vector<float> >(w_field_, w_temp_);
    
    // Initialize pressure fields
    pressure_field_.resize(p_total_cells, 0.0f);
    pressure_temp_.resize(p_total_cells, 0.0f);
    pressure_pair_ = new TexPair<std::vector<float> >(pressure_field_, pressure_temp_);
    
    // Initialize dye fields
    dye_field_.resize(p_total_cells, 0.0f);
    dye_temp_.resize(p_total_cells, 0.0f);
    
    // Set initial dye field if provided
    if (initial_dye_field != nullptr)
    {
        if (initial_dye_field->size() == p_total_cells)
        {
            dye_field_ = *initial_dye_field;
        }
        else
        {
            std::cerr << "Warning: Initial dye field size mismatch. Expected " 
                      << p_total_cells << ", got " << initial_dye_field->size() << std::endl;
        }
    }
    
    dye_pair_ = new TexPair<std::vector<float> >(dye_field_, dye_temp_);
    
    divergence_field_.resize(p_total_cells, 0.0f);
}

void FluidSimulator3D::Initialize(const std::string& output_dir,
                                 const std::string& image_dir,
                                const std::vector<float>* initial_u_field,
                                const std::vector<float>* initial_v_field,
                                const std::vector<float>* initial_w_field,
                                const std::vector<float>* initial_dye_field)
{
    InitializeFields(initial_u_field, initial_v_field, initial_w_field, initial_dye_field);
    
    // Create output directory
    _mkdir(output_dir.c_str());
    _mkdir(image_dir.c_str());
    std::cout << "Output directory: " << output_dir << std::endl;
    std::cout << "Image directory: " << image_dir << std::endl;
    
    // Open output files
    velocity_file_.open((output_dir + "/velocity_data.txt").c_str());
    pressure_file_.open((output_dir + "/pressure_data.txt").c_str());
    dye_file_.open((output_dir + "/dye_data.txt").c_str());
    
    // Write headers
    velocity_file_ << "# Velocity field data: frame time x y z vel_x vel_y vel_z" << std::endl;
    pressure_file_ << "# Pressure field data: frame time x y z pressure" << std::endl;
    dye_file_ << "# Dye field data: frame time x y z dye_concentration" << std::endl;
    
    std::cout << "Grid size: " << params_.width << "x" << params_.height << "x" << params_.depth << std::endl;
    std::cout << "Time step: " << params_.dt << "s" << std::endl;
    std::cout << "Total time: " << params_.total_time << "s" << std::endl;
    std::cout << "Output directory: " << output_dir << std::endl;
}

void FluidSimulator3D::Step(const std::string& image_dir)
{
    // 1. Apply external forces (gravity, user input, etc.)
    ApplyForces();
    
    // 2. Advect velocity field
    AdvectVelocity(params_.advect_method_);
    std::cout << "Velocity field advected." << std::endl;

    // 3. Diffuse velocity field (viscosity)
    DiffuseVelocity();
    std::cout << "Velocity field diffused." << std::endl;

    // 4. Project - Solve pressure to ensure incompressibility
    SolvePressure();
    std::cout << "Pressure field solved." << std::endl;
    
    // Advect and diffuse dye field
    AdvectDye();
    
    // Apply dissipation to dye field
    DissipateDye();
    
    // 5. Apply boundary conditions and inside source constraints
    ApplyBoundaryConditions();
    ApplyInsideSource();
    
    // Update time and frame
    current_time_ += params_.dt;
    current_frame_++;
    
    // Progress output
    if (current_frame_ % 10 == 0)
    {
        SaveFrameData(image_dir);
        std::cout << "Frame " << current_frame_ << ", Time: " 
                  << std::fixed << std::setprecision(3) << current_time_ << "s" << std::endl;
    }
}

void FluidSimulator3D::RunSimulation(const std::string& image_dir)
{
    std::cout << "Starting 3D fluid simulation..." << std::endl;
    
    while (!IsComplete())
    {
        Step(image_dir);
    }
    
    std::cout << "Simulation completed!" << std::endl;
    std::cout << "Total frames: " << current_frame_ << std::endl;
    std::cout << "Final time: " << current_time_ << "s" << std::endl;
}

void FluidSimulator3D::ApplyForces()
{
    if (params_.enable_gravity)
    {
        apply_force_3d<float>(params_.width, params_.height, params_.depth,
                             u_pair_->cur, v_pair_->cur, w_pair_->cur,
                             -9.8f, params_.dt);
    }
}

void FluidSimulator3D::InitializeAMGSolver()
{
    if (amg_initialized_)
        return;
    
    // Create Poisson matrix for pressure solve
    int n = params_.width * params_.height * params_.depth;
    poisson_matrix_ = new FixedSparseMatrix3D<float>(n, n);
    setupPoissonMatrix3D(*poisson_matrix_, params_.width, params_.height, params_.depth);
    
    // Initialize AMG hierarchy (simplified - just allocate placeholders)
    A_L_.clear();
    R_L_.clear();
    P_L_.clear();
    S_L_.clear();
    
    // For now, just use the original matrix as the only level
    A_L_.push_back(poisson_matrix_);
    S_L_.push_back(Vec3i(params_.width, params_.height, params_.depth));
    
    amg_initialized_ = true;
}

void FluidSimulator3D::SolvePressure()
{
    // Initialize AMG solver if not already done
    if (!amg_initialized_) {
        InitializeAMGSolver();
    }
    
    // Calculate divergence of velocity field
    get_divergence_3d<float>(params_.width, params_.height, params_.depth, 
                            u_pair_->cur, v_pair_->cur, w_pair_->cur, divergence_field_);
    
    // Clear pressure field
    std::fill(pressure_pair_->cur.begin(), pressure_pair_->cur.end(), 0.0f);
    
    // Solve pressure using AMG PCG solver
    float residual_out;
    int iterations_out;
    float tolerance = 1e-6f;
    int max_iterations = params_.pressure_iterations;
    
    // Copy divergence to RHS
    std::vector<float> rhs = divergence_field_;
    for (auto& val : rhs) {
        val = -val;  // Negative because we want div(u) = 0
    }
    
    bool converged = AMGPCGSolvePrebuilt3D<float>(
        *poisson_matrix_,
        rhs,
        pressure_pair_->cur,
        A_L_,
        R_L_,
        P_L_,
        S_L_,
        1,
        tolerance,
        max_iterations,
        residual_out,
        iterations_out,
        params_.width,
        params_.height,
        params_.depth,
        false
    );
    
    if (!converged) {
        std::cerr << "Warning: Pressure solver did not converge. Residual: " 
                  << residual_out << " after " << iterations_out << " iterations." << std::endl;
    }
    
    // Subtract pressure gradient from velocity
    subtract_gradient_3d<float>(params_.width, params_.height, params_.depth, 
                               u_pair_->cur, v_pair_->cur, w_pair_->cur, pressure_pair_->cur);
}

void FluidSimulator3D::AdvectVelocity(AdvectMethod method)
{
    switch (method) {
        case AdvectMethod::SemiLagrange:
            advection_velocity_3d<float>(params_.width, params_.height, params_.depth,
                                        u_pair_->cur, v_pair_->cur, w_pair_->cur,
                                        u_pair_->nxt, v_pair_->nxt, w_pair_->nxt,
                                        params_.dt);
            break;
        case AdvectMethod::MacCormack:
            macCormackVelocity_3d<float>(params_.width, params_.height, params_.depth,
                                        u_pair_->cur, v_pair_->cur, w_pair_->cur,
                                        u_pair_->nxt, v_pair_->nxt, w_pair_->nxt,
                                        params_.dt);
            break;
        default:
            std::cerr << "Unknown advection method!" << std::endl;
            break;
    }
    u_pair_->Swap();
    v_pair_->Swap();
    w_pair_->Swap();
}

void FluidSimulator3D::DiffuseVelocity()
{
    if (params_.viscosity <= 0.0f) return;
    
    diffusion_velocity_3d<float>(params_.width, params_.height, params_.depth,
                                 u_pair_->cur, v_pair_->cur, w_pair_->cur,
                                 u_pair_->nxt, v_pair_->nxt, w_pair_->nxt,
                                 params_.viscosity, params_.dt);
    u_pair_->Swap();
    v_pair_->Swap();
    w_pair_->Swap();
}

void FluidSimulator3D::AdvectDye()
{
    advection_dye_3d<float>(params_.width, params_.height, params_.depth,
                           u_pair_->cur, v_pair_->cur, w_pair_->cur,
                           dye_pair_->cur, dye_pair_->nxt, params_.dt);
    dye_pair_->Swap();
}

void FluidSimulator3D::DissipateDye()
{
    if (params_.dissipation > 0.0f) {
        dissipate_3d<float>(params_.width, params_.height, params_.depth, 
                           dye_pair_->cur, params_.dissipation, params_.dt);
    }
}

void FluidSimulator3D::SaveFrameData(const std::string& image_dir)
{
    if (!velocity_file_.is_open() || !pressure_file_.is_open() || !dye_file_.is_open())
        return;

    const std::vector<Vector3f>& velocity_field = GetVelocityField();
    const std::vector<float>& pressure_field = pressure_pair_->cur;
    const std::vector<float>& dye_field = dye_pair_->cur;
    
    // TODO: save data to file
    
    SaveVelocityFieldImage(image_dir, velocity_field, current_frame_);
    
    save_frame_++;
}

void FluidSimulator3D::ApplyBoundaryConditions()
{
    // No-slip boundary conditions for u component
    for (int k = 0; k < params_.depth; ++k) {
        for (int j = 0; j < params_.height; ++j) {
            // Left and right boundaries
            u_field_[IXYZ(0, j, k, params_.width + 1, params_.height)] = 0.0f;
            u_field_[IXYZ(params_.width, j, k, params_.width + 1, params_.height)] = 0.0f;
        }
    }
    
    for (int k = 0; k < params_.depth; ++k) {
        for (int i = 0; i <= params_.width; ++i) {
            // Bottom and top boundaries for u component
            if (i > 0 && i < params_.width) {
                float u_center = 0.5f * (u_field_[IXYZ(i, 0, k, params_.width + 1, params_.height)] + 
                                        u_field_[IXYZ(i, params_.height - 1, k, params_.width + 1, params_.height)]);
                u_field_[IXYZ(i, 0, k, params_.width + 1, params_.height)] = -u_center;
                u_field_[IXYZ(i, params_.height - 1, k, params_.width + 1, params_.height)] = -u_center;
            }
        }
    }
    
    // Similar boundary conditions for v and w components...
    // No-slip boundary conditions for v component
    for (int k = 0; k < params_.depth; ++k) {
        for (int i = 0; i < params_.width; ++i) {
            // Bottom and top boundaries
            v_field_[IXYZ(i, 0, k, params_.width, params_.height + 1)] = 0.0f;
            v_field_[IXYZ(i, params_.height, k, params_.width, params_.height + 1)] = 0.0f;
        }
    }
    
    // No-slip boundary conditions for w component
    for (int j = 0; j < params_.height; ++j) {
        for (int i = 0; i < params_.width; ++i) {
            // Back and front boundaries
            w_field_[IXYZ(i, j, 0, params_.width, params_.height)] = 0.0f;
            w_field_[IXYZ(i, j, params_.depth, params_.width, params_.height)] = 0.0f;
        }
    }
}

void FluidSimulator3D::ApplyInsideSource()
{
    // Apply spherical source conditions
    float cx = params_.source_center_x;
    float cy = params_.source_center_y;
    float cz = params_.source_center_z;
    float radius = params_.source_radius;
    
    // Set u velocity inside sphere
    for (int k = 0; k < params_.depth; ++k) {
        for (int j = 0; j < params_.height; ++j) {
            for (int i = 0; i <= params_.width; ++i) {
                float x = static_cast<float>(i);
                float y = static_cast<float>(j) + 0.5f;
                float z = static_cast<float>(k) + 0.5f;
                
                float dx = x - cx;
                float dy = y - cy;
                float dz = z - cz;
                float distance = std::sqrt(dx*dx + dy*dy + dz*dz);
                
                if (distance <= radius) {
                    u_field_[IXYZ(i, j, k, params_.width + 1, params_.height)] = params_.source_velocity_x;
                }
            }
        }
    }
    
    // Set v velocity inside sphere
    for (int k = 0; k < params_.depth; ++k) {
        for (int j = 0; j <= params_.height; ++j) {
            for (int i = 0; i < params_.width; ++i) {
                float x = static_cast<float>(i) + 0.5f;
                float y = static_cast<float>(j);
                float z = static_cast<float>(k) + 0.5f;
                
                float dx = x - cx;
                float dy = y - cy;
                float dz = z - cz;
                float distance = std::sqrt(dx*dx + dy*dy + dz*dz);
                
                if (distance <= radius) {
                    v_field_[IXYZ(i, j, k, params_.width, params_.height + 1)] = params_.source_velocity_y;
                }
            }
        }
    }
    
    // Set w velocity inside sphere
    for (int k = 0; k <= params_.depth; ++k) {
        for (int j = 0; j < params_.height; ++j) {
            for (int i = 0; i < params_.width; ++i) {
                float x = static_cast<float>(i) + 0.5f;
                float y = static_cast<float>(j) + 0.5f;
                float z = static_cast<float>(k);
                
                float dx = x - cx;
                float dy = y - cy;
                float dz = z - cz;
                float distance = std::sqrt(dx*dx + dy*dy + dz*dz);
                
                if (distance <= radius) {
                    w_field_[IXYZ(i, j, k, params_.width, params_.height)] = params_.source_velocity_z;
                }
            }
        }
    }
}


void FluidSimulator3D::SaveVelocityFieldImage(const std::string& image_dir, const std::vector<Vector3f>& velocity_field, int current_frame)
{
    using namespace Mfree;
    // 可视化xy、xz、yz三个截面
    // 1. xy平面（z = depth / 2）
    int z_xy = params_.depth / 2;
    std::vector<float> rgb_xy(params_.width * params_.height * 3);
    ColorRamp color_ramp;
    float max_v_xy = 0.0f;
    for (int j = 0; j < params_.height; ++j) {
        for (int i = 0; i < params_.width; ++i) {
            int idx = IXYZ(i, j, z_xy, params_.width, params_.height);
            float vx = velocity_field[idx](0);
            float vy = velocity_field[idx](1);
            float mag = std::sqrt(vx * vx + vy * vy);
            max_v_xy = std::max(max_v_xy, mag);
        }
    }
    if (max_v_xy < 1e-6f) max_v_xy = 1.0f;
    for (int j = 0; j < params_.height; ++j) {
        for (int i = 0; i < params_.width; ++i) {
            int idx = IXYZ(i, j, z_xy, params_.width, params_.height);
            int rgb_id = ((params_.height - 1 - j) * params_.width + i) * 3;
            float vx = velocity_field[idx](0);
            float vy = velocity_field[idx](1);
            float mag = std::sqrt(vx * vx + vy * vy) / max_v_xy;
            vec3 color(0, 0, 0);
            color_ramp.set_GLcolor(mag, COLOR_JET, color, false);
            rgb_xy[rgb_id] = color.x;
            rgb_xy[rgb_id + 1] = color.y;
            rgb_xy[rgb_id + 2] = color.z;
        }
    }
    std::string dir_xy = image_dir + "/xy_slice";
    _mkdir(dir_xy.c_str());
    char filename_xy[512];
    sprintf(filename_xy, "%s/velocity_xy_frame_%05d.ppm", dir_xy.c_str(), current_frame);
    std::ofstream ppm_xy(filename_xy, std::ios::out | std::ios::binary);
    if (ppm_xy.is_open()) {
        ppm_xy << "P6\n" << params_.width << " " << params_.height << "\n255\n";
        for (int i = 0; i < params_.width * params_.height * 3; ++i) {
            ppm_xy << (unsigned char)(rgb_xy[i] * 255);
        }
        ppm_xy.close();
    }

    // 2. xz平面（y = height / 2）
    int y_xz = params_.height / 2;
    std::vector<float> rgb_xz(params_.width * params_.depth * 3);
    float max_v_xz = 0.0f;
    for (int k = 0; k < params_.depth; ++k) {
        for (int i = 0; i < params_.width; ++i) {
            int idx = IXYZ(i, y_xz, k, params_.width, params_.height);
            float vx = velocity_field[idx](0);
            float vz = velocity_field[idx](2);
            float mag = std::sqrt(vx * vx + vz * vz);
            max_v_xz = std::max(max_v_xz, mag);
        }
    }
    if (max_v_xz < 1e-6f) max_v_xz = 1.0f;
    for (int k = 0; k < params_.depth; ++k) {
        for (int i = 0; i < params_.width; ++i) {
            int idx = IXYZ(i, y_xz, k, params_.width, params_.height);
            int rgb_id = ((params_.depth - 1 - k) * params_.width + i) * 3;
            float vx = velocity_field[idx](0);
            float vz = velocity_field[idx](2);
            float mag = std::sqrt(vx * vx + vz * vz) / max_v_xz;
            vec3 color(0, 0, 0);
            color_ramp.set_GLcolor(mag, COLOR_JET, color, false);
            rgb_xz[rgb_id] = color.x;
            rgb_xz[rgb_id + 1] = color.y;
            rgb_xz[rgb_id + 2] = color.z;
        }
    }
    std::string dir_xz = image_dir + "/xz_slice";
    _mkdir(dir_xz.c_str());
    char filename_xz[512];
    sprintf(filename_xz, "%s/velocity_xz_frame_%05d.ppm", dir_xz.c_str(), current_frame);
    std::ofstream ppm_xz(filename_xz, std::ios::out | std::ios::binary);
    if (ppm_xz.is_open()) {
        ppm_xz << "P6\n" << params_.width << " " << params_.depth << "\n255\n";
        for (int i = 0; i < params_.width * params_.depth * 3; ++i) {
            ppm_xz << (unsigned char)(rgb_xz[i] * 255);
        }
        ppm_xz.close();
    }

    // 3. yz平面（x = width / 2）
    int x_yz = params_.width / 2;
    std::vector<float> rgb_yz(params_.height * params_.depth * 3);
    float max_v_yz = 0.0f;
    for (int k = 0; k < params_.depth; ++k) {
        for (int j = 0; j < params_.height; ++j) {
            int idx = IXYZ(x_yz, j, k, params_.width, params_.height);
            float vy = velocity_field[idx](1);
            float vz = velocity_field[idx](2);
            float mag = std::sqrt(vy * vy + vz * vz);
            max_v_yz = std::max(max_v_yz, mag);
        }
    }
    if (max_v_yz < 1e-6f) max_v_yz = 1.0f;
    for (int k = 0; k < params_.depth; ++k) {
        for (int j = 0; j < params_.height; ++j) {
            int idx = IXYZ(x_yz, j, k, params_.width, params_.height);
            int rgb_id = ((params_.depth - 1 - k) * params_.height + j) * 3;
            float vy = velocity_field[idx](1);
            float vz = velocity_field[idx](2);
            float mag = std::sqrt(vy * vy + vz * vz) / max_v_yz;
            vec3 color(0, 0, 0);
            color_ramp.set_GLcolor(mag, COLOR_JET, color, false);
            rgb_yz[rgb_id] = color.x;
            rgb_yz[rgb_id + 1] = color.y;
            rgb_yz[rgb_id + 2] = color.z;
        }
    }
    std::string dir_yz = image_dir + "/yz_slice";
    _mkdir(dir_yz.c_str());
    char filename_yz[512];
    sprintf(filename_yz, "%s/velocity_yz_frame_%05d.ppm", dir_yz.c_str(), current_frame);
    std::ofstream ppm_yz(filename_yz, std::ios::out | std::ios::binary);
    if (ppm_yz.is_open()) {
        ppm_yz << "P6\n" << params_.height << " " << params_.depth << "\n255\n";
        for (int i = 0; i < params_.height * params_.depth * 3; ++i) {
            ppm_yz << (unsigned char)(rgb_yz[i] * 255);
        }
        ppm_yz.close();
    }
}