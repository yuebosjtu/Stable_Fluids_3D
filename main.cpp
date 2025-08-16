#include "fluid_simulator.h"
#include "utils.h"
#include <iostream>
#include <vector>
#include <cmath>

int main()
{
    // Create 3D simulation parameters with spherical source
    FluidSimulator3D::SimulationParams params;
    params.width = 128;
    params.height = 256;
    params.depth = 128;
    params.dt = 0.1f;
    params.total_time = 100.0f;
    params.enable_gravity = false;  // Disable gravity to see pure source effect
    params.pressure_iterations = 600;

    // Configure spherical source in lower part of domain
    params.source_center_x = 64.0f;
    params.source_center_y = 32.0f;
    params.source_center_z = 64.0f;
    params.source_radius = 8.0f;        // 8 grid units radius
    params.source_velocity_x = 0.0f;    // No horizontal velocity
    params.source_velocity_y = 5.0f;    // Upward velocity magnitude
    params.source_velocity_z = 0.0f;    // No Z direction velocity

    params.advect_method_ = FluidSimulator3D::AdvectMethod::SemiLagrange;

    // Create 3D fluid simulator
    FluidSimulator3D simulator(params);
    
    // Create initial velocity field with spherical source
    std::vector<float> initial_u_field, initial_v_field, initial_w_field;
    CreateSphericalSourceField(params.width, params.height, params.depth,
        params.source_center_x, params.source_center_y, params.source_center_z,
        params.source_radius, params.source_velocity_x, params.source_velocity_y, params.source_velocity_z,
        initial_u_field, initial_v_field, initial_w_field);

    // Create initial dye field
    std::vector<float> initial_dye_field;
    CreateSphericalDyeField(params.width, params.height, params.depth,
        params.source_center_x, params.source_center_y, params.source_center_z,
        params.source_radius, 1.0f, initial_dye_field);

    std::string output_dir = "fluid_output_3d";
    std::string img_dir = "velocity_field_3d";

    // Initialize simulation with custom fields
    simulator.Initialize(output_dir, img_dir, &initial_u_field, &initial_v_field, 
                        &initial_w_field, &initial_dye_field);

    // Run simulation
    simulator.RunSimulation(img_dir);

    std::cout << "\n=== 3D Simulation completed! ===" << std::endl;
    std::cout << "The simulation shows:" << std::endl;
    std::cout << "- A spherical source at (" << params.source_center_x << ", " 
              << params.source_center_y << ", " << params.source_center_z << ")" << std::endl;
    std::cout << "- Source radius: " << params.source_radius << " grid units" << std::endl;
    std::cout << "- Source velocity: (" << params.source_velocity_x << ", " 
              << params.source_velocity_y << ", " << params.source_velocity_z << ") units/s" << std::endl;
    std::cout << "- No-slip boundary conditions on all walls" << std::endl;
    std::cout << "- Constant velocity maintained inside the source region" << std::endl;
    std::cout << "- VTK files saved for 3D visualization in ParaView" << std::endl;
    
    return 0;
}
