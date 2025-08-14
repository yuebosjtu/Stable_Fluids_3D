#ifndef UTILS_3D_H
#define UTILS_3D_H

#include <iostream>
#include <vector>
#include <cmath>

/**
 * @brief Create spherical source velocity field for 3D simulation
 * @param width Grid width
 * @param height Grid height
 * @param depth Grid depth
 * @param source_center_x Sphere center x coordinate
 * @param source_center_y Sphere center y coordinate
 * @param source_center_z Sphere center z coordinate
 * @param source_radius Sphere radius
 * @param source_velocity_x X component of velocity inside sphere
 * @param source_velocity_y Y component of velocity inside sphere
 * @param source_velocity_z Z component of velocity inside sphere
 * @param u_field Output u velocity field (will be resized)
 * @param v_field Output v velocity field (will be resized)
 * @param w_field Output w velocity field (will be resized)
 */
void CreateSphericalSourceField(int width, int height, int depth,
                               float source_center_x, float source_center_y, float source_center_z,
                               float source_radius, float source_velocity_x, float source_velocity_y, float source_velocity_z,
                               std::vector<float>& u_field, std::vector<float>& v_field, std::vector<float>& w_field)
{
    int u_total_cells = (width + 1) * height * depth;
    int v_total_cells = width * (height + 1) * depth;
    int w_total_cells = width * height * (depth + 1);

    u_field.resize(u_total_cells, 0.0f);
    v_field.resize(v_total_cells, 0.0f);
    w_field.resize(w_total_cells, 0.0f);
    
    // Set u velocity in spherical region
    for (int k = 0; k < depth; ++k)
    {
        for (int j = 0; j < height; ++j)
        {
            for (int i = 0; i <= width; ++i)
            {
                // u component is stored at (i, j+0.5, k+0.5)
                float x = static_cast<float>(i);
                float y = static_cast<float>(j) + 0.5f;
                float z = static_cast<float>(k) + 0.5f;
                
                float dx = x - source_center_x;
                float dy = y - source_center_y;
                float dz = z - source_center_z;
                float distance = std::sqrt(dx*dx + dy*dy + dz*dz);
                
                if (distance <= source_radius)
                {
                    u_field[(width + 1) * height * k + (width + 1) * j + i] = source_velocity_x;
                }
            }
        }
    }
    
    // Set v velocity in spherical region
    for (int k = 0; k < depth; ++k)
    {
        for (int j = 0; j <= height; ++j)
        {
            for (int i = 0; i < width; ++i)
            {
                // v component is stored at (i+0.5, j, k+0.5)
                float x = static_cast<float>(i) + 0.5f;
                float y = static_cast<float>(j);
                float z = static_cast<float>(k) + 0.5f;
                
                float dx = x - source_center_x;
                float dy = y - source_center_y;
                float dz = z - source_center_z;
                float distance = std::sqrt(dx*dx + dy*dy + dz*dz);
                
                if (distance <= source_radius)
                {
                    v_field[width * (height + 1) * k + width * j + i] = source_velocity_y;
                }
            }
        }
    }
    
    // Set w velocity in spherical region
    for (int k = 0; k <= depth; ++k)
    {
        for (int j = 0; j < height; ++j)
        {
            for (int i = 0; i < width; ++i)
            {
                // w component is stored at (i+0.5, j+0.5, k)
                float x = static_cast<float>(i) + 0.5f;
                float y = static_cast<float>(j) + 0.5f;
                float z = static_cast<float>(k);
                
                float dx = x - source_center_x;
                float dy = y - source_center_y;
                float dz = z - source_center_z;
                float distance = std::sqrt(dx*dx + dy*dy + dz*dz);
                
                if (distance <= source_radius)
                {
                    w_field[width * height * k + width * j + i] = source_velocity_z;
                }
            }
        }
    }

    std::cout << "Created spherical source field at (" << source_center_x << ", " 
              << source_center_y << ", " << source_center_z << ") with radius " << source_radius 
              << " and velocity (" << source_velocity_x << ", " << source_velocity_y 
              << ", " << source_velocity_z << ")" << std::endl;
}

/**
 * @brief Create initial dye field with spherical distribution
 */
void CreateSphericalDyeField(int width, int height, int depth,
                            float center_x, float center_y, float center_z,
                            float radius, float dye_concentration,
                            std::vector<float>& dye_field)
{
    int total_cells = width * height * depth;
    dye_field.resize(total_cells, 0.0f);
    
    for (int k = 0; k < depth; ++k)
    {
        for (int j = 0; j < height; ++j)
        {
            for (int i = 0; i < width; ++i)
            {
                // Dye is stored at cell centers (i+0.5, j+0.5, k+0.5)
                float x = static_cast<float>(i) + 0.5f;
                float y = static_cast<float>(j) + 0.5f;
                float z = static_cast<float>(k) + 0.5f;
                
                float dx = x - center_x;
                float dy = y - center_y;
                float dz = z - center_z;
                float distance = std::sqrt(dx*dx + dy*dy + dz*dz);
                
                if (distance <= radius)
                {
                    // Smooth falloff from center
                    float falloff = 1.0f - (distance / radius);
                    dye_field[width * height * k + width * j + i] = dye_concentration * falloff * falloff;
                }
            }
        }
    }
    
    std::cout << "Created spherical dye field at (" << center_x << ", " 
              << center_y << ", " << center_z << ") with radius " << radius 
              << " and concentration " << dye_concentration << std::endl;
}

#endif // UTILS_3D_H
