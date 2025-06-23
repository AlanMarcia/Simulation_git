#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <algorithm> // For std::min, std::max
#include <omp.h>     // For OpenMP

#if defined(_WIN32)
    #include <direct.h> // For _mkdir
    #define MKDIR(path) _mkdir(path)
    #define STAT_STRUCT _stat
    #define STAT_FUNC _stat
#else
    #include <sys/stat.h> // For mkdir and stat
    #include <sys/types.h> // For mode_t
    #define MKDIR(path) mkdir(path, 0755) // 0755 permissions
    #define STAT_STRUCT stat
    #define STAT_FUNC stat
#endif

// Helper function to save a 2D vector to a CSV file
void saveToCSV(const std::vector<std::vector<double>>& data, const std::string& filename) {
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }
    for (size_t j = 0; j < data[0].size(); ++j) { // Iterate over y first for typical CSV/plot orientation
        for (size_t i = 0; i < data.size(); ++i) { // Then x
            outfile << data[i][j] << (i == data.size() - 1 ? "" : ",");
        }
        outfile << std::endl;
    }
    outfile.close();
    std::cout << "Data saved to " << filename << std::endl;
}

// Helper function to save a 1D vector to a CSV file (single column)
void saveCoordinatesToCSV(const std::vector<double>& coords, const std::string& filename) {
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }
    for (size_t i = 0; i < coords.size(); ++i) {
        outfile << coords[i] << std::endl;
    }
    outfile.close();
    std::cout << "Coordinates saved to " << filename << std::endl;
}

// Helper function to save geometry parameters
void saveGeometryParamsToCSV(const std::string& filename,
                             double h_val,
                             double x_fs, double x_sl,
                             double y_sibh, double y_th, // y_si_base_height, y_teeth_height
                             double y_vgt,
                             double x_tw, double x_ts,   // x_teeth_width, x_teeth_spacing
                             double H_tot) {
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }
    outfile << std::fixed << std::setprecision(10); // Ensure precision
    outfile << "h," << h_val << std::endl;
    outfile << "x_free_space," << x_fs << std::endl;
    outfile << "x_structure_len," << x_sl << std::endl;
    outfile << "y_si_base_height," << y_sibh << std::endl;
    outfile << "y_teeth_height," << y_th << std::endl;
    outfile << "y_vacuum_gap_thick," << y_vgt << std::endl;
    outfile << "x_teeth_width," << x_tw << std::endl;
    outfile << "x_teeth_spacing," << x_ts << std::endl;
    outfile << "H_total," << H_tot << std::endl;
    outfile.close();
    std::cout << "Geometry parameters saved to " << filename << std::endl;
}


int main() {
    // --- Parameters ---
    const double h = 0.5; // Grid spacing in micrometers (µm)    // Dimensions in µm
    const double L_total = 320.0;
    // H_total = 2 * (y_si_base_height + y_teeth_height) + y_vacuum_gap_thick + vacuum_padding
    // H_total = 2 * (5.0 + 5.0) + 10.0 + 20.0 = 2 * 10.0 + 10.0 + 20.0 = 20.0 + 30.0 = 50.0
    const double H_total = 50.0; // Updated: added 20 μm (10 μm top + 10 μm bottom vacuum padding)

    const double x_free_space = 10.0;
    const double x_structure_len = 300.0;
    
    // New geometry parameters for silicon layers
    const double y_si_base_height = 5.0;   // µm
    const double y_teeth_height = 5.0;     // µm
    const double x_teeth_width = 10.0;     // µm
    const double x_teeth_spacing = 10.0;   // µm
    const double tooth_period = x_teeth_width + x_teeth_spacing;

    const double y_vacuum_gap_thick = 10.0; // µm

    // SOR parameters
    const double omega = 1.8; // Relaxation factor
    const double tolerance = 1e-4; // Convergence tolerance
    const int max_iterations = 500000;

    // Material properties (relative permittivity)
    const double eps_si = 11.7;
    const double eps_vac = 1.0;

    // Boundary conditions
    const double V_left = 0.0;  // Volts
    const double V_right = -1000.0; // Volts

    // Create output folder
    const std::string output_folder = "geometria_Denti_uguali";

    // Attempt to create the output directory if it doesn't exist
    struct STAT_STRUCT info;
    if (STAT_FUNC(output_folder.c_str(), &info) != 0) { // Check if directory exists
        if (MKDIR(output_folder.c_str()) == 0) {
            std::cout << "Output directory '" << output_folder << "' created." << std::endl;
        } else {
            std::cerr << "Error: Could not create output directory '" << output_folder << "'. Please create it manually." << std::endl;
            // Optionally, exit if directory creation is critical and failed
            // return 1; 
        }
    } else if (!(info.st_mode & S_IFDIR)) { // Check if it's a directory
        std::cerr << "Error: '" << output_folder << "' exists but is not a directory. Please remove it or rename it." << std::endl;
        // Optionally, exit
        // return 1;
    } else {
        std::cout << "Output directory '" << output_folder << "' already exists." << std::endl;
    }
    // The following line requires C++17 (or later) and the <filesystem> header.
    // It was commented out to resolve compilation errors if not using a C++17 compliant compiler/flags.
    // An alternative pre-C++17 directory creation attempt has been added above.
    // std::filesystem::create_directory(output_folder);

    // Save geometry parameters before extensive calculations
    saveGeometryParamsToCSV(output_folder + "/geometry_params.csv", 
                            h, 
                            x_free_space, x_structure_len, 
                            y_si_base_height, y_teeth_height, 
                            y_vacuum_gap_thick,
                            x_teeth_width, x_teeth_spacing,
                            H_total);

    // --- Grid Setup ---
    const int Nx = static_cast<int>(L_total / h) + 1;
    const int Ny = static_cast<int>(H_total / h) + 1;

    std::vector<double> x_coords(Nx);
    std::vector<double> y_coords(Ny);
    for (int i = 0; i < Nx; ++i) x_coords[i] = i * h;
    for (int j = 0; j < Ny; ++j) y_coords[j] = j * h;

    std::vector<std::vector<double>> V(Nx, std::vector<double>(Ny, 0.0));
    std::vector<std::vector<double>> eps_r(Nx, std::vector<double>(Ny, eps_vac));
    std::vector<std::vector<bool>> fixed_potential_mask(Nx, std::vector<bool>(Ny, false));

    // --- Define Material Regions ---
    const int idx_x_struct_start = static_cast<int>(x_free_space / h);
    const int idx_x_struct_end = static_cast<int>((x_free_space + x_structure_len) / h);

    // Y-indices for the layers
    const int idx_y_bot_si_base_end = static_cast<int>(y_si_base_height / h);
    const int idx_y_bot_si_teeth_end = static_cast<int>((y_si_base_height + y_teeth_height) / h);
    
    const int idx_y_vac_end = static_cast<int>((y_si_base_height + y_teeth_height + y_vacuum_gap_thick) / h);
    // This is where the top Si structure (teeth pointing down) begins.
    
    const int idx_y_top_si_teeth_start = idx_y_vac_end; // Start of top Si teeth (bottom edge of top Si structure)
    const int idx_y_top_si_base_start = static_cast<int>((y_si_base_height + y_teeth_height + y_vacuum_gap_thick + y_teeth_height) / h);
    // The top Si base goes up to Ny-1

    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            eps_r[i][j] = eps_vac; // Default to vacuum

            if (i >= idx_x_struct_start && i <= idx_x_struct_end) {
                // Bottom Silicon Layer
                // Base
                if (j >= 0 && j <= idx_y_bot_si_base_end) {
                    eps_r[i][j] = eps_si;
                }
                // Teeth (on top of the base)
                if (j > idx_y_bot_si_base_end && j <= idx_y_bot_si_teeth_end) {
                    double x_coord_in_structure = (i - idx_x_struct_start) * h;
                    // Ensure x_coord_in_structure is non-negative for fmod
                    if (x_coord_in_structure < 0) x_coord_in_structure = 0; 
                    double pos_in_period = fmod(x_coord_in_structure, tooth_period);
                    if (pos_in_period < x_teeth_width) {
                        eps_r[i][j] = eps_si;
                    }
                }

                // Top Silicon Layer
                // Base (topmost part)
                if (j >= idx_y_top_si_base_start && j < Ny) {
                    eps_r[i][j] = eps_si;
                }
                // Teeth (hanging below the base)
                if (j >= idx_y_top_si_teeth_start && j < idx_y_top_si_base_start) {
                     double x_coord_in_structure = (i - idx_x_struct_start) * h;
                     if (x_coord_in_structure < 0) x_coord_in_structure = 0;
                     double pos_in_period = fmod(x_coord_in_structure, tooth_period);
                     if (pos_in_period < x_teeth_width) {
                        eps_r[i][j] = eps_si;
                    }
                }
            }
        }
    }
    
    std::cout << "Grid size: Nx=" << Nx << ", Ny=" << Ny << std::endl;
    std::cout << "Structure x-indices: " << idx_x_struct_start << " to " << idx_x_struct_end << std::endl;
    std::cout << "Bottom Si base y-indices: 0 to " << idx_y_bot_si_base_end << std::endl;
    std::cout << "Bottom Si teeth y-indices: " << idx_y_bot_si_base_end + 1 << " to " << idx_y_bot_si_teeth_end << std::endl;
    std::cout << "Vacuum gap y-indices (approx): " << idx_y_bot_si_teeth_end + 1 << " to " << idx_y_top_si_teeth_start -1 << std::endl;
    std::cout << "Top Si teeth y-indices: " << idx_y_top_si_teeth_start << " to " << idx_y_top_si_base_start - 1 << std::endl;
    std::cout << "Top Si base y-indices: " << idx_y_top_si_base_start << " to " << Ny - 1 << std::endl;


    // --- Boundary Conditions ---
    // Apply to all parts of the silicon structure at the left and right ends of the defined structure length
    for (int j = 0; j < Ny; ++j) {
        // Check if the point (idx_x_struct_start, j) is silicon
        if (eps_r[idx_x_struct_start][j] == eps_si) {
            V[idx_x_struct_start][j] = V_left;
            fixed_potential_mask[idx_x_struct_start][j] = true;
        }
        // Check if the point (idx_x_struct_end, j) is silicon
        if (eps_r[idx_x_struct_end][j] == eps_si) {
            V[idx_x_struct_end][j] = V_right;
            fixed_potential_mask[idx_x_struct_end][j] = true;
        }
    }

    // --- SOR Iteration ---
    for (int iteration = 0; iteration < max_iterations; ++iteration) {
        double max_diff_iter = 0.0; // Reset for each iteration, OpenMP reduction will use this

        // Phase 1: Update "red" points ((i+j) % 2 == 0)
        #pragma omp parallel for reduction(max:max_diff_iter)
        for (int i = 1; i < Nx - 1; ++i) {
            for (int j = 1; j < Ny - 1; ++j) {
                if ((i + j) % 2 == 0) { // "Red" points
                    if (fixed_potential_mask[i][j]) {
                        continue;
                    }

                    double eps_ij = eps_r[i][j];
                    double eps_e = (eps_ij + eps_r[i+1][j]) / 2.0;
                    double eps_w = (eps_ij + eps_r[i-1][j]) / 2.0;
                    double eps_n = (eps_ij + eps_r[i][j+1]) / 2.0;
                    double eps_s = (eps_ij + eps_r[i][j-1]) / 2.0;
                    
                    double sum_eps = eps_e + eps_w + eps_n + eps_s;
                    if (sum_eps == 0) continue;

                    double val_GS = (eps_e * V[i+1][j] + eps_w * V[i-1][j] +
                                     eps_n * V[i][j+1] + eps_s * V[i][j-1]) / sum_eps;
                    
                    double v_old_ij = V[i][j];
                    V[i][j] = (1.0 - omega) * v_old_ij + omega * val_GS;
                    double current_point_diff = std::abs(V[i][j] - v_old_ij);
                    if (current_point_diff > max_diff_iter) {
                        max_diff_iter = current_point_diff;
                    }
                }
            }
        }

        // Phase 2: Update "black" points ((i+j) % 2 != 0)
        #pragma omp parallel for reduction(max:max_diff_iter)
        for (int i = 1; i < Nx - 1; ++i) {
            for (int j = 1; j < Ny - 1; ++j) {
                if ((i + j) % 2 != 0) { // "Black" points
                    if (fixed_potential_mask[i][j]) {
                        continue;
                    }

                    double eps_ij = eps_r[i][j];
                    double eps_e = (eps_ij + eps_r[i+1][j]) / 2.0;
                    double eps_w = (eps_ij + eps_r[i-1][j]) / 2.0;
                    double eps_n = (eps_ij + eps_r[i][j+1]) / 2.0;
                    double eps_s = (eps_ij + eps_r[i][j-1]) / 2.0;
                    
                    double sum_eps = eps_e + eps_w + eps_n + eps_s;
                    if (sum_eps == 0) continue;

                    double val_GS = (eps_e * V[i+1][j] + eps_w * V[i-1][j] +
                                     eps_n * V[i][j+1] + eps_s * V[i][j-1]) / sum_eps;
                    
                    double v_old_ij = V[i][j];
                    V[i][j] = (1.0 - omega) * v_old_ij + omega * val_GS;
                    double current_point_diff = std::abs(V[i][j] - v_old_ij);
                    if (current_point_diff > max_diff_iter) {
                        max_diff_iter = current_point_diff;
                    }
                }
            }
        }

        // Apply Neumann BCs to outer simulation boundaries (dV/dn = 0 -> V_boundary = V_neighbor_inside)
        #pragma omp parallel for
        for (int j_idx = 0; j_idx < Ny; ++j_idx) {
            if (!fixed_potential_mask[0][j_idx]) V[0][j_idx] = V[1][j_idx];
            if (!fixed_potential_mask[Nx - 1][j_idx]) V[Nx - 1][j_idx] = V[Nx - 2][j_idx];
        }
        #pragma omp parallel for
        for (int i_idx = 0; i_idx < Nx; ++i_idx) {
            if (!fixed_potential_mask[i_idx][0]) V[i_idx][0] = V[i_idx][1];
            if (!fixed_potential_mask[i_idx][Ny - 1]) V[i_idx][Ny - 1] = V[i_idx][Ny - 2];
        }
        
        // Re-apply fixed potentials (important if fixed BCs are on the Neumann boundary itself)
        // This is generally good practice, though in this specific setup, fixed BCs are not on the outermost edges.
        for (int i = 0; i < Nx; ++i) {
            for (int j = 0; j < Ny; ++j) {
                if (fixed_potential_mask[i][j]) {
                    // The value is already set, this ensures it's not changed by Neumann logic if on edge.
                    // For this problem, fixed points are idx_x_struct_start/end which are not 0 or Nx-1.
                }
            }
        }


        if ((iteration + 1) % 100 == 0) {
            std::cout << "Iteration " << iteration + 1 << ", Max Potential Change: " << std::scientific << max_diff_iter << std::fixed << std::endl;
        }

        if (max_diff_iter < tolerance) {
            std::cout << "Converged after " << iteration + 1 << " iterations." << std::endl;
            break;
        }
        if (iteration == max_iterations - 1) {
             std::cout << "Max iterations (" << max_iterations << ") reached. Max diff: " << std::scientific << max_diff_iter << std::fixed << std::endl;
        }
    }

    // --- Calculate Electric Field ---
    std::vector<std::vector<double>> Ex(Nx, std::vector<double>(Ny, 0.0));
    std::vector<std::vector<double>> Ey(Nx, std::vector<double>(Ny, 0.0));

    for (int i = 1; i < Nx - 1; ++i) {
        for (int j = 0; j < Ny; ++j) {
            Ex[i][j] = -(V[i+1][j] - V[i-1][j]) / (2 * h);
        }
    }
    for (int i = 0; i < Nx; ++i) {
        for (int j = 1; j < Ny - 1; ++j) {
            Ey[i][j] = -(V[i][j+1] - V[i][j-1]) / (2 * h);
        }
    }

    // Forward/backward difference for boundaries
    for (int j = 0; j < Ny; ++j) {
        Ex[0][j] = -(V[1][j] - V[0][j]) / h;
        Ex[Nx-1][j] = -(V[Nx-1][j] - V[Nx-2][j]) / h;
    }
    for (int i = 0; i < Nx; ++i) {
        Ey[i][0] = -(V[i][1] - V[i][0]) / h;
        Ey[i][Ny-1] = -(V[i][Ny-1] - V[i][Ny-2]) / h;
    }

    // --- Output Results to CSV ---
    saveToCSV(V, output_folder + "/potential.csv");
    saveToCSV(Ex, output_folder + "/electric_field_x.csv");
    saveToCSV(Ey, output_folder + "/electric_field_y.csv");
    saveToCSV(eps_r, output_folder + "/permittivity.csv"); // Save permittivity map for verification/plotting
    saveCoordinatesToCSV(x_coords, output_folder + "/x_coordinates.csv");
    saveCoordinatesToCSV(y_coords, output_folder + "/y_coordinates.csv");

    std::cout << "\n--- Results Summary ---" << std::endl;
    std::cout << "Potential V, Electric fields Ex, Ey, Permittivity eps_r, and coordinates saved to CSV files." << std::endl;
    std::cout << "You can use external tools (e.g., Python with Matplotlib, Gnuplot, Excel) to plot these CSV files." << std::endl;

    return 0;
}
