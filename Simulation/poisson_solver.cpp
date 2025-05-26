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

// Helper struct to return info about a point in the structure
struct PointMaterialInfo {
    bool is_silicon_tooth_region; // True if x_target_rel falls into a tooth's x-span
    double current_dig_depth_at_x; // Dig depth if x_target_rel is in a space's x-span
    double current_tooth_height_at_x; // Height of the tooth if in a tooth region
};

// Helper function to determine material characteristics based on x-position
PointMaterialInfo get_point_x_info(double x_target_rel, double total_struct_len,
                                   double tooth_w,
                                   double initial_space_w, double space_w_increment,
                                   double initial_dig_d, double dig_d_increment,
                                   double si_base_h,
                                   double initial_tooth_h, double tooth_h_decrement) {
    PointMaterialInfo info = {false, 0.0, 0.0}; // Default: not a tooth, no digging, no tooth height
    double current_x_pos = 0.0;
    double current_space_w = initial_space_w;
    double current_dig_d = initial_dig_d;
    double current_tooth_h = initial_tooth_h;
    bool segment_is_tooth = true; // Assume structure starts with a tooth at its left edge

    // Ensure x_target_rel is within the structure's bounds for processing
    if (x_target_rel < 0 || x_target_rel >= total_struct_len) {
        return info; 
    }

    while (current_x_pos < total_struct_len) {
        if (segment_is_tooth) {
            double tooth_end_x = current_x_pos + tooth_w;
            if (x_target_rel >= current_x_pos && x_target_rel < tooth_end_x) {
                info.is_silicon_tooth_region = true;
                info.current_dig_depth_at_x = 0.0; // No digging in teeth
                info.current_tooth_height_at_x = current_tooth_h;
                return info;
            }
            current_x_pos = tooth_end_x;
            segment_is_tooth = false; // Next segment is a space
        } else { // Segment is space
            double space_end_x = current_x_pos + current_space_w;
            if (x_target_rel >= current_x_pos && x_target_rel < space_end_x) {
                info.is_silicon_tooth_region = false;
                info.current_dig_depth_at_x = current_dig_d;
                info.current_tooth_height_at_x = 0.0; // Not in a tooth
                return info;
            }
            current_x_pos = space_end_x;
            // Update parameters for the next space and subsequent tooth
            current_space_w += space_w_increment;
            current_dig_d += dig_d_increment;
            if (current_dig_d > si_base_h) { // Cap dig depth
                current_dig_d = si_base_h;
            }
            current_tooth_h = std::max(0.0, current_tooth_h - tooth_h_decrement); // Decrease height for next tooth
            segment_is_tooth = true; // Next segment is a tooth
        }
    }
    return info;
}

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
                             double y_sibh, 
                             double init_y_th, double y_th_dec, // Modified y_teeth_height
                             double y_vgt,
                             double const_x_tooth_w,
                             double init_x_space_w, double x_space_w_inc,
                             double init_y_dig_d, double y_dig_d_inc,
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
    // outfile << "y_teeth_height," << y_th << std::endl; // Old
    outfile << "initial_y_teeth_height," << init_y_th << std::endl;
    outfile << "y_teeth_height_decrement," << y_th_dec << std::endl;
    outfile << "y_vacuum_gap_thick," << y_vgt << std::endl;
    outfile << "x_teeth_width," << const_x_tooth_w << std::endl; 
    outfile << "initial_x_spacing_width," << init_x_space_w << std::endl;
    outfile << "x_spacing_width_increment," << x_space_w_inc << std::endl;
    outfile << "initial_y_spacing_dig_depth," << init_y_dig_d << std::endl;
    outfile << "y_spacing_dig_depth_increment," << y_dig_d_inc << std::endl;
    outfile << "H_total," << H_tot << std::endl;
    outfile.close();
    std::cout << "Geometry parameters saved to " << filename << std::endl;
}


int main(int argc, char* argv[]) { // Modified main signature
    // --- Parameters ---
    const double h = 0.5; // Grid spacing in micrometers (µm)

    // Dimensions in µm
    const double L_total = 320.0;
  

    const double x_free_space = 10.0;
    const double x_structure_len = 300.0;
    
    // Silicon layer parameters
    const double y_si_base_height = 10.0;   // µm
    // const double y_teeth_height = 10.0;     // µm // Replaced
    const double initial_y_teeth_height_um = 10.0; // Height of the leftmost tooth
    const double y_teeth_height_decrement_um = 1.0; // Decrease in height for each tooth to the right
    
    // New geometry parameters
    const double x_teeth_width_um = 10.0; // Constant width for all teeth
    
    const double initial_x_spacing_width_um = 10.0; // Width of the leftmost space
    const double x_spacing_width_increment_um = 2.0; // Increment for space width (L to R)
    
    const double initial_y_spacing_dig_depth_um = 0.0; // Digging depth for the leftmost space (0 means no initial dig)
    const double y_spacing_dig_depth_increment_um = 1.0; // Increment for digging depth (L to R)

    const double y_vacuum_gap_thick = 10.0; // µm
      const double H_total = 55 ; // Total height of the structure
    // SOR parameters
    const double omega = 1.8; // Relaxation factor
    const double tolerance = 1e-5; // Convergence tolerance
    const int max_iterations = 500000;

    // Material properties (relative permittivity)
    const double eps_si = 11.7;
    const double eps_vac = 1.0;

    // Boundary conditions
    const double V_left = 0.0;  // Volts
    const double V_right = -1000.0; // Volts

    // Create output folder
    std::string output_folder_name = "geometria_Denti_sfasati_profondi_5um_default"; // Default name
    if (argc > 1) {
        output_folder_name = argv[1]; // Use the first command-line argument as folder name
        std::cout << "Output folder specified: " << output_folder_name << std::endl;
    } else {
        std::cout << "No output folder specified, using default: " << output_folder_name << std::endl;
    }
    const std::string output_folder = output_folder_name;


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
                            y_si_base_height, 
                            initial_y_teeth_height_um, y_teeth_height_decrement_um,
                            y_vacuum_gap_thick,
                            x_teeth_width_um, 
                            initial_x_spacing_width_um, x_spacing_width_increment_um,
                            initial_y_spacing_dig_depth_um, y_spacing_dig_depth_increment_um,
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
    // const int idx_x_struct_start = static_cast<int>(x_free_space / h); // Not directly used in loop below
    // const int idx_x_struct_end = static_cast<int>((x_free_space + x_structure_len) / h); // Not directly used

    // Y-indices for the layers (conceptual, not directly used in the refined loop)
    // const int idx_y_bot_si_base_end = static_cast<int>(y_si_base_height / h);
    // const int idx_y_bot_si_teeth_end = static_cast<int>((y_si_base_height + y_teeth_height) / h);
    // const int idx_y_vac_end = static_cast<int>((y_si_base_height + y_teeth_height + y_vacuum_gap_thick) / h);
    // const int idx_y_top_si_teeth_start = idx_y_vac_end; 
    // const int idx_y_top_si_base_start = static_cast<int>((y_si_base_height + y_teeth_height + y_vacuum_gap_thick + y_teeth_height) / h);

    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            eps_r[i][j] = eps_vac; // Default to vacuum

            double x_abs = i * h;
            double y_abs = j * h;

            // Check if the point is within the x-range of the structured region
            if (x_abs >= x_free_space && x_abs < (x_free_space + x_structure_len)) {
                double x_coord_in_structure = x_abs - x_free_space;

                PointMaterialInfo x_info = get_point_x_info(x_coord_in_structure, x_structure_len,
                                                            x_teeth_width_um,
                                                            initial_x_spacing_width_um, x_spacing_width_increment_um,
                                                            initial_y_spacing_dig_depth_um, y_spacing_dig_depth_increment_um,
                                                            y_si_base_height,
                                                            initial_y_teeth_height_um, y_teeth_height_decrement_um);
                
                // Define y-boundaries for clarity based on current tooth height
                double current_tooth_h = x_info.current_tooth_height_at_x; // This is 0 if not in a tooth region

                double y_bottom_base_top = y_si_base_height;
                double y_bottom_teeth_top = y_si_base_height + (x_info.is_silicon_tooth_region ? current_tooth_h : 0.0);
                
                // The vacuum gap top is now dynamic based on the bottom tooth height
                // If not in a tooth region, the "effective" tooth height for gap calculation is 0 at that x.
                // However, the structure definition below handles this by checking is_silicon_tooth_region.
                // The y_vacuum_gap_thick parameter defines the nominal separation.
                // The top of the vacuum gap is H_total - y_si_base_height - (top_tooth_height)
                // Let's define effective boundaries for the current x_info
                
                double y_actual_bottom_teeth_tip = y_si_base_height + (x_info.is_silicon_tooth_region ? current_tooth_h : 0.0);
                double y_actual_top_teeth_tip = H_total - y_si_base_height - (x_info.is_silicon_tooth_region ? current_tooth_h : 0.0);


                // Bottom Layer Processing
                // Bottom Base region: 0 <= y_abs < y_bottom_base_top
                if (y_abs >= 0 && y_abs < y_bottom_base_top) {
                    if (x_info.is_silicon_tooth_region) {
                        eps_r[i][j] = eps_si; 
                    } else { 
                        if (y_abs < (y_bottom_base_top - x_info.current_dig_depth_at_x)) {
                            eps_r[i][j] = eps_si; 
                        } 
                    }
                }
                // Bottom Teeth region: y_bottom_base_top <= y_abs < y_actual_bottom_teeth_tip
                else if (x_info.is_silicon_tooth_region && current_tooth_h > 0 &&
                         y_abs >= y_bottom_base_top && y_abs < y_actual_bottom_teeth_tip) {
                    eps_r[i][j] = eps_si; 
                }
                // Vacuum Gap: y_actual_bottom_teeth_tip <= y_abs < y_actual_top_teeth_tip
                // This region remains eps_vac by default if not part of top/bottom teeth.

                // Top Teeth Region: y_actual_top_teeth_tip <= y_abs < (H_total - y_si_base_height)
                else if (x_info.is_silicon_tooth_region && current_tooth_h > 0 &&
                         y_abs >= y_actual_top_teeth_tip && y_abs < (H_total - y_si_base_height)) {
                     eps_r[i][j] = eps_si; // This is a top tooth
                }
                // Top Base Region: (H_total - y_si_base_height) <= y_abs < H_total
                else if (y_abs >= (H_total - y_si_base_height) && y_abs < H_total) {
                    double y_top_base_bottom_surface = H_total - y_si_base_height;
                    if (x_info.is_silicon_tooth_region) {
                        eps_r[i][j] = eps_si; 
                    } else { 
                        // Digging is from the "bottom" of the top base (y_top_base_bottom_surface), upwards.
                        // Silicon remains if y_abs is beyond (i.e. above) the dug portion.
                        if (y_abs >= (y_top_base_bottom_surface + x_info.current_dig_depth_at_x)) {
                            eps_r[i][j] = eps_si; 
                        } 
                    }
                }
            }
        }
    }
    
    std::cout << "Grid size: Nx=" << Nx << ", Ny=" << Ny << std::endl;
    // Remove old detailed couts for indices as they are less relevant with the new method
    // std::cout << "Structure x-indices: " << idx_x_struct_start << " to " << idx_x_struct_end << std::endl;
    // std::cout << "Bottom Si base y-indices: 0 to " << idx_y_bot_si_base_end << std::endl;
    // std::cout << "Bottom Si teeth y-indices: " << idx_y_bot_si_base_end + 1 << " to " << idx_y_bot_si_teeth_end << std::endl;
    // std::cout << "Vacuum gap y-indices (approx): " << idx_y_bot_si_teeth_end + 1 << " to " << idx_y_top_si_teeth_start -1 << std::endl;
    // std::cout << "Top Si teeth y-indices: " << idx_y_top_si_teeth_start << " to " << idx_y_top_si_base_start - 1 << std::endl;
    // std::cout << "Top Si base y-indices: " << idx_y_top_si_base_start << " to " << Ny - 1 << std::endl;


    // --- Boundary Conditions ---
    // Apply to all parts of the silicon structure at the left and right ends of the defined structure length
    // Need to determine the actual start and end indices of the structure for applying BCs
    const int actual_idx_x_struct_start = static_cast<int>(x_free_space / h);
    const int actual_idx_x_struct_end = static_cast<int>((x_free_space + x_structure_len - h) / h); // -h to be inclusive of last point within structure

    for (int j = 0; j < Ny; ++j) {
        // Check if the point (actual_idx_x_struct_start, j) is silicon
        if (actual_idx_x_struct_start >=0 && actual_idx_x_struct_start < Nx && eps_r[actual_idx_x_struct_start][j] == eps_si) {
            V[actual_idx_x_struct_start][j] = V_left;
            fixed_potential_mask[actual_idx_x_struct_start][j] = true;
        }
        // Check if the point (actual_idx_x_struct_end, j) is silicon
        if (actual_idx_x_struct_end >= 0 && actual_idx_x_struct_end < Nx && eps_r[actual_idx_x_struct_end][j] == eps_si) {
            V[actual_idx_x_struct_end][j] = V_right;
            fixed_potential_mask[actual_idx_x_struct_end][j] = true;
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
