#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <algorithm> // For std::min, std::max
#include <omp.h>     // For OpenMP
#include "geometry_definitions.h" // Includi il nuovo header

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

// Helper function to save a 2D vector to a CSV file with optimized I/O
void saveToCSV(const std::vector<std::vector<double>>& data, const std::string& filename) {
    std::ofstream outfile(filename, std::ios::out | std::ios::binary);
    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }
    
    // Set buffer size to improve I/O performance
    const int BUFFER_SIZE = 16384; // 16KB buffer
    std::vector<char> buffer(BUFFER_SIZE);
    outfile.rdbuf()->pubsetbuf(buffer.data(), BUFFER_SIZE);
    
    // Use string stream for faster composition
    std::ostringstream rowStream;
    rowStream << std::setprecision(10); // Set precision once
    
    // Pre-calculate sizes for efficiency
    size_t nx = data.size();
    if (nx == 0) return;
    size_t ny = data[0].size();
    if (ny == 0) return;
    
    // Write data row by row (easier to parse for analysis tools)
    for (size_t j = 0; j < ny; ++j) {
        rowStream.str(""); // Clear stream
        rowStream.clear();  // Clear flags
        
        for (size_t i = 0; i < nx; ++i) {
            rowStream << data[i][j];
            if (i < nx - 1) rowStream << ",";
        }
        rowStream << "\n";
        
        // Write the entire row at once
        outfile << rowStream.str();
    }
    
    outfile.close();
    std::cout << "Data saved to " << filename << std::endl;
}

// Helper function to save a 1D vector to a CSV file with optimized I/O
void saveCoordinatesToCSV(const std::vector<double>& coords, const std::string& filename) {
    std::ofstream outfile(filename, std::ios::out | std::ios::binary);
    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }
    
    // Set buffer size
    const int BUFFER_SIZE = 8192; // 8KB buffer
    std::vector<char> buffer(BUFFER_SIZE);
    outfile.rdbuf()->pubsetbuf(buffer.data(), BUFFER_SIZE);
    
    // Set precision once
    outfile << std::setprecision(10);
    
    // Write all data at once with a stringstream
    std::ostringstream dataStream;
    for (size_t i = 0; i < coords.size(); ++i) {
        dataStream << coords[i] << "\n";
    }
    outfile << dataStream.str();
    
    outfile.close();
    std::cout << "Coordinates saved to " << filename << std::endl;
}

// La funzione saveGeometryParamsToCSV è stata spostata in geometry_definitions.cpp e rinominata saveGeometryParams

int main(int argc, char* argv[]) { 
    // --- Command Line Argument Parsing ---
    std::string output_folder_name = "default_output";
    if (argc > 1) {
        output_folder_name = argv[1];
    } else {
        std::cout << "No output folder specified, using default: " << output_folder_name << std::endl;
    }

    std::string geometry_type_str = "piana"; // Default geometry
    if (argc > 2) {
        geometry_type_str = argv[2];
    } else {
        std::cout << "No geometry type specified, using default: " << geometry_type_str << std::endl;
    }
    
    GeometryType current_geometry_type = stringToGeometryType(geometry_type_str);
    if (current_geometry_type == GeometryType::UNKNOWN) {
        std::cerr << "Error: Unknown geometry type '" << geometry_type_str << "'. Exiting." << std::endl;
        return 1;
    }
    outfile << std::fixed << std::setprecision(10); // Ensure precision
    outfile << "h," << h_val << std::endl;
    outfile << "x_free_space," << x_fs << std::endl;
    outfile << "x_structure_len," << x_sl << std::endl;
    outfile << "y_si_layer_thick," << y_slt << std::endl;
    outfile << "y_vacuum_gap_thick," << y_vgt << std::endl;
    outfile << "H_total," << H_tot << std::endl;
    outfile.close();
    std::cout << "Geometry parameters saved to " << filename << std::endl;
}


int main() {
    // --- Parameters ---
    const double h = 0.5; // Grid spacing in micrometers (µm)

    // Dimensions in µm
    const double L_total = 320.0;
    const double H_total = 30.0;

    const double x_free_space = 10.0;
    const double x_structure_len = 300.0;
    const double y_si_layer_thick = 10.0;
    const double y_vacuum_gap_thick = 10.0;

    // SOR parameters
    const double omega = 1.8; // Relaxation factor
    const double tolerance = 1e-5; // Convergence tolerance
    const int max_iterations = 50000;

    // Material properties (relative permittivity)
    const double eps_si = 11.7;
    const double eps_vac = 1.0;

    // Dichiarazione delle strutture dei parametri
    GeometryConfig geom_config;
    PianaSpecificParams piana_geom_params;
    DentiSfasatiProfondiSpecificParams denti_geom_params;
    DentiUgualiSpecificParams du_geom_params; // Added

    // Inizializzazione dei parametri basata sul tipo di geometria
    if (current_geometry_type == GeometryType::PIANA) {
        initializePianaGeometry(geom_config, piana_geom_params, common_h_param, eps_sio2_mat, eps_vac_mat);
    } else if (current_geometry_type == GeometryType::DENTI_SFASATI_PROFONDI) {
        initializeDentiSfasatiProfondiGeometry(geom_config, denti_geom_params, common_h_param, eps_sio2_mat, eps_vac_mat);
    } else if (current_geometry_type == GeometryType::DENTI_UGUALI) { // Added
        initializeDentiUgualiGeometry(geom_config, du_geom_params, common_h_param, eps_sio2_mat, eps_vac_mat);
    }

    // Attempt to create the output directory  <- This old block is replaced by the logic above.
    // struct STAT_STRUCT info;
    // if (STAT_FUNC(output_folder.c_str(), &info) != 0) {
    // MKDIR(output_folder.c_str());
    // }

    // Save geometry parameters
    if (current_geometry_type == GeometryType::PIANA) {
        saveGeometryParams(output_folder + "/geometry_params.csv", current_geometry_type, geom_config, &piana_geom_params, nullptr, nullptr);
    } else if (current_geometry_type == GeometryType::DENTI_SFASATI_PROFONDI) {
        saveGeometryParams(output_folder + "/geometry_params.csv", current_geometry_type, geom_config, nullptr, &denti_geom_params, nullptr);    } else if (current_geometry_type == GeometryType::DENTI_UGUALI) { // Added
        saveGeometryParams(output_folder + "/geometry_params.csv", current_geometry_type, geom_config, nullptr, nullptr, &du_geom_params);
    } else if (current_geometry_type == GeometryType::DENTI_SFASATI_PROFONDI_NM) {
        saveGeometryParams(output_folder + "/geometry_params.csv", current_geometry_type, geom_config, nullptr, &denti_geom_params, nullptr);
    } else {
        std::cerr << "Error: Unsupported geometry type for saving parameters." << std::endl;
        return 1; // Exit if geometry type is unsupported
    }
    // --- Grid Setup ---
    const int Nx = static_cast<int>(geom_config.L_total / geom_config.h) + 1;
    const int Ny = static_cast<int>(geom_config.H_total_val / geom_config.h) + 1;

    // Pre-allocate vectors with correct size
    std::vector<double> x_coords(Nx);
    std::vector<double> y_coords(Ny);
    
    // Generate coordinates more efficiently
    #pragma omp parallel sections
    {
        #pragma omp section
        {
            // Fill x coordinates with a single multiplication instead of repeated multiplications
            double x = 0.0;
            for (int i = 0; i < Nx; ++i, x += geom_config.h) {
                x_coords[i] = x;
            }
        }
        
        #pragma omp section
        {
            // Fill y coordinates similarly
            double y = 0.0;
            for (int j = 0; j < Ny; ++j, y += geom_config.h) {
                y_coords[j] = y;
            }
        }
    }
    
    // Allocate main computation arrays with efficient initialization
    // Reserve memory for all arrays at once to help with cache locality
    std::vector<std::vector<double>> V(Nx, std::vector<double>(Ny, 0.0));
    std::vector<std::vector<double>> eps_r(Nx, std::vector<double>(Ny, geom_config.eps_vacuum_val));
    std::vector<std::vector<bool>> fixed_potential_mask(Nx, std::vector<bool>(Ny, false));

    // --- Define Material Regions ---
    if (current_geometry_type == GeometryType::PIANA) {
        setupPianaPermittivity(eps_r, geom_config, piana_geom_params, Nx, Ny);
    } else if (current_geometry_type == GeometryType::DENTI_SFASATI_PROFONDI) {
        setupDentiSfasatiProfondiPermittivity(eps_r, geom_config, denti_geom_params, Nx, Ny);
    } else if (current_geometry_type == GeometryType::DENTI_UGUALI) { // Added
        setupDentiUgualiPermittivity(eps_r, geom_config, du_geom_params, Nx, Ny);
    }
    
    std::cout << "Grid size: Nx=" << Nx << ", Ny=" << Ny << ", Grid spacing h=" << geom_config.h << " µm" << std::endl;

    // --- Boundary Conditions ---
    setupBoundaryConditions(V, fixed_potential_mask, eps_r, geom_config, V_left_bc, V_right_bc, Nx, Ny);

    // --- SOR Iteration ---
    double tolerance_sor = geom_config.current_tolerance;
    double max_diff_iter = 0.0;
    
    // Pre-compute 1.0 - omega_sor once outside the loop
    const double one_minus_omega = 1.0 - omega_sor;
    
    for (int iteration = 0; iteration < max_iter_sor; ++iteration) {
        max_diff_iter = 0.0;

        // Phase 1: Update "red" points
        #pragma omp parallel for reduction(max:max_diff_iter) schedule(static)
        for (int i = 1; i < Nx - 1; ++i) {
            // Cache the column values to improve memory access patterns
            double eps_i_j, eps_ip1_j, eps_im1_j, eps_i_jp1, eps_i_jm1;
            double v_i_j, v_ip1_j, v_im1_j, v_i_jp1, v_i_jm1;
            
            for (int j = 1; j < Ny - 1; ++j) {
                // Process only "red" points and non-fixed points
                if ((i + j) % 2 == 0 && !fixed_potential_mask[i][j]) {
                    // Cache values for current position and neighbors
                    eps_i_j = eps_r[i][j];
                    eps_ip1_j = eps_r[i+1][j];
                    eps_im1_j = eps_r[i-1][j];
                    eps_i_jp1 = eps_r[i][j+1];
                    eps_i_jm1 = eps_r[i][j-1];
                    
                    v_i_j = V[i][j];
                    v_ip1_j = V[i+1][j];
                    v_im1_j = V[i-1][j];
                    v_i_jp1 = V[i][j+1];
                    v_i_jm1 = V[i][j-1];
                    
                    // Pre-compute average permittivity at interfaces
                    double eps_e = (eps_i_j + eps_ip1_j) * 0.5;
                    double eps_w = (eps_i_j + eps_im1_j) * 0.5;
                    double eps_n = (eps_i_j + eps_i_jp1) * 0.5;
                    double eps_s = (eps_i_j + eps_i_jm1) * 0.5;
                    
                    double sum_eps = eps_e + eps_w + eps_n + eps_s;
                    if (sum_eps > 0) {
                        double val_GS = (eps_e * v_ip1_j + eps_w * v_im1_j +
                                        eps_n * v_i_jp1 + eps_s * v_i_jm1) / sum_eps;
                        
                        // SOR update
                        double v_new = one_minus_omega * v_i_j + omega_sor * val_GS;
                        double current_point_diff = std::abs(v_new - v_i_j);
                        
                        // Update in place
                        V[i][j] = v_new;
                        
                        max_diff_iter = std::max(max_diff_iter, current_point_diff);
                    }
                }
            }
        }

        // Phase 2: Update "black" points
        #pragma omp parallel for reduction(max:max_diff_iter) schedule(static)
        for (int i = 1; i < Nx - 1; ++i) {
            // Cache the column values
            double eps_i_j, eps_ip1_j, eps_im1_j, eps_i_jp1, eps_i_jm1;
            double v_i_j, v_ip1_j, v_im1_j, v_i_jp1, v_i_jm1;
            
            for (int j = 1; j < Ny - 1; ++j) {
                // Process only "black" points and non-fixed points
                if ((i + j) % 2 != 0 && !fixed_potential_mask[i][j]) {
                    // Cache values for current position and neighbors
                    eps_i_j = eps_r[i][j];
                    eps_ip1_j = eps_r[i+1][j];
                    eps_im1_j = eps_r[i-1][j];
                    eps_i_jp1 = eps_r[i][j+1];
                    eps_i_jm1 = eps_r[i][j-1];
                    
                    v_i_j = V[i][j];
                    v_ip1_j = V[i+1][j];
                    v_im1_j = V[i-1][j];
                    v_i_jp1 = V[i][j+1];
                    v_i_jm1 = V[i][j-1];
                    
                    // Pre-compute average permittivity at interfaces
                    double eps_e = (eps_i_j + eps_ip1_j) * 0.5;
                    double eps_w = (eps_i_j + eps_im1_j) * 0.5;
                    double eps_n = (eps_i_j + eps_i_jp1) * 0.5;
                    double eps_s = (eps_i_j + eps_i_jm1) * 0.5;
                    
                    double sum_eps = eps_e + eps_w + eps_n + eps_s;
                    if (sum_eps > 0) {
                        double val_GS = (eps_e * v_ip1_j + eps_w * v_im1_j +
                                        eps_n * v_i_jp1 + eps_s * v_i_jm1) / sum_eps;
                        
                        // SOR update
                        double v_new = one_minus_omega * v_i_j + omega_sor * val_GS;
                        double current_point_diff = std::abs(v_new - v_i_j);
                        
                        // Update in place
                        V[i][j] = v_new;
                        
                        max_diff_iter = std::max(max_diff_iter, current_point_diff);
                    }
                }
            }
        }

        // Apply Neumann BCs to outer simulation boundaries more efficiently
        // Left and right boundaries
        #pragma omp parallel for
        for (int j = 0; j < Ny; ++j) {
            // Only update if not a fixed point
            if (!fixed_potential_mask[0][j]) {
                V[0][j] = V[1][j];
            }
            if (!fixed_potential_mask[Nx - 1][j]) {
                V[Nx - 1][j] = V[Nx - 2][j];
            }
        }
        
        // Bottom and top boundaries
        #pragma omp parallel for
        for (int i = 0; i < Nx; ++i) {
            // Only update if not a fixed point
            if (!fixed_potential_mask[i][0]) {
                V[i][0] = V[i][1];
            }
            if (!fixed_potential_mask[i][Ny - 1]) {
                V[i][Ny - 1] = V[i][Ny - 2];
            }
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

        if (max_diff_iter < tolerance_sor) { // Usa tolerance_sor
            std::cout << "Converged after " << iteration + 1 << " iterations." << std::endl;
            break;
        }
        if (iteration == max_iter_sor - 1) { // Usa max_iter_sor
             std::cout << "Max iterations (" << max_iter_sor << ") reached. Max diff: " << std::scientific << max_diff_iter << std::fixed << std::endl;
        }
    }

    // --- Calculate Electric Field more efficiently ---
    // Declare electric field vectors
    std::vector<std::vector<double>> Ex(Nx, std::vector<double>(Ny, 0.0));
    std::vector<std::vector<double>> Ey(Nx, std::vector<double>(Ny, 0.0));
    
    // Pre-compute the reciprocal of 2*h to avoid division in loops
    const double inv_2h = 1.0 / (2.0 * geom_config.h);
    const double inv_h = 1.0 / geom_config.h;
    
    // Central difference for interior points
    #pragma omp parallel for schedule(static)
    for (int i = 1; i < Nx - 1; ++i) {
        for (int j = 0; j < Ny; ++j) {
            // Calculate x-component with central difference
            Ex[i][j] = -(V[i+1][j] - V[i-1][j]) * inv_2h;
        }
    }
    
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < Nx; ++i) {
        for (int j = 1; j < Ny - 1; ++j) {
            // Calculate y-component with central difference
            Ey[i][j] = -(V[i][j+1] - V[i][j-1]) * inv_2h;
        }
    }

    // Forward/backward difference for boundaries
    #pragma omp parallel for
    for (int j = 0; j < Ny; ++j) {
        // Left boundary - forward difference
        Ex[0][j] = -(V[1][j] - V[0][j]) * inv_h;
        // Right boundary - backward difference
        Ex[Nx-1][j] = -(V[Nx-1][j] - V[Nx-2][j]) * inv_h;
    }
    
    // Explicitly set Ey to zero at boundaries (this can be vectorized)
    #pragma omp parallel for simd
    for (int i = 0; i < Nx; ++i) {
        Ey[i][0] = 0.0;
        Ey[i][Ny-1] = 0.0;
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