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
    rowStream << std::scientific << std::setprecision(15); // Scientific notation with high precision
    
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
    
    std::cout << "Selected output folder: " << output_folder_name << std::endl;
    std::cout << "Selected geometry type: " << geometryTypeToString(current_geometry_type) << std::endl;

    const std::string base_simulation_dir = "./";
    const std::string output_folder = base_simulation_dir + output_folder_name; // Full path to final output directory

    // --- Create Output Directories Recursively ---
    // This section ensures that if output_folder_name is "A/B", both "base_simulation_dir/A" 
    // and "base_simulation_dir/A/B" are created if they don't exist.
    std::string path_being_built = base_simulation_dir;
    // Ensure base_simulation_dir exists (typically a given, but good for full robustness check if needed)
    // For this solution, we assume base_simulation_dir exists and is writable.
    // We also ensure base_simulation_dir ends with a slash for proper concatenation.
    if (!path_being_built.empty() && path_being_built.back() != '/') {
        path_being_built += '/';
    }

    if (output_folder_name.empty()) {
        // If output_folder_name is empty, output_folder is just base_simulation_dir.
        // Check if base_simulation_dir itself exists and is a directory.
        struct STAT_STRUCT info;
        if (STAT_FUNC(base_simulation_dir.c_str(), &info) != 0) {
            std::cerr << "Error: Base simulation directory " << base_simulation_dir << " does not exist. Exiting." << std::endl;
            return 1;
        } else if (!(info.st_mode & S_IFDIR)) {
            std::cerr << "Error: Path " << base_simulation_dir << " is not a directory. Exiting." << std::endl;
            return 1;
        }
    } else {
        size_t start = 0;
        size_t end = 0;
        std::string current_segment_to_process = output_folder_name;

        // Remove trailing slash from current_segment_to_process if present, to avoid empty last component
        if (!current_segment_to_process.empty() && current_segment_to_process.back() == '/') {
            current_segment_to_process.pop_back();
        }
        
        if (!current_segment_to_process.empty()) { // Proceed only if there are segments to create
            while (true) {
                // Find the next path separator
                end = current_segment_to_process.find('/', start);
                
                // Extract the current directory component
                std::string dir_component;
                if (end == std::string::npos) { // Last component
                    dir_component = current_segment_to_process.substr(start);
                } else { // Intermediate component
                    dir_component = current_segment_to_process.substr(start, end - start);
                }

                if (!dir_component.empty()) { // Avoid issues with empty components (e.g., "A//B")
                    // Append the component to the path being built
                    // path_being_built should already have a trailing slash from base_simulation_dir or previous iteration
                    path_being_built += dir_component;

                    struct STAT_STRUCT info;
                    if (STAT_FUNC(path_being_built.c_str(), &info) != 0) { // If path component does not exist
                        if (MKDIR(path_being_built.c_str()) != 0) { // Attempt to create it
                            std::cerr << "Error: Could not create directory " << path_being_built << ". Exiting." << std::endl;
                            return 1; // Critical error, cannot proceed
                        }
                        std::cout << "Created directory: " << path_being_built << std::endl;
                    } else if (!(info.st_mode & S_IFDIR)) { // If path exists but is not a directory
                        std::cerr << "Error: Path " << path_being_built << " exists but is not a directory. Exiting." << std::endl;
                        return 1; // Critical error
                    }
                    // Add slash for the next component, only if it's not the last one
                    if (end != std::string::npos) {
                         path_being_built += '/';
                    }
                }
                
                if (end == std::string::npos) {
                    break; // All components processed
                }
                start = end + 1; // Move to the character after the slash
                // If start is at or beyond the end of the string (e.g. "A/B/"), no more components.
                if (start >= current_segment_to_process.length()) {
                    break;
                }
            }
        }
    }
    // At this point, the full 'output_folder' path should exist.

    // --- Common Parameters ---
    const double common_h_param = 0.5; // Grid spacing in micrometers (µm)
    const double V_left_bc = 0.0;  // Volts
    const double V_right_bc = -50000.0; // Volts (50 kV for high-field acceleration)
    
    // Optimal SOR parameter calculation for faster convergence
    // For 2D Laplace with grid spacing h, optimal omega ≈ 2/(1 + π*h/L)
    // We use a slightly more conservative value for stability with material interfaces
    const double omega_sor = 1.92; // Optimized relaxation factor (was 1.8) - closer to optimal ~1.95
    const int max_iter_sor = 5000000; // 
    
    // Material properties
    const double eps_sio2_mat = 3.9;
    const double eps_si_mat = 11.7;
    const double eps_vac_mat = 1.0;

    // Dichiarazione delle strutture dei parametri
    GeometryConfig geom_config;
    PianaSpecificParams piana_geom_params;
    PianaRastremataSpecificParams piana_rastremata_params;
    PianaVariabileSpecificParams piana_variabile_params;
    DentiSfasatiProfondiSpecificParams denti_geom_params;
    DentiUgualiSpecificParams du_geom_params; // Added

    // Inizializzazione dei parametri basata sul tipo di geometria
    if (current_geometry_type == GeometryType::PIANA) {
        initializePianaGeometry(geom_config, piana_geom_params, common_h_param, eps_sio2_mat, eps_vac_mat);
    } else if (current_geometry_type == GeometryType::PIANA_RASTREMATA) {
        initializePianaRastremataGeometry(geom_config, piana_rastremata_params, common_h_param, eps_sio2_mat, eps_vac_mat);
    } else if (current_geometry_type == GeometryType::PIANA_VARIABILE) {
        initializePianaVariabileGeometry(geom_config, piana_variabile_params, common_h_param, eps_sio2_mat, eps_vac_mat);
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
        saveGeometryParams(output_folder + "/geometry_params.csv", current_geometry_type, geom_config, &piana_geom_params, nullptr, nullptr, nullptr, nullptr);
    } else if (current_geometry_type == GeometryType::PIANA_RASTREMATA) {
        saveGeometryParams(output_folder + "/geometry_params.csv", current_geometry_type, geom_config, nullptr, &piana_rastremata_params, nullptr, nullptr, nullptr);
    } else if (current_geometry_type == GeometryType::PIANA_VARIABILE) {
        saveGeometryParams(output_folder + "/geometry_params.csv", current_geometry_type, geom_config, nullptr, nullptr, &piana_variabile_params, nullptr, nullptr);
    } else if (current_geometry_type == GeometryType::DENTI_SFASATI_PROFONDI) {
        saveGeometryParams(output_folder + "/geometry_params.csv", current_geometry_type, geom_config, nullptr, nullptr, nullptr, &denti_geom_params, nullptr);
    } else if (current_geometry_type == GeometryType::DENTI_UGUALI) { // Added
        saveGeometryParams(output_folder + "/geometry_params.csv", current_geometry_type, geom_config, nullptr, nullptr, nullptr, nullptr, &du_geom_params);
    } else if (current_geometry_type == GeometryType::DENTI_SFASATI_PROFONDI_NM) {
        saveGeometryParams(output_folder + "/geometry_params.csv", current_geometry_type, geom_config, nullptr, nullptr, nullptr, &denti_geom_params, nullptr);
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
    } else if (current_geometry_type == GeometryType::PIANA_RASTREMATA) {
        setupPianaRastremataPermittivity(eps_r, geom_config, piana_rastremata_params, Nx, Ny);
    } else if (current_geometry_type == GeometryType::PIANA_VARIABILE) {
        setupPianaVariabilePermittivity(eps_r, geom_config, piana_variabile_params, Nx, Ny);
    } else if (current_geometry_type == GeometryType::DENTI_SFASATI_PROFONDI) {
        setupDentiSfasatiProfondiPermittivity(eps_r, geom_config, denti_geom_params, Nx, Ny);
    } else if (current_geometry_type == GeometryType::DENTI_UGUALI) { // Added
        setupDentiUgualiPermittivity(eps_r, geom_config, du_geom_params, Nx, Ny);
    }
    
    // Add aluminum plates at the boundaries (electrodes)
    addAluminumPlates(eps_r, geom_config, Nx, Ny);
    
    std::cout << "Grid size: Nx=" << Nx << ", Ny=" << Ny << ", Grid spacing h=" << geom_config.h << " µm" << std::endl;

    // --- Boundary Conditions ---
    setupBoundaryConditions(V, fixed_potential_mask, eps_r, geom_config, V_left_bc, V_right_bc, Nx, Ny);
    
    // Apply initial Neumann boundary conditions on top and bottom for symmetry
    std::cout << "Applying Neumann boundary conditions on top/bottom for symmetry..." << std::endl;
    #pragma omp parallel for
    for (int i = 0; i < Nx; ++i) {
        // Bottom boundary: dV/dy = 0
        if (!fixed_potential_mask[i][0]) {
            V[i][0] = V[i][1];
        }
        // Top boundary: dV/dy = 0
        if (!fixed_potential_mask[i][Ny - 1]) {
            V[i][Ny - 1] = V[i][Ny - 2];
        }
    }

    // --- SOR Iteration ---
    double tolerance_sor = geom_config.current_tolerance;
    double max_diff_iter = 0.0;
    
    // Sistema di convergenza avanzato: controlla ultimi 10 valori
    const int convergence_window = 10;
    std::vector<double> last_diffs(convergence_window, 1e6); // Inizializza con valori alti
    int no_improvement_count = 0;
    const double relative_improvement_threshold = 1e-12; // Soglia per miglioramento relativo (0.01%)
    
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
        // Left and right boundaries - improved Neumann conditions
        #pragma omp parallel for
        for (int j = 0; j < Ny; ++j) {
            // Only update if not a fixed point
            if (!fixed_potential_mask[0][j]) {
                // Use linear extrapolation for smoother boundary
                if (j > 0 && j < Ny - 1) {
                    V[0][j] = V[1][j]; // Standard Neumann
                } else {
                    V[0][j] = V[1][j];
                }
            }
            if (!fixed_potential_mask[Nx - 1][j]) {
                // Use linear extrapolation for smoother boundary
                if (j > 0 && j < Ny - 1) {
                    V[Nx - 1][j] = V[Nx - 2][j]; // Standard Neumann
                } else {
                    V[Nx - 1][j] = V[Nx - 2][j];
                }
            }
        }
        
        // Bottom and top boundaries - improved Neumann conditions
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
            std::cout << "Iteration " << iteration + 1 << ", Max Potential Change: " << std::scientific << max_diff_iter << std::fixed;
            if (iteration >= convergence_window) {
                std::cout << " [No improvement count: " << no_improvement_count << "]";
            }
            std::cout << std::endl;
        }
        
        // Warning if convergence is very slow
        if (iteration == 10000 && max_diff_iter > tolerance_sor * 100) {
            std::cout << "WARNING: Slow convergence detected at 10k iterations. Current error: " 
                      << std::scientific << max_diff_iter << " vs target: " << tolerance_sor << std::fixed << std::endl;
            std::cout << "Consider: 1) Reducing voltage, 2) Increasing tolerance, 3) Checking geometry" << std::endl;
        }
        if (iteration == 50000 && max_diff_iter > tolerance_sor * 10) {
            std::cout << "WARNING: Very slow convergence at 50k iterations. Error: " 
                      << std::scientific << max_diff_iter << " vs target: " << tolerance_sor << std::fixed << std::endl;
        }

        // Controllo di convergenza standard
        if (max_diff_iter < tolerance_sor) { // Usa tolerance_sor
            std::cout << "Converged after " << iteration + 1 << " iterations." << std::endl;
            break;
        }
        
        // Sistema di convergenza avanzato: controlla stagnazione nelle ultime 10 iterazioni
        int current_idx = iteration % convergence_window;
        double previous_diff = last_diffs[current_idx];
        last_diffs[current_idx] = max_diff_iter;
        
        // Dopo almeno 'convergence_window' iterazioni, controlla se c'è stagnazione
        if (iteration >= convergence_window) {
            // Calcola il miglioramento relativo rispetto a 10 iterazioni fa
            double relative_improvement = 0.0;
            if (previous_diff > 1e-12) { // Evita divisione per zero
                relative_improvement = (previous_diff - max_diff_iter) / previous_diff;
            }
            
            if (relative_improvement < relative_improvement_threshold) {
                no_improvement_count++;
                if (no_improvement_count >= convergence_window) {
                    std::cout << "Convergenza fermate dopo " << iteration + 1 << " iterazioni: ";
                    std::cout << "nessun miglioramento significativo nelle ultime " << convergence_window << " iterazioni." << std::endl;
                    std::cout << "Miglioramento relativo: " << std::scientific << relative_improvement << std::fixed;
                    std::cout << " (soglia: " << std::scientific << relative_improvement_threshold << std::fixed << ")" << std::endl;
                    break;
                }
            } else {
                no_improvement_count = 0; // Reset counter se c'è miglioramento
            }
        }
        
        if (iteration == max_iter_sor - 1) { // Usa max_iter_sor
             std::cout << "Max iterations (" << max_iter_sor << ") reached. Max diff: " << std::scientific << max_diff_iter << std::fixed << std::endl;
        }
    }

    // --- Calculate Electric Field with improved accuracy at material interfaces ---
    // Declare electric field vectors
    std::vector<std::vector<double>> Ex(Nx, std::vector<double>(Ny, 0.0));
    std::vector<std::vector<double>> Ey(Nx, std::vector<double>(Ny, 0.0));
    
    // Pre-compute the reciprocal of 2*h and h to avoid division in loops
    const double inv_2h = 1.0 / (2.0 * geom_config.h);
    const double inv_h = 1.0 / geom_config.h;
    
    std::cout << "Calculating electric field with interface-aware method..." << std::endl;
    
    // Calculate electric field considering material interfaces
    // Use 4th-order accurate scheme for Ey in uniform regions, 2nd-order at interfaces
    #pragma omp parallel for schedule(static) collapse(2)
    for (int i = 2; i < Nx - 2; ++i) {
        for (int j = 2; j < Ny - 2; ++j) {
            // Check if we're at a material interface for Ex calculation
            double eps_curr = eps_r[i][j];
            double eps_next = eps_r[i+1][j];
            double eps_prev = eps_r[i-1][j];
            
            // Detect material discontinuity (threshold 0.1 to avoid numerical noise)
            bool interface_x = (std::abs(eps_next - eps_curr) > 0.1 || std::abs(eps_prev - eps_curr) > 0.1);
            
            if (interface_x) {
                // At material interfaces, use weighted averaging for better accuracy
                // Weight by permittivity to account for field discontinuity
                double weight_next = eps_next / (eps_next + 1e-12);
                double weight_prev = eps_prev / (eps_prev + 1e-12);
                double weight_sum = weight_next + weight_prev;
                
                if (weight_sum > 0.01) {
                    // Weighted gradient considering material properties
                    Ex[i][j] = -(weight_next * (V[i+1][j] - V[i][j]) / geom_config.h + 
                                 weight_prev * (V[i][j] - V[i-1][j]) / geom_config.h) / weight_sum;
                } else {
                    // Fallback to standard central difference
                    Ex[i][j] = -(V[i+1][j] - V[i-1][j]) * inv_2h;
                }
            } else {
                // Standard central difference for uniform regions
                Ex[i][j] = -(V[i+1][j] - V[i-1][j]) * inv_2h;
            }
            
            // Enhanced treatment for Ey component with 4th-order accuracy
            double eps_up = eps_r[i][j+1];
            double eps_down = eps_r[i][j-1];
            double eps_up2 = eps_r[i][j+2];
            double eps_down2 = eps_r[i][j-2];
            
            bool interface_y = (std::abs(eps_up - eps_curr) > 0.1 || std::abs(eps_down - eps_curr) > 0.1);
            
            if (interface_y) {
                // At interfaces, use 2nd-order weighted scheme
                double weight_up = eps_up / (eps_up + 1e-12);
                double weight_down = eps_down / (eps_down + 1e-12);
                double weight_sum = weight_up + weight_down;
                
                if (weight_sum > 0.01) {
                    Ey[i][j] = -(weight_up * (V[i][j+1] - V[i][j]) / geom_config.h + 
                                 weight_down * (V[i][j] - V[i][j-1]) / geom_config.h) / weight_sum;
                } else {
                    Ey[i][j] = -(V[i][j+1] - V[i][j-1]) * inv_2h;
                }
            } else {
                // Check if 5-point stencil is also in uniform region
                bool uniform_5point = (std::abs(eps_up2 - eps_curr) < 0.1 && 
                                      std::abs(eps_down2 - eps_curr) < 0.1);
                
                if (uniform_5point) {
                    // 4th-order accurate central difference (5-point stencil)
                    // Formula: f'(x) ≈ [-f(x+2h) + 8f(x+h) - 8f(x-h) + f(x-2h)] / (12h)
                    Ey[i][j] = -(-V[i][j+2] + 8.0*V[i][j+1] - 8.0*V[i][j-1] + V[i][j-2]) / (12.0 * geom_config.h);
                } else {
                    // Fallback to 2nd-order if extended stencil crosses interface
                    Ey[i][j] = -(V[i][j+1] - V[i][j-1]) * inv_2h;
                }
            }
        }
    }
    
    // Handle edges and near-edge points with 2nd-order scheme
    #pragma omp parallel for
    for (int i = 1; i < Nx - 1; ++i) {
        for (int j = 1; j < Ny - 1; ++j) {
            // Skip points already computed with 4th-order scheme
            if (i >= 2 && i < Nx - 2 && j >= 2 && j < Ny - 2) continue;
            
            // Use 2nd-order scheme for near-boundary points
            Ex[i][j] = -(V[i+1][j] - V[i-1][j]) * inv_2h;
            Ey[i][j] = -(V[i][j+1] - V[i][j-1]) * inv_2h;
        }
    }
    
    // Handle edges with special care to avoid artifacts
    #pragma omp parallel for
    for (int j = 1; j < Ny - 1; ++j) {
        // Left boundary (i=0) - use forward difference
        Ex[0][j] = -(V[1][j] - V[0][j]) * inv_h;
        Ey[0][j] = -(V[0][j+1] - V[0][j-1]) * inv_2h;
        
        // Right boundary (i=Nx-1) - use backward difference
        Ex[Nx-1][j] = -(V[Nx-1][j] - V[Nx-2][j]) * inv_h;
        Ey[Nx-1][j] = -(V[Nx-1][j+1] - V[Nx-1][j-1]) * inv_2h;
    }
    
    // Top and bottom boundaries
    #pragma omp parallel for
    for (int i = 1; i < Nx - 1; ++i) {
        // Bottom boundary (j=0)
        Ex[i][0] = -(V[i+1][0] - V[i-1][0]) * inv_2h;
        Ey[i][0] = -(V[i][1] - V[i][0]) * inv_h;
        
        // Top boundary (j=Ny-1)
        Ex[i][Ny-1] = -(V[i+1][Ny-1] - V[i-1][Ny-1]) * inv_2h;
        Ey[i][Ny-1] = -(V[i][Ny-1] - V[i][Ny-2]) * inv_h;
    }
    
    // Corner points - use appropriate one-sided differences
    Ex[0][0] = -(V[1][0] - V[0][0]) * inv_h;
    Ey[0][0] = -(V[0][1] - V[0][0]) * inv_h;
    Ex[Nx-1][0] = -(V[Nx-1][0] - V[Nx-2][0]) * inv_h;
    Ey[Nx-1][0] = -(V[Nx-1][1] - V[Nx-1][0]) * inv_h;
    Ex[0][Ny-1] = -(V[1][Ny-1] - V[0][Ny-1]) * inv_h;
    Ey[0][Ny-1] = -(V[0][Ny-1] - V[0][Ny-2]) * inv_h;
    Ex[Nx-1][Ny-1] = -(V[Nx-1][Ny-1] - V[Nx-2][Ny-1]) * inv_h;
    Ey[Nx-1][Ny-1] = -(V[Nx-1][Ny-1] - V[Nx-1][Ny-2]) * inv_h;
    
    std::cout << "Electric field calculation completed." << std::endl;

    // --- Apply selective smoothing to Ey in vacuum regions ---
    std::cout << "Applying selective smoothing to Ey in vacuum regions..." << std::endl;
    
    // Create a temporary copy of Ey for smoothing
    std::vector<std::vector<double>> Ey_smoothed = Ey;
    
    // Gaussian kernel 3x3 (normalized)
    const double kernel[3][3] = {
        {0.0625, 0.125, 0.0625},
        {0.125,  0.25,  0.125},
        {0.0625, 0.125, 0.0625}
    };
    
    // Threshold to identify vacuum (eps_r close to 1.0)
    const double vacuum_threshold = 1.5;
    
    // Apply smoothing only in vacuum regions
    #pragma omp parallel for schedule(static) collapse(2)
    for (int i = 1; i < Nx - 1; ++i) {
        for (int j = 1; j < Ny - 1; ++j) {
            // Check if current point is in vacuum
            if (eps_r[i][j] < vacuum_threshold) {
                double sum = 0.0;
                double weight_sum = 0.0;
                
                // Apply 3x3 Gaussian kernel
                for (int di = -1; di <= 1; ++di) {
                    for (int dj = -1; dj <= 1; ++dj) {
                        int ni = i + di;
                        int nj = j + dj;
                        
                        // Only include neighboring vacuum points in smoothing
                        if (eps_r[ni][nj] < vacuum_threshold) {
                            double w = kernel[di + 1][dj + 1];
                            sum += Ey[ni][nj] * w;
                            weight_sum += w;
                        }
                    }
                }
                
                // Apply smoothed value if we had enough vacuum neighbors
                if (weight_sum > 0.1) {
                    Ey_smoothed[i][j] = sum / weight_sum;
                }
            }
        }
    }
    
    // Replace Ey with smoothed version
    Ey = Ey_smoothed;
    
    std::cout << "Selective smoothing completed." << std::endl;

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