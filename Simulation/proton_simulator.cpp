#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <random>
#include <algorithm> // For std::min, std::max
#include <sstream>   // For std::stringstream
#include <omp.h>     // For OpenMP
#include "geometry_definitions.h" // For addAluminumPlates

// Branch prediction hint macros
#if defined(__GNUC__) || defined(__clang__)
    #define LIKELY(x)    __builtin_expect(!!(x), 1)
    #define UNLIKELY(x)  __builtin_expect(!!(x), 0)
#else
    #define LIKELY(x)    (x)
    #define UNLIKELY(x)  (x)
#endif

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

// --- Constants ---
const double Q_ALPHA = 3.2043532e-19; // Coulombs (SI) - Elementary charge (2 * elementary charge)
const double M_ALPHA = 6.644657e-27;  // kg (SI) - Alpha particle mass
// K_ACCEL = Q_ALPHA / M_ALPHA, since E will be in V/m, acceleration will be m/s^2
const double K_ACCEL = Q_ALPHA / M_ALPHA; // SI units: (C/kg) * (V/m) -> m/s^2

// --- Simulation Parameters ---
const int NUM_ALPHAS = 10000;  // Reduced temporarily for testing
const double TOTAL_SIM_TIME_S = 1e-9;  // Total simulation time in seconds (SI)
// TIME_STEP_S and OUTPUT_TIME_INTERVAL_S will be calculated dynamically in main()
// Initial X position will be set in main after loading geometry, in meters

const double REL_PERMITTIVITY_MATERIAL_THRESHOLD = 1.1; // If eps_r > this, it's material (vacuum is ~1.0)

inline int find_cell_index(const std::vector<double>& coords, double pos);

// --- Helper Functions for Adaptive Time Step Calculation ---
double calculate_time_step(double e_k_keV, double h_grid_m, double cfl_factor = 0.01) {
    const double e_k_J = e_k_keV * 1000.0 * 1.60218e-19; // Convert keV to Joules
    const double v_max = std::sqrt(2 * e_k_J / M_ALPHA); // Maximum velocity (m/s)
    const double dt = cfl_factor * h_grid_m / v_max;
    return dt;
}

double calculate_output_interval(double e_k_keV, double L_device_m, int points_per_trajectory = 100) {
    const double e_k_J = e_k_keV * 1000.0 * 1.60218e-19;
    const double v_proton = std::sqrt(2 * e_k_J / M_ALPHA);
    const double T_transit = L_device_m / v_proton; // Time to traverse the device
    return T_transit / points_per_trajectory;
}

// --- Structures ---
struct Proton {
    int id;          // Unique identifier for the proton
    double x, y;     // Position in meters (SI)
    double vx, vy;   // Velocity in m/s (SI)
    bool active;
    // std::ofstream trajectory_file; // Removed: Will be handled locally
};

struct GeometryParameters {
    double h;         // Grid spacing in meters (SI) from geometry_params.csv
    double x_fs;      // x_free_space in meters (SI)
    double x_sl;      // x_structure_len in meters (SI)
    // Removed: y_vgt, H_tot_geom, y_si_base_height, initial_y_teeth_height
};

// --- Helper Functions ---
void create_directory_if_not_exists(const std::string& path) {
    struct STAT_STRUCT info;
    if (STAT_FUNC(path.c_str(), &info) != 0) {
        if (MKDIR(path.c_str()) == 0) {
            std::cout << "Output directory '" << path << "' created." << std::endl;
        } else {
            std::cerr << "Error: Could not create output directory '" << path << "'. Exiting." << std::endl;
            exit(1);
        }
    } else if (!(info.st_mode & S_IFDIR)) {
        std::cerr << "Error: '" << path << "' exists but is not a directory. Exiting." << std::endl;
        exit(1);
    } else {
        // Directory already exists, no message needed unless for debugging
        // std::cout << "Output directory '" << path << "' already exists." << std::endl;
    }
}

bool load_1d_csv(const std::string& filename, std::vector<double>& data_vec, bool convert_to_meters = false) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return false;
    }
    data_vec.clear();
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        try {
            double val = std::stod(line);
            if (convert_to_meters) {
                data_vec.push_back(val * 1.0e-6); // Convert µm to m
            } else {
                data_vec.push_back(val);
            }
        } catch (const std::invalid_argument& ia) {
            std::cerr << "Invalid argument: " << ia.what() << " in file " << filename << " line: " << line << std::endl;
            file.close();
            return false;
        } catch (const std::out_of_range& oor) {
            std::cerr << "Out of range: " << oor.what() << " in file " << filename << " line: " << line << std::endl;
            file.close();
            return false;
        }
    }
    file.close();
    return true;
}

bool load_field_csv(const std::string& filename, std::vector<std::vector<double>>& field_xy, int Nx, int Ny) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return false;
    }
    field_xy.assign(Nx, std::vector<double>(Ny));
    std::string line;
    int y_idx = 0;
    while (std::getline(file, line) && y_idx < Ny) {
        std::stringstream ss(line);
        std::string cell;
        int x_idx = 0;
        while (std::getline(ss, cell, ',') && x_idx < Nx) {
            try {
                field_xy[x_idx][y_idx] = std::stod(cell) * 1.0e6; // Convert V/µm to V/m
            } catch (const std::invalid_argument& ia) {
                std::cerr << "Invalid argument: " << ia.what() << " in file " << filename << " at (" << x_idx << "," << y_idx << ")" << " value: " << cell << std::endl;
                file.close();
                return false;
            } catch (const std::out_of_range& oor) {
                 std::cerr << "Out of range: " << oor.what() << " in file " << filename << " at (" << x_idx << "," << y_idx << ")" << " value: " << cell << std::endl;
                file.close();
                return false;
            }
            x_idx++;
        }
        if (x_idx != Nx) {
            std::cerr << "Error: Row " << y_idx << " in " << filename << " has " << x_idx << " columns, expected " << Nx << std::endl;
            // file.close(); return false; // Allow if it's just trailing commas on last line
        }
        y_idx++;
    }
    if (y_idx != Ny) {
        std::cerr << "Error: File " << filename << " has " << y_idx << " rows, expected " << Ny << std::endl;
        // file.close(); return false;
    }
    file.close();
    return true;
}

bool load_permittivity_map(const std::string& filename, std::vector<std::vector<double>>& eps_r_map, int Nx, int Ny) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open permittivity map file " << filename << std::endl;
        return false;
    }
    eps_r_map.assign(Nx, std::vector<double>(Ny));
    std::string line;
    int y_idx = 0;
    while (std::getline(file, line) && y_idx < Ny) {
        std::stringstream ss(line);
        std::string cell;
        int x_idx = 0;
        while (std::getline(ss, cell, ',') && x_idx < Nx) {
            try {
                eps_r_map[x_idx][y_idx] = std::stod(cell); // Permittivity values are unitless
            } catch (const std::invalid_argument& ia) {
                std::cerr << "Invalid argument in permittivity map: " << ia.what() << " in file " << filename << " at (" << x_idx << "," << y_idx << ")" << " value: " << cell << std::endl;
                file.close();
                return false;
            } catch (const std::out_of_range& oor) {
                 std::cerr << "Out of range in permittivity map: " << oor.what() << " in file " << filename << " at (" << x_idx << "," << y_idx << ")" << " value: " << cell << std::endl;
                file.close();
                return false;
            }
            x_idx++;
        }
        if (x_idx != Nx) {
            // std::cerr << "Warning: Row " << y_idx << " in permittivity map " << filename << " has " << x_idx << " columns, expected " << Nx << std::endl;
        }
        y_idx++;
    }
    if (y_idx != Ny) {
        // std::cerr << "Warning: Permittivity map file " << filename << " has " << y_idx << " rows, expected " << Ny << std::endl;
    }
    file.close();
    std::cout << "Permittivity map loaded successfully from " << filename << std::endl;
    return true;
}


bool load_geometry_params(const std::string& filename, GeometryParameters& geom) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open geometry parameters file " << filename << std::endl;
        return false;
    }
    std::string line;
    geom.h = 0.0; // Initialize to detect if not loaded
    geom.x_fs = 0.0;
    geom.x_sl = 0.0;
    bool h_loaded = false, x_fs_loaded = false, x_sl_loaded = false;

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string key, value_str;
        std::getline(ss, key, ',');
        std::getline(ss, value_str, ',');
        double value_um = 0.0;
        try {
            value_um = std::stod(value_str); 
        } catch (const std::exception& e) {
            std::cerr << "Error parsing value for key '" << key << "' in geometry_params: " << value_str << " (" << e.what() << ")" << std::endl;
            continue;
        }
        double value_m = value_um * 1.0e-6;     // Convert to meters

        if (key == "h") { geom.h = value_m; h_loaded = true; }
        else if (key == "x_free_space") { geom.x_fs = value_m; x_fs_loaded = true; }
        else if (key == "x_structure_len") { geom.x_sl = value_m; x_sl_loaded = true; }
    }
    file.close();
    
    if (!h_loaded) {
        std::cout << "Info: Grid spacing 'h' not found in " << filename << ". Will attempt to derive from coordinate files." << std::endl;
    }
    if (!x_fs_loaded || !x_sl_loaded) {
        std::cout << "Info: 'x_free_space' or 'x_structure_len' not found in " << filename << ". Will use fallbacks if needed for vacuum search." << std::endl;
    }
    
    std::cout << "Geometry parameters from CSV. h (m): " << geom.h 
              << ", x_fs (m): " << geom.x_fs 
              << ", x_sl (m): " << geom.x_sl << std::endl;
    return true;
}

// Optimized vacuum channel detection function with SIMD opportunities
std::pair<double, double> find_vacuum_channel_from_map(
    const std::vector<std::vector<double>>& eps_r_map,
    const std::vector<double>& x_coords_m, 
    const std::vector<double>& y_coords_m, 
    double h_grid_m_from_coords, 
    double x_search_start_m,     
    double x_search_end_m) {          // Fast fallback for edge cases
    const int Nx_map = eps_r_map.size();
    if (UNLIKELY(Nx_map == 0)) {
        double fallback_height = y_coords_m.back() - y_coords_m.front();
        return {y_coords_m.front() + fallback_height/3.0, y_coords_m.back() - fallback_height/3.0};
    }
    
    const int Ny_map = eps_r_map[0].size();
    if (UNLIKELY(Ny_map == 0 || y_coords_m.empty())) {
        return {0.0, 0.0};
    }

    // Find the exact index at the search position (should be a single column)
    // Use the middle of the search range to get the exact position
    double x_target_m = (x_search_start_m + x_search_end_m) / 2.0;
    int idx_x_target = find_cell_index(x_coords_m, x_target_m);
    
    // Clamp index to be within bounds
    idx_x_target = std::max(0, std::min(idx_x_target, Nx_map - 1));
    
    // Calculate gap at THIS specific x position only, no intersection with other columns
    int i = idx_x_target;
    
    // Scan from bottom to find the LAST (topmost) material cell
    int j_bottom_material_top_idx = -1;
    for (int j = 0; j < Ny_map; ++j) {
        if (eps_r_map[i][j] >= REL_PERMITTIVITY_MATERIAL_THRESHOLD) {
            j_bottom_material_top_idx = j; // Keep updating, we want the last one
        }
    }
    
    // Scan from top to find the FIRST (bottommost) material cell from the top
    int j_top_material_bottom_idx = -1;
    for (int j = Ny_map - 1; j >= 0; --j) {
        if (eps_r_map[i][j] >= REL_PERMITTIVITY_MATERIAL_THRESHOLD) {
            j_top_material_bottom_idx = j; // Keep updating, we want the first one from top
        }
    }
    
    // Calculate vacuum boundaries at this specific x position
    double overall_vac_min_y = (j_bottom_material_top_idx == -1) ? 
                              y_coords_m.front() : 
                              y_coords_m[j_bottom_material_top_idx] + h_grid_m_from_coords;
                              
    double overall_vac_max_y = (j_top_material_bottom_idx == -1) ? 
                           y_coords_m.back() : 
                           y_coords_m[j_top_material_bottom_idx];
    
    // Ensure bounds are respected
    overall_vac_min_y = std::max(overall_vac_min_y, y_coords_m.front());
    overall_vac_max_y = std::min(overall_vac_max_y, y_coords_m.back());
    
    bool found_any_gap_in_any_column = (overall_vac_min_y < overall_vac_max_y);

    if (!found_any_gap_in_any_column || overall_vac_min_y >= overall_vac_max_y) {
        // Fallback if no consistent channel found
        std::cerr << "Warning: Could not determine a consistent vacuum channel from epsilon map in the search range." << std::endl;
        double H_sim_total = y_coords_m.back() - y_coords_m.front();
        overall_vac_min_y = y_coords_m.front() + H_sim_total / 3.0; // Middle third
        overall_vac_max_y = y_coords_m.back() - H_sim_total / 3.0;
        std::cout << "Using fallback vacuum channel: [" << overall_vac_min_y * 1e6 << " um, " << overall_vac_max_y * 1e6 << " um]" << std::endl;
    } else {
        // Apply a small margin if a valid gap was found and is thick enough
        double gap_thickness = overall_vac_max_y - overall_vac_min_y;
        double min_sensible_gap_for_margin = 2.0 * h_grid_m_from_coords; // e.g., 2 grid cells
        if (gap_thickness > min_sensible_gap_for_margin) {
            double margin = h_grid_m_from_coords * 0.1; // 10% of grid spacing
            double new_min_y = overall_vac_min_y + margin;
            double new_max_y = overall_vac_max_y - margin;
            // Ensure margin doesn't invalidate a small but valid gap (must be at least half a cell wide after margin)
            if (new_min_y < new_max_y && (new_max_y - new_min_y) >= (h_grid_m_from_coords * 0.5) ) {
                overall_vac_min_y = new_min_y;
                overall_vac_max_y = new_max_y;
            }
        }
    }
    
    // The old margin logic is replaced by the clearer block above.
    // double gap_thickness = overall_vac_max_y - overall_vac_min_y;
    // if (gap_thickness > 2.0 * h_grid_m_from_coords) { 
    //     overall_vac_min_y += h_grid_m_from_coords * 0.1; 
    //     overall_vac_max_y -= h_grid_m_from_coords * 0.1;
    // }
    //  if (overall_vac_min_y >= overall_vac_max_y) { 
    //     overall_vac_min_y -= h_grid_m_from_coords * 0.1;
    //     overall_vac_max_y += h_grid_m_from_coords * 0.1;
    //     if (!found_any_gap_in_any_column) { 
    //          double H_sim_total = y_coords_m.back() - y_coords_m.front();
    //          overall_vac_min_y = y_coords_m.front() + H_sim_total / 3.0;
    //          overall_vac_max_y = y_coords_m.back() - H_sim_total / 3.0;
    //     }
    //  }

    return {overall_vac_min_y, overall_vac_max_y};
}

std::pair<double, double> get_field_at_point(
    double px, double py,
    const std::vector<double>& x_coords, const std::vector<double>& y_coords,
    const std::vector<std::vector<double>>& Ex_field, const std::vector<std::vector<double>>& Ey_field,
    int Nx, int Ny) {

    // Quick bounds check - avoids unnecessary calculations
    if (px < x_coords.front() || px > x_coords.back() || py < y_coords.front() || py > y_coords.back()) {
        return {0.0, 0.0}; // Outside grid
    }

    const int i = std::max(0, std::min(find_cell_index(x_coords, px), Nx - 2));
    const int j = std::max(0, std::min(find_cell_index(y_coords, py), Ny - 2));

    const double x0 = x_coords[i];
    const double x1 = x_coords[i + 1];
    const double y0 = y_coords[j];
    const double y1 = y_coords[j + 1];

    const double tx = (x1 - x0) > 0.0 ? (px - x0) / (x1 - x0) : 0.0;
    const double ty = (y1 - y0) > 0.0 ? (py - y0) / (y1 - y0) : 0.0;
    const double tx_comp = 1.0 - tx;
    const double ty_comp = 1.0 - ty;

    // Cache field values to improve memory access patterns
    const double ex00 = Ex_field[i][j];
    const double ex10 = Ex_field[i+1][j];
    const double ex01 = Ex_field[i][j+1];
    const double ex11 = Ex_field[i+1][j+1];
    
    const double ey00 = Ey_field[i][j];
    const double ey10 = Ey_field[i+1][j];
    const double ey01 = Ey_field[i][j+1];
    const double ey11 = Ey_field[i+1][j+1];

    // Bilinear interpolation combined for Ex and Ey (reduces computation)
    const double tx_ty_comp = tx_comp * ty_comp;
    const double tx_ty = tx * ty;
    const double tx_ty_comp_y = tx_comp * ty;
    const double tx_ty_comp_x = tx * ty_comp;
    
    const double Ex_interp = tx_ty_comp * ex00 + tx_ty_comp_x * ex10 + 
                            tx_ty_comp_y * ex01 + tx_ty * ex11;
    
    const double Ey_interp = tx_ty_comp * ey00 + tx_ty_comp_x * ey10 + 
                            tx_ty_comp_y * ey01 + tx_ty * ey11;

    return {Ex_interp, Ey_interp};
}

// Optimized acceleration calculation with inlining hint
inline std::pair<double, double> get_acceleration(
    double px, double py,
    const std::vector<double>& x_coords, const std::vector<double>& y_coords,
    const std::vector<std::vector<double>>& Ex_field, const std::vector<std::vector<double>>& Ey_field,
    int Nx, int Ny) {

    std::pair<double, double> E_vec = get_field_at_point(px, py, x_coords, y_coords, Ex_field, Ey_field, Nx, Ny);
    
    // Using the constant directly for multiplication saves a memory lookup
    return {K_ACCEL * E_vec.first, K_ACCEL * E_vec.second}; // (ax, ay) in m/s^2
}

// Optimized with inline hint for better compiler optimization
inline int find_cell_index(const std::vector<double>& coords, double pos) {
    if (pos <= coords.front()) return 0;
    if (pos >= coords.back()) return static_cast<int>(coords.size()) - 2;
    auto upper = std::upper_bound(coords.begin(), coords.end(), pos);
    int idx = static_cast<int>(std::distance(coords.begin(), upper)) - 1;
    return std::max(0, idx);
}

inline bool is_in_material_or_out_of_bounds(
    double px, double py,
    const std::vector<double>& x_coords_m,
    const std::vector<double>& y_coords_m,
    const std::vector<std::vector<double>>& eps_r_map,
    double L_total_sim_m, double H_total_sim_m) {

    if (px <= 0.0 || px >= L_total_sim_m || py <= 0.0 || py >= H_total_sim_m) {
        return true;
    }

    const int Nx_map = static_cast<int>(eps_r_map.size());
    if (UNLIKELY(Nx_map == 0)) return true;
    const int Ny_map = static_cast<int>(eps_r_map[0].size());
    if (UNLIKELY(Ny_map == 0)) return true;

    const int i = std::max(0, std::min(find_cell_index(x_coords_m, px), Nx_map - 1));
    const int j = std::max(0, std::min(find_cell_index(y_coords_m, py), Ny_map - 1));

    return eps_r_map[i][j] >= REL_PERMITTIVITY_MATERIAL_THRESHOLD;
}


int main(int argc, char* argv[]) { // Modified main signature
    
    std::string input_base_folder_name = "geometria_Denti_sfasati_profondi_5um_default"; // Default input folder
    double user_y_min_um = -1.0; // -1 means not specified, will use automatic detection
    double user_y_max_um = -1.0; // -1 means not specified, will use automatic detection
    
    if (argc > 1) {
        input_base_folder_name = argv[1]; // Use the first command-line argument as folder name
        std::cout << "Input data folder specified: " << input_base_folder_name << std::endl;
    } else {
        std::cout << "No input data folder specified, using default: " << input_base_folder_name << std::endl;
    }
    
    // ============================================
    // SIMULATION PARAMETERS (DEFINE ONCE)
    // ============================================
    const double INITIAL_ALPHA_ENERGY_KEV = 50.0; // Initial alpha energy in keV
    
    // Parse optional y_min and y_max parameters (in micrometers)
    // Usage: ./proton_simulator <folder> [y_min_um] [y_max_um]
    if (argc > 2) {
        try {
            user_y_min_um = std::stod(argv[2]);
            std::cout << "User-specified y_min: " << user_y_min_um << " um" << std::endl;
        } catch (...) {
            std::cerr << "Warning: Invalid y_min argument, will use automatic detection." << std::endl;
            user_y_min_um = -1.0;
        }
    }
    
    if (argc > 3) {
        try {
            user_y_max_um = std::stod(argv[3]);
            std::cout << "User-specified y_max: " << user_y_max_um << " um" << std::endl;
        } catch (...) {
            std::cerr << "Warning: Invalid y_max argument, will use automatic detection." << std::endl;
            user_y_max_um = -1.0;
        }
    }
    
    const std::string folder = input_base_folder_name; // Use 'folder' as it was used before, now initialized from CLI or default

    // The line below was: const std::string input_base_folder = folder; 
    // It's redundant now as 'folder' serves this purpose.
    // We will use 'folder' directly where 'input_base_folder' was used.

    // const std::string output_traj_folder = folder + "/proton_trajectories"; // Removed
    // create_directory_if_not_exists(output_traj_folder); // Removed

    const std::string all_trajectories_filename = folder + "/all_proton_trajectories.csv";
    std::ofstream all_trajectories_file_stream;
    
    // Optimize file I/O with larger buffer and binary mode
    const int OUTPUT_BUFFER_SIZE = 1024 * 1024; // 1MB buffer for better I/O performance
    char* output_buffer = new char[OUTPUT_BUFFER_SIZE];
    
    all_trajectories_file_stream.open(all_trajectories_filename, std::ios::out | std::ios::binary);
    all_trajectories_file_stream.rdbuf()->pubsetbuf(output_buffer, OUTPUT_BUFFER_SIZE);
    
    if (!all_trajectories_file_stream.is_open()) {
        std::cerr << "Error: Could not open the consolidated trajectory file: " << all_trajectories_filename << ". Exiting." << std::endl;
        delete[] output_buffer;
        return 1;
    }
    
    // Set numeric precision once for the entire stream
    all_trajectories_file_stream << std::scientific << std::setprecision(12);
    all_trajectories_file_stream << "proton_id,time_s,x_m,y_m,vx_m_per_s,vy_m_per_s\n";


    GeometryParameters geom;
    if (!load_geometry_params(folder + "/geometry_params.csv", geom)) { // Use 'folder'
        std::cerr << "Warning: Problem loading from geometry_params.csv. Proceeding with defaults/derivations." << std::endl;
        // geom.h, geom.x_fs, geom.x_sl will be 0.0 if file not found or keys missing.
    }

    std::vector<double> x_coords, y_coords;
    if (!load_1d_csv(folder + "/x_coordinates.csv", x_coords, true /*convert_to_meters*/) || // Use 'folder'
        !load_1d_csv(folder + "/y_coordinates.csv", y_coords, true /*convert_to_meters*/)) { // Use 'folder'
        std::cerr << "Failed to load coordinates. Exiting." << std::endl;
        return 1;
    }

    int Nx = x_coords.size();
    int Ny = y_coords.size();
    if (Nx == 0 || Ny == 0) {
        std::cerr << "Coordinate data is empty. Exiting." << std::endl;
        return 1;
    }
    double L_total_sim = x_coords.back(); 
    double H_total_sim = y_coords.back(); 

    double h_from_coords_x = (Nx > 1) ? (x_coords[1] - x_coords[0]) : 0.0;
    double h_from_coords_y = (Ny > 1) ? (y_coords[1] - y_coords[0]) : 0.0;
    double h_for_simulation = geom.h; 

    if (geom.h <= 1e-9) { // geom.h not loaded or zero
        if (h_from_coords_x > 1e-9 && h_from_coords_y > 1e-9 && std::abs(h_from_coords_x - h_from_coords_y) < 1e-9) { // Use 1e-9 as effective zero for double
            h_for_simulation = h_from_coords_x;
            std::cout << "Using h derived from coordinate files: " << h_for_simulation * 1e6 << " um." << std::endl;
        } else if (h_from_coords_x > 1e-9) { // If only x is reliable
             h_for_simulation = h_from_coords_x;
             std::cout << "Warning: y-coordinate spacing inconsistent or zero. Using h from x-coordinates: " << h_for_simulation * 1e6 << " um." << std::endl;
        } else if (h_from_coords_y > 1e-9) { // If only y is reliable
             h_for_simulation = h_from_coords_y;
             std::cout << "Warning: x-coordinate spacing inconsistent or zero. Using h from y-coordinates: " << h_for_simulation * 1e6 << " um." << std::endl;
        }
         else {
            std::cerr << "Error: h from geometry_params.csv is invalid, and h cannot be reliably derived from coordinate files. Exiting." << std::endl;
            return 1;
        }
    } else { 
         if ( (h_from_coords_x > 1e-9 && std::abs(geom.h - h_from_coords_x) > 1e-9) || 
              (h_from_coords_y > 1e-9 && std::abs(geom.h - h_from_coords_y) > 1e-9) ) {
            std::cout << "Warning: h from geometry_params.csv (" << geom.h * 1e6 
                      << " um) differs from h derived from coordinates (x: " << h_from_coords_x * 1e6 
                      << " um, y: " << h_from_coords_y * 1e6 
                      << " um). Using h from geometry_params.csv: " << h_for_simulation * 1e6 << " um." << std::endl;
        }
    }
     if (h_for_simulation <= 1e-9) {
        std::cerr << "Error: Final h_for_simulation is invalid (zero or negative). Exiting." << std::endl;
        return 1;
     }

    // --- Calculate Adaptive Time Step Parameters ---
    // Physical parameters for time step calculation
    const double CFL_FACTOR = 0.01; // Safety factor for CFL condition (0.01 = 1% for high precision)
    const int POINTS_PER_TRAJECTORY = 5000; // Number of output points along the entire trajectory
    
    // Calculate adaptive time step based on CFL condition
    const double TIME_STEP_S = calculate_time_step(INITIAL_ALPHA_ENERGY_KEV, h_for_simulation, CFL_FACTOR);
    
    // Calculate output interval based on device transit time
    const double OUTPUT_TIME_INTERVAL_S = calculate_output_interval(INITIAL_ALPHA_ENERGY_KEV, L_total_sim, POINTS_PER_TRAJECTORY);
    
    // Calculate initial proton velocity for verification
    const double e_k_initial_J = INITIAL_ALPHA_ENERGY_KEV * 1000.0 * 1.60218e-19;
    const double v_initial = std::sqrt(2 * e_k_initial_J / M_ALPHA);
    const double T_transit = L_total_sim / v_initial; // Transit time through device
    
    // Print calculated parameters
    std::cout << "\n=== ADAPTIVE TIME STEP PARAMETERS ===" << std::endl;
    std::cout << "Initial alpha energy: " << INITIAL_ALPHA_ENERGY_KEV << " keV" << std::endl;
    std::cout << "Initial alpha velocity: " << v_initial / 1e6 << " × 10^6 m/s" << std::endl;
    std::cout << "Device length: " << L_total_sim * 1e6 << " µm (" << L_total_sim * 1e3 << " mm)" << std::endl;
    std::cout << "Grid spacing (h): " << h_for_simulation * 1e6 << " µm" << std::endl;
    std::cout << "CFL factor: " << CFL_FACTOR << std::endl;
    std::cout << "Calculated TIME_STEP: " << TIME_STEP_S << " s (" << TIME_STEP_S * 1e15 << " fs)" << std::endl;
    std::cout << "Expected transit time: " << T_transit * 1e9 << " ns" << std::endl;
    std::cout << "OUTPUT_TIME_INTERVAL: " << OUTPUT_TIME_INTERVAL_S << " s (" << OUTPUT_TIME_INTERVAL_S * 1e12 << " ps)" << std::endl;
    std::cout << "Total simulation time: " << TOTAL_SIM_TIME_S << " s (" << TOTAL_SIM_TIME_S * 1e9 << " ns)" << std::endl;
    std::cout << "Target output points per trajectory: " << POINTS_PER_TRAJECTORY << std::endl;
    std::cout << "Estimated total steps: " << static_cast<int>(TOTAL_SIM_TIME_S / TIME_STEP_S) << std::endl;
    std::cout << "Estimated output points per proton: " << static_cast<int>(TOTAL_SIM_TIME_S / OUTPUT_TIME_INTERVAL_S) << std::endl;
    std::cout << "======================================\n" << std::endl;

    double initial_x_position_m = 5e-6;
    
    std::vector<std::vector<double>> Ex_field, Ey_field; 
    std::vector<std::vector<double>> eps_r_map_data;     

    if (!load_field_csv(folder + "/electric_field_x.csv", Ex_field, Nx, Ny) || // Use 'folder'
        !load_field_csv(folder + "/electric_field_y.csv", Ey_field, Nx, Ny) || // Use 'folder'
        !load_permittivity_map(folder + "/permittivity.csv", eps_r_map_data, Nx, Ny)) { // Use 'folder'
        std::cerr << "Failed to load electric field or permittivity data. Exiting." << std::endl;
        return 1;
    }
    std::cout << "Data loaded successfully. Nx=" << Nx << ", Ny=" << Ny << std::endl;
    
    // Apply aluminum plates to permittivity map (if needed, for proton collision detection)
    // Note: The field data is pre-computed with aluminum plates already considered
    // This is just to ensure the eps_r_map used for material detection includes aluminum
    GeometryConfig dummy_config;
    dummy_config.h = h_for_simulation;
    dummy_config.x_free_space = geom.x_fs;
    dummy_config.x_structure_len = geom.x_sl;
    dummy_config.aluminum_plate_thickness = 1.0e-6; // 1 µm
    dummy_config.eps_aluminum = 1000.0;
    
    // Only apply if we have valid geometry parameters
    if (geom.x_fs > 0 && geom.x_sl > 0) {
        std::cout << "\n=== Applying Aluminum Plates to Loaded Permittivity Map ===" << std::endl;
        addAluminumPlates(eps_r_map_data, dummy_config, Nx, Ny);
    }
    
    // Calculate the vacuum gap specifically at the initial proton position
    // This ensures protons start at the correct gap height at the beginning of the structure
    double x_search_for_gap = initial_x_position_m;
    double x_search_start_for_gap = x_search_for_gap - 0.5 * h_for_simulation; // Small window around initial position
    double x_search_end_for_gap = x_search_for_gap + 0.5 * h_for_simulation;
    
    // Ensure search range is within bounds
    x_search_start_for_gap = std::max(x_search_start_for_gap, 0.0);
    x_search_end_for_gap = std::min(x_search_end_for_gap, L_total_sim);
    
    if (x_search_start_for_gap >= x_search_end_for_gap) {
        std::cerr << "Error: Invalid search range for vacuum gap calculation. Exiting." << std::endl;
        return 1;
    }
    
    std::cout << "Calculating vacuum channel gap at initial proton position x = " 
              << initial_x_position_m * 1e6 << " um (search range: [" 
              << x_search_start_for_gap * 1e6 << ", " << x_search_end_for_gap * 1e6 << "] um)." << std::endl;
    
    double h_coords_for_vacuum_find = (h_from_coords_y > 1e-9) ? h_from_coords_y : h_from_coords_x;
    if (h_coords_for_vacuum_find <= 1e-9) h_coords_for_vacuum_find = h_for_simulation; // Last resort

    double vacuum_gap_start_y, vacuum_gap_end_y;
    
    // Check if user specified y range manually
    if (user_y_min_um > 0 && user_y_max_um > 0) {
        // Use user-specified range
        vacuum_gap_start_y = user_y_min_um * 1e-6; // Convert um to m
        vacuum_gap_end_y = user_y_max_um * 1e-6;   // Convert um to m
        std::cout << "Using user-specified y-range: [" << user_y_min_um << ", " << user_y_max_um << "] um" << std::endl;
    } else {
        // Automatic detection from permittivity map at the beginning of the structure
        std::pair<double, double> vacuum_channel = find_vacuum_channel_from_map(
            eps_r_map_data, x_coords, y_coords, h_coords_for_vacuum_find,
            x_search_start_for_gap, x_search_end_for_gap
        );
        double gap_center_y = (vacuum_channel.first + vacuum_channel.second) / 2.0;
        double emission_width_m = 10.0e-6; // 10 micrometers total width
        
        vacuum_gap_start_y = gap_center_y - emission_width_m / 2.0;
        vacuum_gap_end_y = gap_center_y + emission_width_m / 2.0;
        
        std::cout << "Auto-detected gap center at y = " << gap_center_y*1e6 << " um" << std::endl;
        std::cout << "Emission window: ±5 um around center = [" 
                  << vacuum_gap_start_y*1e6 << ", " << vacuum_gap_end_y*1e6 << "] um" << std::endl;
    }

    if (vacuum_gap_start_y >= vacuum_gap_end_y - (h_for_simulation * 0.5) ) { // Ensure gap is at least half a cell thick
        std::cerr << "Critical Error: Determined vacuum gap is too small or invalid [" << vacuum_gap_start_y << ", " << vacuum_gap_end_y << "]. Exiting." << std::endl;
        return 1;
    }
    
    std::vector<Proton> protons(NUM_ALPHAS);

    // Allow user to specify RNG seed via command line (fourth argument), else use random_device
    // Usage: ./proton_simulator <folder> [y_min_um] [y_max_um] [rng_seed]
    unsigned int rng_seed = std::random_device{}();
    if (argc > 4) {
        try {
            rng_seed = static_cast<unsigned int>(std::stoul(argv[4]));
            std::cout << "Using user-specified RNG seed: " << rng_seed << std::endl;
        } catch (...) {
            std::cerr << "Warning: Invalid RNG seed argument, using random_device." << std::endl;
        }
    } else {
        std::cout << "No RNG seed specified, using random_device: " << rng_seed << std::endl;
    }
    std::mt19937 rng(rng_seed);

    // Calculate emission range
    double y_emit_min, y_emit_max;
    
    if (user_y_min_um > 0 && user_y_max_um > 0) {
        // User specified values - use them directly without margins
        y_emit_min = vacuum_gap_start_y;
        y_emit_max = vacuum_gap_end_y;
        std::cout << "Using user-specified y-emission range directly: [" 
                  << y_emit_min*1e6 << ", " << y_emit_max*1e6 << "] um" << std::endl;
    } else {
        // Auto mode: already calculated as ±5 µm around center
        y_emit_min = vacuum_gap_start_y;
        y_emit_max = vacuum_gap_end_y;
    }
    
    // Calculate initial proton velocity from the global energy constant
    const double e_k_J = INITIAL_ALPHA_ENERGY_KEV * 1000.0 * 1.60218e-19; // Convert keV -> eV -> Joules
    const double v_total = std::sqrt(2 * e_k_J / M_ALPHA); // Initial velocity in m/s
    
    std::cout << "\n=== FINAL ALPHA EMISSION PARAMETERS ===" << std::endl;
    std::cout << "Initial X position: " << initial_x_position_m * 1e6 << " um" << std::endl;
    std::cout << "Initial Y range: [" << y_emit_min * 1e6 << ", " << y_emit_max * 1e6 << "] um" << std::endl;
    std::cout << "Y range width: " << (y_emit_max - y_emit_min) * 1e6 << " um" << std::endl;
    std::cout << "Initial velocity: " << v_total << " m/s (" << INITIAL_ALPHA_ENERGY_KEV << " keV)" << std::endl;
    std::cout << "Angular spread: +/- 5 degrees (10 degree total range)" << std::endl;
    std::cout << "========================================\n" << std::endl;
    
    std::uniform_real_distribution<double> dist_y(y_emit_min, y_emit_max);
    std::uniform_real_distribution<double> dist_angle(-2.0 * M_PI / 180.0, 2.0 * M_PI / 180.0); // +/- 5 degrees in radians

    for (int i = 0; i < NUM_ALPHAS; ++i) {
        protons[i].id = i; // Assign ID
        protons[i].x = initial_x_position_m; 
        protons[i].y = dist_y(rng); // y in meters

        // Random initial direction within +/- 40 degrees
        double angle = dist_angle(rng);
        protons[i].vx = v_total * std::cos(angle); // x-component of velocity
        protons[i].vy = v_total * std::sin(angle); // y-component of velocity
        protons[i].active = true;

        // Write initial state to the single file
        all_trajectories_file_stream << protons[i].id << "," << 0.0 << "," << protons[i].x << "," << protons[i].y << "," << protons[i].vx << "," << protons[i].vy << "\n";
    }

    double current_time = 0.0;
    int num_steps = static_cast<int>(TOTAL_SIM_TIME_S / TIME_STEP_S);
    int output_every_n_steps = static_cast<int>(std::max(1.0, OUTPUT_TIME_INTERVAL_S / TIME_STEP_S));
    int active_protons_count = NUM_ALPHAS;
    int protons_reached_end_successfully = 0; // Counter for successful protons

    std::cout << "Starting simulation for " << NUM_ALPHAS << " alphas..." << std::endl;
    std::cout << "Total steps: " << num_steps << ", Outputting every " << output_every_n_steps << " steps." << std::endl;

    // Main simulation loop with optimized iteration strategy
    for (int step = 0; step < num_steps; ++step) {
        // Early exit if no active protons
        if (UNLIKELY(active_protons_count == 0)) {
            std::cout << "All protons inactive. Stopping simulation early at step " << step << "." << std::endl;
            break;
        }

        // Calculate output needed only once per step
        const bool needs_output_this_step = ((step + 1) % output_every_n_steps == 0);
        const double step_time = current_time + TIME_STEP_S;

        // Parallel region with optimized scheduling for better load balancing
        #pragma omp parallel for schedule(dynamic, 32)
        for (int i = 0; i < NUM_ALPHAS; ++i) {
            // Skip inactive protons immediately
            if (!protons[i].active) continue;

            // Reference to current proton for better readability and avoids extra copying
            Proton& p = protons[i];

            // Local cache of acceleration and RK4 variables
            double px = p.x, py = p.y;
            double vx = p.vx, vy = p.vy;
            double ax1, ay1, ax2, ay2, ax3, ay3, ax4, ay4;
            double k1x, k1y, k1vx, k1vy;
            double k2x, k2y, k2vx, k2vy;
            double k3x, k3y, k3vx, k3vy;
            double k4x, k4y, k4vx, k4vy;

            // k1 computation - original position and velocity
            std::pair<double, double> a1 = get_acceleration(px, py, x_coords, y_coords, Ex_field, Ey_field, Nx, Ny);
            ax1 = a1.first;
            ay1 = a1.second;
            k1vx = TIME_STEP_S * ax1;
            k1vy = TIME_STEP_S * ay1;
            k1x = TIME_STEP_S * vx;
            k1y = TIME_STEP_S * vy;

            // k2 computation - midpoint using k1
            std::pair<double, double> a2 = get_acceleration(px + k1x*0.5, py + k1y*0.5, x_coords, y_coords, Ex_field, Ey_field, Nx, Ny);
            ax2 = a2.first;
            ay2 = a2.second;
            k2vx = TIME_STEP_S * ax2;
            k2vy = TIME_STEP_S * ay2;
            k2x = TIME_STEP_S * (vx + k1vx*0.5);
            k2y = TIME_STEP_S * (vy + k1vy*0.5);

            // k3 computation - midpoint using k2
            std::pair<double, double> a3 = get_acceleration(px + k2x*0.5, py + k2y*0.5, x_coords, y_coords, Ex_field, Ey_field, Nx, Ny);
            ax3 = a3.first;
            ay3 = a3.second;
            k3vx = TIME_STEP_S * ax3;
            k3vy = TIME_STEP_S * ay3;
            k3x = TIME_STEP_S * (vx + k2vx*0.5);
            k3y = TIME_STEP_S * (vy + k2vy*0.5);

            // k4 computation - endpoint using k3
            std::pair<double, double> a4 = get_acceleration(px + k3x, py + k3y, x_coords, y_coords, Ex_field, Ey_field, Nx, Ny);
            ax4 = a4.first;
            ay4 = a4.second;
            k4vx = TIME_STEP_S * ax4;
            k4vy = TIME_STEP_S * ay4;
            k4x = TIME_STEP_S * (vx + k3vx);
            k4y = TIME_STEP_S * (vy + k3vy);

            // Optimized RK4 update with constant multipliers pre-calculated
            const double one_sixth = 1.0/6.0;
            const double one_third = 1.0/3.0;
            px += (k1x + k4x) * one_sixth + (k2x + k3x) * one_third;
            py += (k1y + k4y) * one_sixth + (k2y + k3y) * one_third;
            vx += (k1vx + k4vx) * one_sixth + (k2vx + k3vx) * one_third;
            vy += (k1vy + k4vy) * one_sixth + (k2vy + k3vy) * one_third;

            // Update proton state with optimized values
            p.x = px;
            p.y = py;
            p.vx = vx;
            p.vy = vy;


            // Check if proton reached the end successfully - optimized branch prediction
            if (UNLIKELY(px >= L_total_sim)) {
                if (p.active) {
                    p.active = false;
                    
                    // Atomic operations for thread-safety, but minimized for better performance
                    #pragma omp atomic update
                    protons_reached_end_successfully++;
                    
                    #pragma omp atomic update
                    active_protons_count--;
                    
                    // Use string buffer for better I/O performance
                    #pragma omp critical (trajectory_file_write)
                    {
                        all_trajectories_file_stream 
                            << p.id << ',' 
                            << step_time << ',' 
                            << px << ',' 
                            << py << ',' 
                            << vx << ',' 
                            << vy << '\n';
                    }
                }
                continue; 
            }

            // Check for material collision or out of bounds - less frequent condition
            if (p.active && is_in_material_or_out_of_bounds(px, py, x_coords, y_coords, eps_r_map_data, L_total_sim, H_total_sim)) {
                p.active = false;
                
                #pragma omp atomic update
                active_protons_count--;
                
                #pragma omp critical (trajectory_file_write)
                {
                    all_trajectories_file_stream 
                        << p.id << ',' 
                        << step_time << ',' 
                        << px << ',' 
                        << py << ',' 
                        << vx << ',' 
                        << vy << '\n';
                }
            }

            // Regular output if still active and it's an output step - less frequent condition
            else if (p.active && needs_output_this_step) {
                #pragma omp critical (trajectory_file_write)
                {
                    all_trajectories_file_stream 
                        << p.id << ',' 
                        << step_time << ',' 
                        << px << ',' 
                        << py << ',' 
                        << vx << ',' 
                        << vy << '\n';
                }
            }
        } // End of parallel for loop over protons

        // Update time for next step
        current_time += TIME_STEP_S;
        
        // Progress reporting with optimized modulo calculation and flush only when needed
        if ((step + 1) % std::max(1, num_steps / 100) == 0 || step == num_steps - 1) { 
            double progress = static_cast<double>(step + 1) / num_steps * 100.0;
            std::cout << "\rSimulation progress: " << std::fixed << std::setprecision(2) 
                     << progress << "% (" << active_protons_count << " active protons)" << std::flush;
        }
    }
    std::cout << "\nSimulation finished." << std::endl;

    // Ensure all data is flushed and file is properly closed
    all_trajectories_file_stream.flush();
    all_trajectories_file_stream.close();
    delete[] output_buffer; // Clean up the buffer

    double success_percentage = 0.0;
    if (NUM_ALPHAS > 0) {
        success_percentage = (static_cast<double>(protons_reached_end_successfully) / NUM_ALPHAS) * 100.0;
    }
    std::cout << protons_reached_end_successfully << " out of " << NUM_ALPHAS << " alphas reached the end successfully." << std::endl;
    std::cout << "Success percentage: " << std::fixed << std::setprecision(2) << success_percentage << "%" << std::endl;

    std::cout << "All proton trajectories saved to '" << all_trajectories_filename << "' (data in SI units)." << std::endl;
    return 0;
}

