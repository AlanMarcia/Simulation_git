#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <random>
#include <algorithm> // For std::min, std::max
#include <sstream>   // For std::stringstream

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
const double Q_PROTON = 1.60217663e-19; // Coulombs (SI)
const double M_PROTON = 1.67262192e-27; // kg (SI)
// K_ACCEL = Q_PROTON / M_PROTON, since E will be in V/m, acceleration will be m/s^2
const double K_ACCEL = Q_PROTON / M_PROTON; // SI units: (C/kg) * (V/m) -> m/s^2

// --- Simulation Parameters ---
const int NUM_PROTONS = 10000;
const double TIME_STEP_S = 1e-12;       // Time step in seconds (SI)
const double TOTAL_SIM_TIME_S = 1e-8;  // Total simulation time in seconds (SI)
const double OUTPUT_TIME_INTERVAL_S = 1e-10; // Interval for writing trajectory data (SI)
// Initial X position will be set in main after loading geometry, in meters

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

std::pair<double, double> find_vacuum_channel_from_map(
    const std::vector<std::vector<double>>& eps_r_map,
    const std::vector<double>& x_coords_m, 
    const std::vector<double>& y_coords_m, 
    double h_grid_m_from_coords, 
    double x_search_start_m,     
    double x_search_end_m) {      

    const double EPS_SI_SIM = 11.7;
    const double EPS_VAC_SIM = 1.0;
    const double material_threshold = (EPS_SI_SIM + EPS_VAC_SIM) / 2.0;

    int Nx_map = eps_r_map.size();
    if (Nx_map == 0) return {y_coords_m.front() + (y_coords_m.back()-y_coords_m.front())/3.0, y_coords_m.back() - (y_coords_m.back()-y_coords_m.front())/3.0};
    int Ny_map = eps_r_map[0].size();
    if (Ny_map == 0 || y_coords_m.empty()) return {0.0, 0.0}; // Should not happen with earlier checks

    double x_coord_spacing = (x_coords_m.size() > 1) ? (x_coords_m[1] - x_coords_m[0]) : 1.0;
    if (x_coord_spacing <= 0) x_coord_spacing = 1.0; // Avoid division by zero

    int idx_x_start_search = static_cast<int>(std::max(0.0, x_search_start_m / x_coord_spacing ));
    int idx_x_end_search = static_cast<int>(std::min((double)Nx_map - 1, x_search_end_m / x_coord_spacing ));
    
    idx_x_start_search = std::min(idx_x_start_search, Nx_map -1);
    idx_x_end_search = std::min(idx_x_end_search, Nx_map-1);

    if (idx_x_start_search < 0) idx_x_start_search = 0;
    if (idx_x_end_search < 0) idx_x_end_search = Nx_map -1;


    if (idx_x_start_search > idx_x_end_search || Nx_map == 0) { 
        std::cout << "Warning: Invalid x-search range for vacuum. Analyzing central 50% of x-domain." << std::endl;
        idx_x_start_search = std::max(0, Nx_map / 4); 
        idx_x_end_search = std::min(Nx_map - 1, Nx_map * 3 / 4);
        if (idx_x_start_search > idx_x_end_search && Nx_map > 0) idx_x_end_search = idx_x_start_search;
    }
    if (idx_x_start_search > idx_x_end_search && Nx_map == 0) {idx_x_start_search = 0; idx_x_end_search = -1;}


    double overall_vac_min_y = y_coords_m.back(); 
    double overall_vac_max_y = y_coords_m.front();
    bool found_any_gap_in_any_column = false;

    for (int i = idx_x_start_search; i <= idx_x_end_search; ++i) {
        int j_bottom_material_top_idx = -1; // Index of the highest cell of bottom material
        for (int j = 0; j < Ny_map; ++j) {
            if (eps_r_map[i][j] >= material_threshold) {
                j_bottom_material_top_idx = j;
            } else {
                break; // First vacuum cell found, so material below ends at j-1 or this is start of domain
            }
        }

        int j_top_material_bottom_idx = -1; // Index of the lowest cell of top material
        for (int j = Ny_map - 1; j >= 0; --j) {
            if (eps_r_map[i][j] >= material_threshold) {
                j_top_material_bottom_idx = j;
            } else {
                break; // First vacuum cell found from top
            }
        }
        
        double col_vac_start_y, col_vac_end_y;

        if (j_bottom_material_top_idx == -1) { // All vacuum from bottom
            col_vac_start_y = y_coords_m.front();
        } else {
            col_vac_start_y = y_coords_m[j_bottom_material_top_idx] + h_grid_m_from_coords; // Start of vacuum is top of this material cell
        }

        if (j_top_material_bottom_idx == -1) { // All vacuum to top
            col_vac_end_y = y_coords_m.back();
        } else {
            col_vac_end_y = y_coords_m[j_top_material_bottom_idx]; // End of vacuum is bottom of this material cell
        }
        col_vac_start_y = std::max(col_vac_start_y, y_coords_m.front());
        col_vac_end_y = std::min(col_vac_end_y, y_coords_m.back());


        if (col_vac_start_y < col_vac_end_y) { // Valid gap in this column
            if (!found_any_gap_in_any_column) {
                overall_vac_min_y = col_vac_start_y;
                overall_vac_max_y = col_vac_end_y;
                found_any_gap_in_any_column = true;
            } else {
                // Intersection: take the highest min_y and lowest max_y
                overall_vac_min_y = std::max(overall_vac_min_y, col_vac_start_y);
                overall_vac_max_y = std::min(overall_vac_max_y, col_vac_end_y);
            }
        }
    }

    if (!found_any_gap_in_any_column || overall_vac_min_y >= overall_vac_max_y) {
        std::cerr << "Warning: Could not determine a consistent vacuum channel from epsilon map." << std::endl;
        double H_sim_total = y_coords_m.back() - y_coords_m.front();
        overall_vac_min_y = y_coords_m.front() + H_sim_total / 3.0;
        overall_vac_max_y = y_coords_m.back() - H_sim_total / 3.0;
        std::cout << "Using fallback vacuum channel: [" << overall_vac_min_y << ", " << overall_vac_max_y << "]" << std::endl;
    }
    
    double gap_thickness = overall_vac_max_y - overall_vac_min_y;
    if (gap_thickness > 2.0 * h_grid_m_from_coords) { 
        overall_vac_min_y += h_grid_m_from_coords * 0.1; 
        overall_vac_max_y -= h_grid_m_from_coords * 0.1;
    }
     if (overall_vac_min_y >= overall_vac_max_y) { 
        // Revert margin if it made gap invalid
        overall_vac_min_y -= h_grid_m_from_coords * 0.1;
        overall_vac_max_y += h_grid_m_from_coords * 0.1;
        // Clamp to original bounds if necessary
        if (!found_any_gap_in_any_column) { // If still no gap, use broad fallback
             double H_sim_total = y_coords_m.back() - y_coords_m.front();
             overall_vac_min_y = y_coords_m.front() + H_sim_total / 3.0;
             overall_vac_max_y = y_coords_m.back() - H_sim_total / 3.0;
        }
     }

    return {overall_vac_min_y, overall_vac_max_y};
}

std::pair<double, double> get_field_at_point(
    double px, double py, // px, py in meters
    const std::vector<double>& x_coords, const std::vector<double>& y_coords, // x_coords, y_coords in meters
    const std::vector<std::vector<double>>& Ex_field, const std::vector<std::vector<double>>& Ey_field, // Fields in V/m
    double h_grid, int Nx, int Ny) { // h_grid in meters

    if (px < x_coords.front() || px > x_coords.back() || py < y_coords.front() || py > y_coords.back()) {
        return {0.0, 0.0}; // Outside grid
    }

    int i_idx = static_cast<int>(std::floor(px / h_grid));
    int j_idx = static_cast<int>(std::floor(py / h_grid));

    // Clamp indices to be within valid range for interpolation
    i_idx = std::max(0, std::min(i_idx, Nx - 2));
    j_idx = std::max(0, std::min(j_idx, Ny - 2));

    double x0 = x_coords[i_idx];
    double y0 = y_coords[j_idx];

    double tx = (px - x0) / h_grid;
    double ty = (py - y0) / h_grid;

    tx = std::max(0.0, std::min(1.0, tx)); // Ensure weights are in [0,1]
    ty = std::max(0.0, std::min(1.0, ty));

    // Bilinear interpolation for Ex
    double Ex_interp = (1 - tx) * (1 - ty) * Ex_field[i_idx][j_idx] +
                       tx * (1 - ty) * Ex_field[i_idx + 1][j_idx] +
                       (1 - tx) * ty * Ex_field[i_idx][j_idx + 1] +
                       tx * ty * Ex_field[i_idx + 1][j_idx + 1];

    // Bilinear interpolation for Ey
    double Ey_interp = (1 - tx) * (1 - ty) * Ey_field[i_idx][j_idx] +
                       tx * (1 - ty) * Ey_field[i_idx + 1][j_idx] +
                       (1 - tx) * ty * Ey_field[i_idx][j_idx + 1] +
                       tx * ty * Ey_field[i_idx + 1][j_idx + 1]; // Ey_field is V/m

    return {Ex_interp, Ey_interp}; // Returns V/m
}

std::pair<double, double> get_acceleration(
    double px, double py, // px, py in meters
    const std::vector<double>& x_coords, const std::vector<double>& y_coords, // x_coords, y_coords in meters
    const std::vector<std::vector<double>>& Ex_field, const std::vector<std::vector<double>>& Ey_field, // Fields in V/m
    double h_grid, int Nx, int Ny) { // h_grid in meters
    
    std::pair<double, double> E_vec = get_field_at_point(px, py, x_coords, y_coords, Ex_field, Ey_field, h_grid, Nx, Ny); // E_vec is V/m
    return {K_ACCEL * E_vec.first, K_ACCEL * E_vec.second}; // (ax, ay) in m/s^2
}

bool is_in_material_or_out_of_bounds(double px, double py, // px, py in meters
                                      const std::vector<double>& x_coords_m, // in meters
                                      const std::vector<double>& y_coords_m, // in meters
                                      const std::vector<std::vector<double>>& eps_r_map,
                                      double h_grid_m, // grid spacing in meters
                                      double L_total_sim_m, double H_total_sim_m) {
    if (px < 0.0 || px >= L_total_sim_m || py < 0.0 || py >= H_total_sim_m) {
        return true; // Out of simulation box
    }

    int i_idx = static_cast<int>(std::floor(px / h_grid_m));
    int j_idx = static_cast<int>(std::floor(py / h_grid_m));

    // Ensure indices are within the bounds of the eps_r_map
    // The map dimensions are Nx and Ny from the loaded coordinates
    int Nx_map = eps_r_map.size();
    if (Nx_map == 0) return true; // Should not happen if loaded correctly
    int Ny_map = eps_r_map[0].size();
    if (Ny_map == 0) return true; // Should not happen

    i_idx = std::max(0, std::min(i_idx, Nx_map - 1));
    j_idx = std::max(0, std::min(j_idx, Ny_map - 1));

    double permittivity_at_point = eps_r_map[i_idx][j_idx];

    // Define a threshold to distinguish silicon from vacuum
    const double EPS_SI_SIM = 11.7; // Typical relative permittivity of Silicon
    const double EPS_VAC_SIM = 1.0;  // Relative permittivity of Vacuum
    const double material_threshold = (EPS_SI_SIM + EPS_VAC_SIM) / 2.0;

    if (permittivity_at_point >= material_threshold) {
        return true; // Proton is in a material region
    }

    return false; // Proton is in vacuum or free space within bounds
}


int main() {
    std::cout << std::fixed << std::setprecision(6);

    const std::string input_base_folder = "geometria_Denti_sfasati_profondi_10kV"; // Input folder
    // const std::string output_traj_folder = input_base_folder + "/proton_trajectories"; // Removed
    // create_directory_if_not_exists(output_traj_folder); // Removed

    const std::string all_trajectories_filename = input_base_folder + "/all_proton_trajectories.csv";
    std::ofstream all_trajectories_file_stream;
    all_trajectories_file_stream.open(all_trajectories_filename, std::ios::out);
    if (!all_trajectories_file_stream.is_open()) {
        std::cerr << "Error: Could not open the consolidated trajectory file: " << all_trajectories_filename << ". Exiting." << std::endl;
        return 1;
    }
    all_trajectories_file_stream << std::scientific << std::setprecision(8);
    all_trajectories_file_stream << "proton_id,time_s,x_m,y_m,vx_m_per_s,vy_m_per_s\n";


    GeometryParameters geom;
    if (!load_geometry_params(input_base_folder + "/geometry_params.csv", geom)) {
        std::cerr << "Warning: Problem loading from geometry_params.csv. Proceeding with defaults/derivations." << std::endl;
        // geom.h, geom.x_fs, geom.x_sl will be 0.0 if file not found or keys missing.
    }

    std::vector<double> x_coords, y_coords;
    if (!load_1d_csv(input_base_folder + "/x_coordinates.csv", x_coords, true /*convert_to_meters*/) ||
        !load_1d_csv(input_base_folder + "/y_coordinates.csv", y_coords, true /*convert_to_meters*/)) {
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


    double initial_x_position_m=5e-6;
    


    std::vector<std::vector<double>> Ex_field, Ey_field; 
    std::vector<std::vector<double>> eps_r_map_data;     

    if (!load_field_csv(input_base_folder + "/electric_field_x.csv", Ex_field, Nx, Ny) ||
        !load_field_csv(input_base_folder + "/electric_field_y.csv", Ey_field, Nx, Ny) ||
        !load_permittivity_map(input_base_folder + "/permittivity.csv", eps_r_map_data, Nx, Ny)) {
        std::cerr << "Failed to load electric field or permittivity data. Exiting." << std::endl;
        return 1;
    }
    std::cout << "Data loaded successfully. Nx=" << Nx << ", Ny=" << Ny << std::endl;
    
    // Determine the x-range to search for the vacuum channel based on geometry parameters
    // This range should ideally cover the main structural part to find the "inner gap".
    double x_search_start_for_gap = geom.x_fs;
    double x_search_end_for_gap = geom.x_fs + geom.x_sl;

    if (geom.x_fs <= 1e-9 || geom.x_sl <= 1e-9 || x_search_start_for_gap < 0 || x_search_end_for_gap <= x_search_start_for_gap || x_search_end_for_gap > L_total_sim ) { 
        std::cout << "Warning: x_fs or x_sl from CSV is zero, not loaded, or defines an invalid/out-of-bounds range." << std::endl;
        std::cout << "Analyzing central 80% of x-domain for vacuum gap instead." << std::endl;
        x_search_start_for_gap = L_total_sim * 0.1;
        x_search_end_for_gap = L_total_sim * 0.9;
        // Ensure this fallback range is valid
        if (x_search_start_for_gap >= x_search_end_for_gap) {
            x_search_start_for_gap = 0.0;
            x_search_end_for_gap = L_total_sim;
        }
    }
    
    std::cout << "Searching for vacuum channel for proton initialization by analyzing x-range: [" 
              << x_search_start_for_gap * 1e6 << ", " << x_search_end_for_gap * 1e6 << "] um." << std::endl;
    
    double h_coords_for_vacuum_find = (h_from_coords_y > 1e-9) ? h_from_coords_y : h_from_coords_x;
    if (h_coords_for_vacuum_find <= 1e-9) h_coords_for_vacuum_find = h_for_simulation; // Last resort

    std::pair<double, double> vacuum_channel = find_vacuum_channel_from_map(
        eps_r_map_data, x_coords, y_coords, h_coords_for_vacuum_find,
        x_search_start_for_gap, x_search_end_for_gap
    );
    double vacuum_gap_start_y = vacuum_channel.first;
    double vacuum_gap_end_y = vacuum_channel.second;

    if (vacuum_gap_start_y >= vacuum_gap_end_y - (h_for_simulation * 0.5) ) { // Ensure gap is at least half a cell thick
        std::cerr << "Critical Error: Determined vacuum gap is too small or invalid [" << vacuum_gap_start_y << ", " << vacuum_gap_end_y << "]. Exiting." << std::endl;
        return 1;
    }
    std::cout << "Proton initial y-distribution range (m) from epsilon_map: [" << vacuum_gap_start_y << ", " << vacuum_gap_end_y << "]" << std::endl;
    
    std::vector<Proton> protons(NUM_PROTONS);
    std::mt19937 rng(std::random_device{}()); // Corrected initialization
    std::uniform_real_distribution<double> dist_y(20e-6,30e-6 ); 


    for (int i = 0; i < NUM_PROTONS; ++i) {
        protons[i].id = i; // Assign ID
        protons[i].x = initial_x_position_m; 
        protons[i].y = dist_y(rng); // y in meters
        protons[i].vx = 0.0; // Initial velocity in m/s
        protons[i].vy = 0.0; // Initial velocity in m/s
        protons[i].active = true;

        // Write initial state to the single file
        all_trajectories_file_stream << protons[i].id << "," << 0.0 << "," << protons[i].x << "," << protons[i].y << "," << protons[i].vx << "," << protons[i].vy << "\n";
    }

    double current_time = 0.0;
    int num_steps = static_cast<int>(TOTAL_SIM_TIME_S / TIME_STEP_S);
    int output_every_n_steps = static_cast<int>(std::max(1.0, OUTPUT_TIME_INTERVAL_S / TIME_STEP_S));
    int active_protons_count = NUM_PROTONS;
    int protons_reached_end_successfully = 0; // Counter for successful protons

    std::cout << "Starting simulation for " << NUM_PROTONS << " protons..." << std::endl;
    std::cout << "Total steps: " << num_steps << ", Outputting every " << output_every_n_steps << " steps." << std::endl;

    for (int step = 0; step < num_steps; ++step) {
        if (active_protons_count == 0) {
            std::cout << "All protons inactive. Stopping simulation early at step " << step << "." << std::endl;
            break;
        }

        bool needs_output_this_step = ((step + 1) % output_every_n_steps == 0);

        for (int i = 0; i < NUM_PROTONS; ++i) {
            if (!protons[i].active) continue;

            Proton& p = protons[i];
            // std::string traj_filename_for_proton = output_traj_folder + "/proton_" + std::to_string(i) + "_trajectory.csv"; // Removed
            // std::ofstream traj_file_stream_loop; // Removed

            // RK4 step
            std::pair<double, double> a1, a2, a3, a4;
            double k1x, k1y, k1vx, k1vy;
            double k2x, k2y, k2vx, k2vy;
            double k3x, k3y, k3vx, k3vy;
            double k4x, k4y, k4vx, k4vy;

            // k1
            a1 = get_acceleration(p.x, p.y, x_coords, y_coords, Ex_field, Ey_field, h_for_simulation, Nx, Ny);
            k1vx = TIME_STEP_S * a1.first;
            k1vy = TIME_STEP_S * a1.second;
            k1x = TIME_STEP_S * p.vx;
            k1y = TIME_STEP_S * p.vy;

            // k2
            a2 = get_acceleration(p.x + k1x/2.0, p.y + k1y/2.0, x_coords, y_coords, Ex_field, Ey_field, h_for_simulation, Nx, Ny);
            k2vx = TIME_STEP_S * a2.first;
            k2vy = TIME_STEP_S * a2.second;
            k2x = TIME_STEP_S * (p.vx + k1vx/2.0);
            k2y = TIME_STEP_S * (p.vy + k1vy/2.0);

            // k3
            a3 = get_acceleration(p.x + k2x/2.0, p.y + k2y/2.0, x_coords, y_coords, Ex_field, Ey_field, h_for_simulation, Nx, Ny);
            k3vx = TIME_STEP_S * a3.first;
            k3vy = TIME_STEP_S * a3.second;
            k3x = TIME_STEP_S * (p.vx + k2vx/2.0);
            k3y = TIME_STEP_S * (p.vy + k2vy/2.0);

            // k4
            a4 = get_acceleration(p.x + k3x, p.y + k3y, x_coords, y_coords, Ex_field, Ey_field, h_for_simulation, Nx, Ny);
            k4vx = TIME_STEP_S * a4.first;
            k4vy = TIME_STEP_S * a4.second;
            k4x = TIME_STEP_S * (p.vx + k3vx);
            k4y = TIME_STEP_S * (p.vy + k3vy);

            p.x += (k1x + 2.0*k2x + 2.0*k3x + k4x) / 6.0;
            p.y += (k1y + 2.0*k2y + 2.0*k3y + k4y) / 6.0;
            p.vx += (k1vx + 2.0*k2vx + 2.0*k3vx + k4vx) / 6.0;
            p.vy += (k1vy + 2.0*k2vy + 2.0*k3vy + k4vy) / 6.0;

            // bool became_inactive_this_substep = false; // Not strictly needed anymore for file closing logic

            // Check if proton reached the end successfully
            if (p.x >= L_total_sim) { // L_total_sim is in meters
                protons_reached_end_successfully++;
                p.active = false;
                active_protons_count--;
                // became_inactive_this_substep = true; // Not strictly needed
                // Write final position to the single file
                all_trajectories_file_stream << p.id << "," << (current_time + TIME_STEP_S) << "," << p.x << "," << p.y << "," << p.vx << "," << p.vy << "\n";
                continue; 
            }

            // Check if proton hit material or went out of other bounds (all checks in meters)
            if (is_in_material_or_out_of_bounds(p.x, p.y, x_coords, y_coords, eps_r_map_data, h_for_simulation, L_total_sim, H_total_sim)) {
                p.active = false;
                active_protons_count--;
                // became_inactive_this_substep = true; // Not strictly needed
                 // Write final position before deactivating due to collision/OOB to the single file
                all_trajectories_file_stream << p.id << "," << (current_time + TIME_STEP_S) << "," << p.x << "," << p.y << "," << p.vx << "," << p.vy << "\n";
            }

            // Regular output if still active and it's an output step
            if (p.active && needs_output_this_step) {
                all_trajectories_file_stream << p.id << "," << (current_time + TIME_STEP_S) << "," << p.x << "," << p.y << "," << p.vx << "," << p.vy << "\n";
            }
        }

        current_time += TIME_STEP_S;
        
        // Progress reporting remains the same
        if ((step + 1) % std::max(1, num_steps / 100) == 0 || step == num_steps -1 ) { 
             double progress = static_cast<double>(step + 1) / num_steps * 100.0;
             std::cout << "\rSimulation progress: " << std::fixed << std::setprecision(2) << progress << "% (" << active_protons_count << " active protons)" << std::flush;
        }
    }
    std::cout << "\nSimulation finished." << std::endl;

    all_trajectories_file_stream.close(); // Close the single trajectory file

    double success_percentage = 0.0;
    if (NUM_PROTONS > 0) {
        success_percentage = (static_cast<double>(protons_reached_end_successfully) / NUM_PROTONS) * 100.0;
    }
    std::cout << protons_reached_end_successfully << " out of " << NUM_PROTONS << " protons reached the end successfully." << std::endl;
    std::cout << "Success percentage: " << std::fixed << std::setprecision(2) << success_percentage << "%" << std::endl;

    std::cout << "All proton trajectories saved to '" << all_trajectories_filename << "' (data in SI units)." << std::endl;
    return 0;
}

