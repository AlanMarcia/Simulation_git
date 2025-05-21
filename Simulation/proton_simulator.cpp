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
const int NUM_PROTONS = 1000;
const double TIME_STEP_S = 1e-12;       // Time step in seconds (SI)
const double TOTAL_SIM_TIME_S = 1e-8;  // Total simulation time in seconds (SI)
const double OUTPUT_TIME_INTERVAL_S = 1e-10; // Interval for writing trajectory data (SI)
// Initial X position will be set in main after loading geometry, in meters

// --- Structures ---
struct Proton {
    double x, y;     // Position in meters (SI)
    double vx, vy;   // Velocity in m/s (SI)
    bool active;
    std::ofstream trajectory_file;
};

struct GeometryParameters {
    double h;         // Grid spacing in meters (SI)
    double x_fs;      // x_free_space in meters (SI)
    double x_sl;      // x_structure_len in meters (SI)
    double y_slt;     // y_si_layer_thick in meters (SI)
    double y_vgt;     // y_vacuum_gap_thick in meters (SI)
    double H_tot_geom; // H_total (from geometry file) in meters (SI)
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

bool load_geometry_params(const std::string& filename, GeometryParameters& geom) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return false;
    }
    std::string line;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string key, value_str;
        std::getline(ss, key, ',');
        std::getline(ss, value_str, ',');
        double value_um = std::stod(value_str); // Value from CSV is in µm
        double value_m = value_um * 1.0e-6;     // Convert to meters

        if (key == "h") geom.h = value_m;
        else if (key == "x_free_space") geom.x_fs = value_m;
        else if (key == "x_structure_len") geom.x_sl = value_m;
        else if (key == "y_si_layer_thick") geom.y_slt = value_m;
        else if (key == "y_vacuum_gap_thick") geom.y_vgt = value_m;
        else if (key == "H_total") geom.H_tot_geom = value_m;
    }
    file.close();
    return true;
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
                                      const GeometryParameters& geom, // geom parameters in meters
                                      double L_total_sim, double H_total_sim) { // L_total_sim, H_total_sim in meters
    if (px < 0.0 || px > L_total_sim || py < 0.0 || py > H_total_sim) {
        return true; // Out of simulation box
    }

    double x_struct_start = geom.x_fs;
    double x_struct_end = geom.x_fs + geom.x_sl;

    if (px >= x_struct_start && px <= x_struct_end) {
        double y_si_bottom_layer_end = geom.y_slt;
        double y_top_si_layer_start = geom.y_slt + geom.y_vgt;

        if (py <= y_si_bottom_layer_end) return true; // In bottom Si
        if (py >= y_top_si_layer_start) return true;  // In top Si
    }
    return false; // In vacuum or free space outside structure
}


int main() {
    std::cout << std::fixed << std::setprecision(6);

    const std::string input_base_folder = "geometria_piana";
    const std::string output_traj_folder = input_base_folder + "/proton_trajectories";
    create_directory_if_not_exists(output_traj_folder);

    GeometryParameters geom;
    if (!load_geometry_params(input_base_folder + "/geometry_params.csv", geom)) {
        std::cerr << "Failed to load geometry parameters. Exiting." << std::endl;
        return 1;
    }

    // Set initial X position based on loaded geometry to be inside the active region
    const double initial_x_position_m = geom.x_fs + 1.0e-6; // 1 µm (converted to m) into the structure start
    std::cout << "Protons will start at x = " << initial_x_position_m * 1e6 << " um." << std::endl; // Display in um for clarity

    std::vector<double> x_coords, y_coords;
    // Load coordinates and convert them from µm (in CSV) to m for simulation
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
    double L_total_sim = x_coords.back(); // Now in meters
    double H_total_sim = y_coords.back(); // Now in meters

     // Verify h from geom matches h from coordinates if needed (all in meters now)
    if (Nx > 1 && Ny > 1 && (std::abs(geom.h - (x_coords[1] - x_coords[0])) > 1e-12 || std::abs(geom.h - (y_coords[1] - y_coords[0])) > 1e-12)) {
        std::cout << "Warning: h from geometry_params.csv (" << geom.h * 1e6 // Display in um
                  << " um) differs from h derived from coordinates (" << (x_coords[1]-x_coords[0]) * 1e6 << " um, " // Display in um
                  << (y_coords[1]-y_coords[0]) * 1e6 << " um). Using h from geometry_params.csv for interpolation." << std::endl; // Display in um
    }


    std::vector<std::vector<double>> Ex_field, Ey_field; // Will be loaded in V/m
    if (!load_field_csv(input_base_folder + "/electric_field_x.csv", Ex_field, Nx, Ny) ||
        !load_field_csv(input_base_folder + "/electric_field_y.csv", Ey_field, Nx, Ny)) {
        std::cerr << "Failed to load electric field data. Exiting." << std::endl;
        return 1;
    }
    std::cout << "Data loaded successfully. Nx=" << Nx << ", Ny=" << Ny << std::endl;

    std::vector<Proton> protons(NUM_PROTONS);
    std::mt19937 rng(std::random_device{}());
    // Ensure y-distribution is within the vacuum gap correctly defined by geometry (all in meters)
    double vacuum_gap_start_y = geom.y_slt;
    double vacuum_gap_end_y = geom.y_slt + geom.y_vgt;
    if (vacuum_gap_start_y >= vacuum_gap_end_y) {
        std::cerr << "Error: Vacuum gap has zero or negative thickness. Check y_slt and y_vgt (after conversion to meters)." << std::endl;
        std::cerr << "y_slt (m): " << geom.y_slt << ", y_vgt (m): " << geom.y_vgt << std::endl;
        // Fallback or safe default if parameters are problematic
        if (geom.H_tot_geom > 0.0) { // if H_total_geom is available and positive (in meters)
             vacuum_gap_start_y = geom.H_tot_geom / 3.0;
             vacuum_gap_end_y = geom.H_tot_geom * 2.0 / 3.0;
             std::cout << "Warning: Using fallback y-distribution for protons due to problematic gap params." << std::endl;
        } else { // Absolute fallback
            vacuum_gap_start_y = 1.0;
            vacuum_gap_end_y = 2.0;
             std::cout << "Warning: Using ABSOLUTE fallback y-distribution for protons." << std::endl;
        }
         if (vacuum_gap_start_y >= vacuum_gap_end_y) { // Final check on fallback
            std::cerr << "Critical Error: Fallback y-distribution also invalid. Exiting." << std::endl;
            return 1;
         }
    }
    std::uniform_real_distribution<double> dist_y(vacuum_gap_start_y, vacuum_gap_end_y); // Range in meters


    for (int i = 0; i < NUM_PROTONS; ++i) {
        protons[i].x = 5e-6; // Initial x position in meters (5 µm)
        //initial_x_position_m; 
        protons[i].y = dist_y(rng); // y in meters
        protons[i].vx = 0.0; // Initial velocity in m/s
        protons[i].vy = 0.0; // Initial velocity in m/s
        protons[i].active = true;

        std::string traj_filename = output_traj_folder + "/proton_" + std::to_string(i) + "_trajectory.csv";
        protons[i].trajectory_file.open(traj_filename);
        if (!protons[i].trajectory_file.is_open()) {
            std::cerr << "Error: Could not open trajectory file " << traj_filename << std::endl;
            protons[i].active = false; // Cannot save, so don't simulate
            continue;
        }
        protons[i].trajectory_file << std::scientific << std::setprecision(8); // Use scientific for SI potentially small/large numbers
        protons[i].trajectory_file << "time_s,x_m,y_m,vx_m_per_s,vy_m_per_s\n";
        protons[i].trajectory_file << 0.0 << "," << protons[i].x << "," << protons[i].y << "," << protons[i].vx << "," << protons[i].vy << "\n";
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

        for (int i = 0; i < NUM_PROTONS; ++i) {
            if (!protons[i].active) continue;

            Proton& p = protons[i];

            // RK4 step
            std::pair<double, double> a1, a2, a3, a4;
            double k1x, k1y, k1vx, k1vy;
            double k2x, k2y, k2vx, k2vy;
            double k3x, k3y, k3vx, k3vy;
            double k4x, k4y, k4vx, k4vy;

            // k1
            a1 = get_acceleration(p.x, p.y, x_coords, y_coords, Ex_field, Ey_field, geom.h, Nx, Ny);
            k1vx = TIME_STEP_S * a1.first;
            k1vy = TIME_STEP_S * a1.second;
            k1x = TIME_STEP_S * p.vx;
            k1y = TIME_STEP_S * p.vy;

            // k2
            a2 = get_acceleration(p.x + k1x/2.0, p.y + k1y/2.0, x_coords, y_coords, Ex_field, Ey_field, geom.h, Nx, Ny);
            k2vx = TIME_STEP_S * a2.first;
            k2vy = TIME_STEP_S * a2.second;
            k2x = TIME_STEP_S * (p.vx + k1vx/2.0);
            k2y = TIME_STEP_S * (p.vy + k1vy/2.0);

            // k3
            a3 = get_acceleration(p.x + k2x/2.0, p.y + k2y/2.0, x_coords, y_coords, Ex_field, Ey_field, geom.h, Nx, Ny);
            k3vx = TIME_STEP_S * a3.first;
            k3vy = TIME_STEP_S * a3.second;
            k3x = TIME_STEP_S * (p.vx + k2vx/2.0);
            k3y = TIME_STEP_S * (p.vy + k2vy/2.0);

            // k4
            a4 = get_acceleration(p.x + k3x, p.y + k3y, x_coords, y_coords, Ex_field, Ey_field, geom.h, Nx, Ny);
            k4vx = TIME_STEP_S * a4.first;
            k4vy = TIME_STEP_S * a4.second;
            k4x = TIME_STEP_S * (p.vx + k3vx);
            k4y = TIME_STEP_S * (p.vy + k3vy);

            p.x += (k1x + 2.0*k2x + 2.0*k3x + k4x) / 6.0;
            p.y += (k1y + 2.0*k2y + 2.0*k3y + k4y) / 6.0;
            p.vx += (k1vx + 2.0*k2vx + 2.0*k3vx + k4vx) / 6.0;
            p.vy += (k1vy + 2.0*k2vy + 2.0*k3vy + k4vy) / 6.0;

            // Check if proton reached the end successfully
            if (p.x >= L_total_sim) { // L_total_sim is in meters
                protons_reached_end_successfully++;
                p.active = false;
                if (p.trajectory_file.is_open()) { 
                    p.trajectory_file << (current_time + TIME_STEP_S) << "," << p.x << "," << p.y << "," << p.vx << "," << p.vy << "\n";
                    p.trajectory_file.close();
                }
                active_protons_count--;
                continue; 
            }

            // Check if proton hit material or went out of other bounds (all checks in meters)
            if (is_in_material_or_out_of_bounds(p.x, p.y, geom, L_total_sim, H_total_sim)) {
                p.active = false;
                if (p.trajectory_file.is_open()) {
                    p.trajectory_file.close();
                }
                active_protons_count--;
            }
        }

        current_time += TIME_STEP_S;

        if ((step + 1) % output_every_n_steps == 0) {
            for (int i = 0; i < NUM_PROTONS; ++i) {
                if (protons[i].active) {
                    Proton& p = protons[i];
                    // Output data in SI units (meters, m/s)
                    p.trajectory_file << current_time << "," << p.x << "," << p.y << "," << p.vx << "," << p.vy << "\n";
                }
            }
        }
        
        if ((step + 1) % std::max(1, num_steps / 100) == 0 || step == num_steps -1 ) { // Print progress roughly every 1%
             double progress = static_cast<double>(step + 1) / num_steps * 100.0;
             std::cout << "\rSimulation progress: " << std::fixed << std::setprecision(2) << progress << "% (" << active_protons_count << " active protons)" << std::flush;
        }
    }
    std::cout << "\nSimulation finished." << std::endl;

    for (int i = 0; i < NUM_PROTONS; ++i) {
        if (protons[i].trajectory_file.is_open()) {
            protons[i].trajectory_file.close();
        }
    }

    double success_percentage = 0.0;
    if (NUM_PROTONS > 0) {
        success_percentage = (static_cast<double>(protons_reached_end_successfully) / NUM_PROTONS) * 100.0;
    }
    std::cout << protons_reached_end_successfully << " out of " << NUM_PROTONS << " protons reached the end successfully." << std::endl;
    std::cout << "Success percentage: " << std::fixed << std::setprecision(2) << success_percentage << "%" << std::endl;

    std::cout << "Proton trajectories saved to '" << output_traj_folder << "' (data in SI units)." << std::endl;
    return 0;
}

