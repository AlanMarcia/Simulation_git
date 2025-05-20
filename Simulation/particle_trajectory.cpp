#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <random>
#include <iomanip>
#include <map>
#include <algorithm> // For std::upper_bound, std::stod

// Physical Constants
const double ELEMENTARY_CHARGE = 1.602176634e-19; // C
const double PROTON_MASS_KG = 1.67262192369e-27; // kg
const double ALPHA_PARTICLE_MASS_KG = 6.6446573357e-27; // kg (approx 4 * proton mass, more precisely)
const double M_PER_UM = 1.0e-6;
const double UM_PER_M = 1.0e6;
const double V_PER_M_FROM_V_PER_UM = 1.0e6;

enum class ParticleType {
    PROTON,
    ALPHA
};

struct Vector2D {
    double x, y;

    Vector2D(double _x = 0.0, double _y = 0.0) : x(_x), y(_y) {}

    double magnitude() const {
        return std::sqrt(x * x + y * y);
    }
};

struct Particle {
    Vector2D pos_m; // Position in meters
    Vector2D vel_m_s; // Velocity in m/s
    Vector2D acc_m_s2; // Acceleration in m/s^2
    double mass_kg;
    double charge_C;
    double time_s;

    // For storing trajectory history
    std::vector<double> history_time_s;
    std::vector<Vector2D> history_pos_um;
    std::vector<Vector2D> history_vel_m_s;
    std::vector<Vector2D> history_acc_m_s2;
    std::vector<double> history_energy_eV;

    Particle(ParticleType type) : time_s(0.0) {
        if (type == ParticleType::PROTON) {
            mass_kg = PROTON_MASS_KG;
            charge_C = ELEMENTARY_CHARGE;
        } else { // ALPHA
            mass_kg = ALPHA_PARTICLE_MASS_KG;
            charge_C = 2 * ELEMENTARY_CHARGE;
        }
    }

    void record_state() {
        history_time_s.push_back(time_s);
        history_pos_um.push_back(Vector2D(pos_m.x * UM_PER_M, pos_m.y * UM_PER_M));
        history_vel_m_s.push_back(vel_m_s);
        history_acc_m_s2.push_back(acc_m_s2);
        double kinetic_energy_J = 0.5 * mass_kg * vel_m_s.magnitude() * vel_m_s.magnitude();
        history_energy_eV.push_back(kinetic_energy_J / ELEMENTARY_CHARGE);
    }
};

// CSV Loading functions
std::vector<double> load_1d_csv_cpp(const std::string& filename) {
    std::vector<double> data;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return data;
    }
    std::string line;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string cell;
        while (std::getline(ss, cell, ',')) {
            try {
                data.push_back(std::stod(cell));
            } catch (const std::invalid_argument& ia) {
                std::cerr << "Invalid argument: " << ia.what() << " for cell " << cell << std::endl;
            } catch (const std::out_of_range& oor) {
                std::cerr << "Out of Range error: " << oor.what() << " for cell " << cell << std::endl;
            }
        }
    }
    file.close();
    return data;
}

std::vector<std::vector<double>> load_2d_csv_cpp(const std::string& filename) {
    std::vector<std::vector<double>> data;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return data;
    }
    std::string line;
    while (std::getline(file, line)) {
        std::vector<double> row;
        std::stringstream ss(line);
        std::string cell;
        while (std::getline(ss, cell, ',')) {
             try {
                row.push_back(std::stod(cell));
            } catch (const std::invalid_argument& ia) {
                std::cerr << "Invalid argument: " << ia.what() << " for cell " << cell << std::endl;
            } catch (const std::out_of_range& oor) {
                std::cerr << "Out of Range error: " << oor.what() << " for cell " << cell << std::endl;
            }
        }
        if (!row.empty()) {
            data.push_back(row);
        }
    }
    file.close();
    return data;
}

std::map<std::string, double> load_geometry_params_cpp(const std::string& filename) {
    std::map<std::string, double> params;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open geometry params file " << filename << std::endl;
        return params;
    }
    std::string line;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string key, value_str;
        if (std::getline(ss, key, ',') && std::getline(ss, value_str, ',')) {
            try {
                params[key] = std::stod(value_str);
            } catch (const std::invalid_argument& ia) {
                std::cerr << "Invalid argument in geometry params: " << ia.what() << " for key " << key << std::endl;
            } catch (const std::out_of_range& oor) {
                std::cerr << "Out of Range error in geometry params: " << oor.what() << " for key " << key << std::endl;
            }
        }
    }
    file.close();
    // Check for required keys
    std::vector<std::string> required_keys = {"y_si_layer_thick", "y_vacuum_gap_thick", "H_total", "x_free_space", "x_structure_len"};
    for(const auto& key : required_keys) {
        if(params.find(key) == params.end()) {
            std::cerr << "Error: Missing required key in geometry_params.csv: " << key << std::endl;
            params.clear(); // Indicate error by returning empty map
            return params;
        }
    }
    return params;
}


// Bilinear Interpolation for Electric Field
Vector2D interpolate_field(double x_m, double y_m,
                           const std::vector<double>& x_coords_m, const std::vector<double>& y_coords_m,
                           const std::vector<std::vector<double>>& Ex_V_m, // Ex_V_m[ix][iy]
                           const std::vector<std::vector<double>>& Ey_V_m) { // Ey_V_m[ix][iy]
    if (x_coords_m.empty() || y_coords_m.empty() || Ex_V_m.empty() || Ey_V_m.empty()) {
        return Vector2D(0,0);
    }
    // Find indices for x
    auto it_x = std::upper_bound(x_coords_m.begin(), x_coords_m.end(), x_m);
    if (it_x == x_coords_m.begin() || it_x == x_coords_m.end()) { // x_m is outside grid
        return Vector2D(0, 0); // Or handle as error/boundary condition
    }
    size_t ix2 = std::distance(x_coords_m.begin(), it_x);
    size_t ix1 = ix2 - 1;

    // Find indices for y
    auto it_y = std::upper_bound(y_coords_m.begin(), y_coords_m.end(), y_m);
    if (it_y == y_coords_m.begin() || it_y == y_coords_m.end()) { // y_m is outside grid
        return Vector2D(0, 0);
    }
    size_t iy2 = std::distance(y_coords_m.begin(), it_y);
    size_t iy1 = iy2 - 1;

    if (ix1 >= Ex_V_m.size() || ix2 >= Ex_V_m.size() || 
        iy1 >= Ex_V_m[0].size() || iy2 >= Ex_V_m[0].size() ||
        ix1 >= Ey_V_m.size() || ix2 >= Ey_V_m.size() ||
        iy1 >= Ey_V_m[0].size() || iy2 >= Ey_V_m[0].size()) {
        std::cerr << "Interpolation index out of bounds." << std::endl;
        return Vector2D(0,0);
    }


    double x1 = x_coords_m[ix1];
    double x2 = x_coords_m[ix2];
    double y1 = y_coords_m[iy1];
    double y2 = y_coords_m[iy2];

    double Ex_q11 = Ex_V_m[ix1][iy1];
    double Ex_q12 = Ex_V_m[ix1][iy2];
    double Ex_q21 = Ex_V_m[ix2][iy1];
    double Ex_q22 = Ex_V_m[ix2][iy2];

    double Ey_q11 = Ey_V_m[ix1][iy1];
    double Ey_q12 = Ey_V_m[ix1][iy2];
    double Ey_q21 = Ey_V_m[ix2][iy1];
    double Ey_q22 = Ey_V_m[ix2][iy2];

    double interpolated_Ex = 0.0, interpolated_Ey = 0.0;

    if ((x2 - x1) == 0 || (y2 - y1) == 0) { // Avoid division by zero if grid points are identical
        if ((x2-x1)==0 && (y2-y1)!=0) { // Vertical line
             interpolated_Ex = Ex_q11 + (Ex_q12 - Ex_q11) * (y_m - y1) / (y2 - y1);
             interpolated_Ey = Ey_q11 + (Ey_q12 - Ey_q11) * (y_m - y1) / (y2 - y1);
        } else if ((y2-y1)==0 && (x2-x1)!=0) { // Horizontal line
             interpolated_Ex = Ex_q11 + (Ex_q21 - Ex_q11) * (x_m - x1) / (x2 - x1);
             interpolated_Ey = Ey_q11 + (Ey_q21 - Ey_q11) * (x_m - x1) / (x2 - x1);
        } else { // Single point or invalid
            interpolated_Ex = Ex_q11;
            interpolated_Ey = Ey_q11;
        }
    } else {
        // Interpolate in x direction first
        double Ex_r1 = ((x2 - x_m) / (x2 - x1)) * Ex_q11 + ((x_m - x1) / (x2 - x1)) * Ex_q21;
        double Ex_r2 = ((x2 - x_m) / (x2 - x1)) * Ex_q12 + ((x_m - x1) / (x2 - x1)) * Ex_q22;
        // Interpolate in y direction
        interpolated_Ex = ((y2 - y_m) / (y2 - y1)) * Ex_r1 + ((y_m - y1) / (y2 - y1)) * Ex_r2;

        // Same for Ey
        double Ey_r1 = ((x2 - x_m) / (x2 - x1)) * Ey_q11 + ((x_m - x1) / (x2 - x1)) * Ey_q21;
        double Ey_r2 = ((x2 - x_m) / (x2 - x1)) * Ey_q12 + ((x_m - x1) / (x2 - x1)) * Ey_q22;
        interpolated_Ey = ((y2 - y_m) / (y2 - y1)) * Ey_r1 + ((y_m - y1) / (y2 - y1)) * Ey_r2;
    }
    
    return Vector2D(interpolated_Ex, interpolated_Ey);
}

Vector2D get_acceleration_m_s2(const Vector2D& pos_m, double charge_C, double mass_kg,
                               const std::vector<double>& x_coords_m, const std::vector<double>& y_coords_m,
                               const std::vector<std::vector<double>>& Ex_V_m,
                               const std::vector<std::vector<double>>& Ey_V_m) {
    Vector2D E_field_V_m = interpolate_field(pos_m.x, pos_m.y, x_coords_m, y_coords_m, Ex_V_m, Ey_V_m);
    return Vector2D((charge_C / mass_kg) * E_field_V_m.x, (charge_C / mass_kg) * E_field_V_m.y);
}

// RK4 Step
void rk4_step(Particle& p, double dt_s,
              const std::vector<double>& x_coords_m, const std::vector<double>& y_coords_m,
              const std::vector<std::vector<double>>& Ex_V_m,
              const std::vector<std::vector<double>>& Ey_V_m) {

    Vector2D k1_pos, k1_vel, k2_pos, k2_vel, k3_pos, k3_vel, k4_pos, k4_vel;
    Vector2D current_acc;

    // k1
    current_acc = get_acceleration_m_s2(p.pos_m, p.charge_C, p.mass_kg, x_coords_m, y_coords_m, Ex_V_m, Ey_V_m);
    k1_pos = Vector2D(p.vel_m_s.x * dt_s, p.vel_m_s.y * dt_s);
    k1_vel = Vector2D(current_acc.x * dt_s, current_acc.y * dt_s);

    // k2
    Vector2D pos_k2 = Vector2D(p.pos_m.x + k1_pos.x / 2.0, p.pos_m.y + k1_pos.y / 2.0);
    Vector2D vel_k2 = Vector2D(p.vel_m_s.x + k1_vel.x / 2.0, p.vel_m_s.y + k1_vel.y / 2.0);
    current_acc = get_acceleration_m_s2(pos_k2, p.charge_C, p.mass_kg, x_coords_m, y_coords_m, Ex_V_m, Ey_V_m);
    k2_pos = Vector2D(vel_k2.x * dt_s, vel_k2.y * dt_s);
    k2_vel = Vector2D(current_acc.x * dt_s, current_acc.y * dt_s);

    // k3
    Vector2D pos_k3 = Vector2D(p.pos_m.x + k2_pos.x / 2.0, p.pos_m.y + k2_pos.y / 2.0);
    Vector2D vel_k3 = Vector2D(p.vel_m_s.x + k2_vel.x / 2.0, p.vel_m_s.y + k2_vel.y / 2.0);
    current_acc = get_acceleration_m_s2(pos_k3, p.charge_C, p.mass_kg, x_coords_m, y_coords_m, Ex_V_m, Ey_V_m);
    k3_pos = Vector2D(vel_k3.x * dt_s, vel_k3.y * dt_s);
    k3_vel = Vector2D(current_acc.x * dt_s, current_acc.y * dt_s);

    // k4
    Vector2D pos_k4 = Vector2D(p.pos_m.x + k3_pos.x, p.pos_m.y + k3_pos.y);
    Vector2D vel_k4 = Vector2D(p.vel_m_s.x + k3_vel.x, p.vel_m_s.y + k3_vel.y);
    current_acc = get_acceleration_m_s2(pos_k4, p.charge_C, p.mass_kg, x_coords_m, y_coords_m, Ex_V_m, Ey_V_m);
    k4_pos = Vector2D(vel_k4.x * dt_s, vel_k4.y * dt_s);
    k4_vel = Vector2D(current_acc.x * dt_s, current_acc.y * dt_s);

    // Update position and velocity
    p.pos_m.x += (k1_pos.x + 2.0 * k2_pos.x + 2.0 * k3_pos.x + k4_pos.x) / 6.0;
    p.pos_m.y += (k1_pos.y + 2.0 * k2_pos.y + 2.0 * k3_pos.y + k4_pos.y) / 6.0;
    p.vel_m_s.x += (k1_vel.x + 2.0 * k2_vel.x + 2.0 * k3_vel.x + k4_vel.x) / 6.0;
    p.vel_m_s.y += (k1_vel.y + 2.0 * k2_vel.y + 2.0 * k3_vel.y + k4_vel.y) / 6.0;
    
    p.time_s += dt_s;
    // Update acceleration based on new position for storage
    p.acc_m_s2 = get_acceleration_m_s2(p.pos_m, p.charge_C, p.mass_kg, x_coords_m, y_coords_m, Ex_V_m, Ey_V_m);
}


void save_trajectory(const Particle& p, const std::string& filename) {
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }

    outfile << std::fixed << std::setprecision(6);
    outfile << "time_s,pos_x_um,pos_y_um,vel_x_m_s,vel_y_m_s,vel_mag_m_s,acc_x_m_s2,acc_y_m_s2,acc_mag_m_s2,energy_eV\n";

    for (size_t i = 0; i < p.history_time_s.size(); ++i) {
        outfile << p.history_time_s[i] << ","
                << p.history_pos_um[i].x << ","
                << p.history_pos_um[i].y << ","
                << p.history_vel_m_s[i].x << ","
                << p.history_vel_m_s[i].y << ","
                << p.history_vel_m_s[i].magnitude() << ","
                << p.history_acc_m_s2[i].x << ","
                << p.history_acc_m_s2[i].y << ","
                << p.history_acc_m_s2[i].magnitude() << ","
                << p.history_energy_eV[i] << "\n";
    }
    outfile.close();
    std::cout << "Trajectory saved to " << filename << std::endl;
}


int main() {
    // --- Simulation Parameters ---
    ParticleType particle_choice = ParticleType::PROTON; // Or ParticleType::ALPHA
    double dt_s = 1.0e-12; // Time step in seconds (e.g., 1 picosecond)
    double max_sim_time_s = 1.0e-8; // Maximum simulation time (e.g., 10 nanoseconds)
    // --- End Simulation Parameters ---


    // Load data
    std::vector<double> x_coords_um = load_1d_csv_cpp("x_coordinates.csv");
    std::vector<double> y_coords_um = load_1d_csv_cpp("y_coordinates.csv");
    std::vector<std::vector<double>> Ex_V_um = load_2d_csv_cpp("electric_field_x.csv"); // Ex[ix][iy]
    std::vector<std::vector<double>> Ey_V_um = load_2d_csv_cpp("electric_field_y.csv"); // Ey[ix][iy]
    std::map<std::string, double> geo_params = load_geometry_params_cpp("geometry_params.csv");

    if (x_coords_um.empty() || y_coords_um.empty() || Ex_V_um.empty() || Ey_V_um.empty() || geo_params.empty()) {
        std::cerr << "Failed to load one or more data files. Aborting." << std::endl;
        return 1;
    }
    
    // Convert coordinates and E-field to SI units (meters, V/m)
    std::vector<double> x_coords_m(x_coords_um.size());
    for(size_t i=0; i<x_coords_um.size(); ++i) x_coords_m[i] = x_coords_um[i] * M_PER_UM;

    std::vector<double> y_coords_m(y_coords_um.size());
    for(size_t i=0; i<y_coords_um.size(); ++i) y_coords_m[i] = y_coords_um[i] * M_PER_UM;

    std::vector<std::vector<double>> Ex_V_m(Ex_V_um.size(), std::vector<double>(Ex_V_um[0].size()));
    std::vector<std::vector<double>> Ey_V_m(Ey_V_um.size(), std::vector<double>(Ey_V_um[0].size()));

    for(size_t i=0; i<Ex_V_um.size(); ++i) {
        for(size_t j=0; j<Ex_V_um[i].size(); ++j) {
            Ex_V_m[i][j] = Ex_V_um[i][j] * V_PER_M_FROM_V_PER_UM;
        }
    }
    for(size_t i=0; i<Ey_V_um.size(); ++i) {
        for(size_t j=0; j<Ey_V_um[i].size(); ++j) {
            Ey_V_m[i][j] = Ey_V_um[i][j] * V_PER_M_FROM_V_PER_UM;
        }
    }

    // Initialize Particle
    Particle particle(particle_choice);
    particle.pos_m.x = 0.0; // Starts at x=0 meters
    particle.vel_m_s = Vector2D(0.0, 0.0); // Starts from rest

    // Random y in vacuum gap
    double y_si_thick_um = geo_params["y_si_layer_thick"];
    double y_vac_gap_thick_um = geo_params["y_vacuum_gap_thick"];
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    // Ensure y_start is within the vacuum gap, slightly offset from boundaries
    double y_offset_um = 0.01 * y_vac_gap_thick_um; // Small offset to avoid starting exactly on boundary
    double y_start_um = y_si_thick_um + y_offset_um + dis(gen) * (y_vac_gap_thick_um - 2 * y_offset_um);
    particle.pos_m.y = y_start_um * M_PER_UM;

    // Boundary definitions in meters
    double y_vac_bottom_m = y_si_thick_um * M_PER_UM;
    double y_vac_top_m = (y_si_thick_um + y_vac_gap_thick_um) * M_PER_UM;
    double sim_domain_x_min_m = x_coords_m.front();
    double sim_domain_x_max_m = x_coords_m.back();
    double sim_domain_y_min_m = y_coords_m.front();
    double sim_domain_y_max_m = y_coords_m.back();

    std::cout << "Starting simulation for " << (particle_choice == ParticleType::PROTON ? "Proton" : "Alpha Particle") << std::endl;
    std::cout << "Initial position: x = " << particle.pos_m.x * UM_PER_M << " um, y = " << particle.pos_m.y * UM_PER_M << " um" << std::endl;
    
    // Initial acceleration and state recording
    particle.acc_m_s2 = get_acceleration_m_s2(particle.pos_m, particle.charge_C, particle.mass_kg, x_coords_m, y_coords_m, Ex_V_m, Ey_V_m);
    particle.record_state();

    // Simulation Loop
    int step_count = 0;
    while (particle.time_s < max_sim_time_s) {
        rk4_step(particle, dt_s, x_coords_m, y_coords_m, Ex_V_m, Ey_V_m);
        particle.record_state();
        step_count++;

        // Boundary checks
        if (particle.pos_m.y <= y_vac_bottom_m || particle.pos_m.y >= y_vac_top_m) {
            std::cout << "Particle hit Si layer at t = " << particle.time_s << " s. Position: (" 
                      << particle.pos_m.x*UM_PER_M << ", " << particle.pos_m.y*UM_PER_M << ") um." << std::endl;
            break;
        }
        if (particle.pos_m.x < sim_domain_x_min_m || particle.pos_m.x > sim_domain_x_max_m ||
            particle.pos_m.y < sim_domain_y_min_m || particle.pos_m.y > sim_domain_y_max_m) {
            std::cout << "Particle left simulation domain at t = " << particle.time_s << " s. Position: (" 
                      << particle.pos_m.x*UM_PER_M << ", " << particle.pos_m.y*UM_PER_M << ") um." << std::endl;
            break;
        }
        if (step_count % 1000 == 0) { // Optional: print progress
             std::cout << "Time: " << particle.time_s << " s, Pos (um): (" << particle.pos_m.x*UM_PER_M << ", " << particle.pos_m.y*UM_PER_M 
                       << "), Vel (m/s): (" << particle.vel_m_s.x << ", " << particle.vel_m_s.y << ")" << std::endl;
        }
    }
    
    if (particle.time_s >= max_sim_time_s) {
        std::cout << "Maximum simulation time reached." << std::endl;
    }
    std::cout << "Simulation finished after " << step_count << " steps." << std::endl;

    // Save trajectory
    std::string output_filename = (particle_choice == ParticleType::PROTON) ? "trajectory_proton.csv" : "trajectory_alpha.csv";
    save_trajectory(particle, output_filename);

    return 0;
}

