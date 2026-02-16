#include<iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <random>
#include <algorithm>
#include <sstream>
#include <omp.h>
#include <filesystem>
#include <mutex>

namespace fs = std::filesystem;

// Mutex per sincronizzazione dell'output e file operations
std::mutex cout_mutex;
std::mutex file_mutex;

// Struttura per parametri della simulazione
struct SimulationParams {
    double dt = 1e-12;
    double max_time = 3e-8;
    double charge = 1.602e-19;
    double initial_energy = 10000;
    double mass = 1.672e-27;
    int n_protons = 10;
    double start_x = -265e-6;
    double start_y_min = -400e-6;
    double start_y_max = 400e-6;
};

// Struttura per dati della griglia di campo
struct FieldGrid {
    std::vector<double> grid_x;
    std::vector<double> grid_y;
    std::vector<double> Ex_grid;
    std::vector<double> Ey_grid;
    size_t Nx;
    size_t Ny;
};

// Funzione per leggere il file di campo
FieldGrid readFieldFile(const std::string& filename) {
    FieldGrid grid;
    
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return grid;
    }

    std::string line;
    std::getline(file, line);
    std::getline(file, line);

    std::vector<double> Ex_values;
    std::vector<double> Ey_values;
    std::vector<double> X_values;
    std::vector<double> Y_values;
    
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::vector<double> row_values;
        double value;
        while (iss >> value) {
            row_values.push_back(value);
        }
        if (row_values.size() >= 6) {
            X_values.push_back(row_values[0]);  // già in metri (SI)
            Y_values.push_back(row_values[1]);  // già in metri (SI)
            Ex_values.push_back(row_values[3]);
            Ey_values.push_back(row_values[4]);
        }
    }
    file.close();

    if (Ex_values.empty() || Ey_values.empty()) {
        std::cerr << "No field data loaded from " << filename << std::endl;
        return grid;
    }

    // Ricostruisci la griglia 2D
    grid.grid_x = X_values;
    grid.grid_y = Y_values;
    std::sort(grid.grid_x.begin(), grid.grid_x.end());
    grid.grid_x.erase(std::unique(grid.grid_x.begin(), grid.grid_x.end()), grid.grid_x.end());
    std::sort(grid.grid_y.begin(), grid.grid_y.end());
    grid.grid_y.erase(std::unique(grid.grid_y.begin(), grid.grid_y.end()), grid.grid_y.end());

    grid.Nx = grid.grid_x.size();
    grid.Ny = grid.grid_y.size();

    if (grid.Nx == 0 || grid.Ny == 0 || Ex_values.size() != grid.Nx * grid.Ny) {
        std::cerr << "Grid dimension mismatch in " << filename << std::endl;
        return grid;
    }

    // Costruisci griglie Ex/Ey indicizzate
    grid.Ex_grid.assign(grid.Nx * grid.Ny, 0.0);
    grid.Ey_grid.assign(grid.Nx * grid.Ny, 0.0);
    
    for (size_t k = 0; k < Ex_values.size(); ++k) {
        double x = X_values[k];
        double y = Y_values[k];
        auto itx = std::lower_bound(grid.grid_x.begin(), grid.grid_x.end(), x);
        auto ity = std::lower_bound(grid.grid_y.begin(), grid.grid_y.end(), y);
        if (itx == grid.grid_x.end() || ity == grid.grid_y.end()) continue;
        
        size_t ix = static_cast<size_t>(std::distance(grid.grid_x.begin(), itx));
        size_t iy = static_cast<size_t>(std::distance(grid.grid_y.begin(), ity));
        grid.Ex_grid[iy * grid.Nx + ix] = Ex_values[k];
        grid.Ey_grid[iy * grid.Nx + ix] = Ey_values[k];
    }

    return grid;
}

// Funzione per simulare un singolo protone
struct ProtonTrajectory {
    int proton_id;
    std::vector<std::tuple<double, double, double, double, double, double>> data; // time, pos_x, pos_y, vel_x, vel_y, ...
};

void simulateProton(int proton_id, const FieldGrid& grid, const SimulationParams& params,
                    std::vector<ProtonTrajectory>& trajectories, unsigned int seed) {
    ProtonTrajectory traj;
    traj.proton_id = proton_id;

    std::mt19937 gen(seed + proton_id);
    std::uniform_real_distribution<> dis_y(params.start_y_min, params.start_y_max);

    double pos_x = params.start_x;
    double pos_y = dis_y(gen);
    double vel_magnitude = std::sqrt((2 * params.initial_energy * 1.602e-19) / params.mass);
    double vel_x = vel_magnitude;
    double vel_y = 0.0;
    double time = 0.0;

    while (pos_x < 5.7e-3 && time < params.max_time)  {
        double Ex = 0.0;
        double Ey = 0.0;

        if (pos_x >= grid.grid_x.front() && pos_x <= grid.grid_x.back() &&
            pos_y >= grid.grid_y.front() && pos_y <= grid.grid_y.back()) {
            
            auto it_x = std::lower_bound(grid.grid_x.begin(), grid.grid_x.end(), pos_x);
            auto it_y = std::lower_bound(grid.grid_y.begin(), grid.grid_y.end(), pos_y);

            if (it_x != grid.grid_x.end() && it_y != grid.grid_y.end()) {
                size_t ix = static_cast<size_t>(std::distance(grid.grid_x.begin(), it_x));
                size_t iy = static_cast<size_t>(std::distance(grid.grid_y.begin(), it_y));
                ix = std::min(ix, grid.Nx - 1);
                iy = std::min(iy, grid.Ny - 1);
                size_t index = iy * grid.Nx + ix;

                Ex = grid.Ex_grid[index];
                Ey = grid.Ey_grid[index];
            }
        }

        double acc_x = (params.charge / params.mass) * Ex;
        double acc_y = (params.charge / params.mass) * Ey;

        // RK4 Integration
        double k1_vx = acc_x * params.dt;
        double k1_vy = acc_y * params.dt;
        double k1_px = vel_x * params.dt;
        double k1_py = vel_y * params.dt;

        double k2_vx = acc_x * params.dt;
        double k2_vy = acc_y * params.dt;
        double k2_px = (vel_x + k1_vx / 2.0) * params.dt;
        double k2_py = (vel_y + k1_vy / 2.0) * params.dt;

        double k3_vx = acc_x * params.dt;
        double k3_vy = acc_y * params.dt;
        double k3_px = (vel_x + k2_vx / 2.0) * params.dt;
        double k3_py = (vel_y + k2_vy / 2.0) * params.dt;

        double k4_vx = acc_x * params.dt;
        double k4_vy = acc_y * params.dt;
        double k4_px = (vel_x + k3_vx) * params.dt;
        double k4_py = (vel_y + k3_vy) * params.dt;

        vel_x += (k1_vx + 2.0 * k2_vx + 2.0 * k3_vx + k4_vx) / 6.0;
        vel_y += (k1_vy + 2.0 * k2_vy + 2.0 * k3_vy + k4_vy) / 6.0;
        pos_x += (k1_px + 2.0 * k2_px + 2.0 * k3_px + k4_px) / 6.0;
        pos_y += (k1_py + 2.0 * k2_py + 2.0 * k3_py + k4_py) / 6.0;
        time += params.dt;

        traj.data.push_back({time, pos_x, pos_y, vel_x, vel_y, 0.0});
    }

    trajectories[proton_id] = traj;
}

// Funzione per salvare le traiettorie su file
void saveTrajectories(const std::string& output_path, const std::vector<ProtonTrajectory>& trajectories) {
    std::ofstream traj_file(output_path);
    traj_file << "proton_id,time_s,pos_x_m,pos_y_m,vel_x_m_s,vel_y_m_s\n";

    for (const auto& traj : trajectories) {
        for (const auto& point : traj.data) {
            traj_file << traj.proton_id << ","
                      << std::get<0>(point) << ","
                      << std::get<1>(point) << ","
                      << std::get<2>(point) << ","
                      << std::get<3>(point) << ","
                      << std::get<4>(point) << "\n";
        }
    }
    traj_file.close();
}

// Funzione principale per processare un singolo file di campo
void processFieldFile(const std::string& field_file_path, const std::string& base_output_dir,
                      const SimulationParams& params, int file_index) {
    std::string filename = fs::path(field_file_path).filename().string();
    std::string output_dir = base_output_dir + "/" + filename.substr(0, filename.find_last_of('.'));

    // Crea la cartella di output
    fs::create_directories(output_dir);

    {
        std::lock_guard<std::mutex> lock(cout_mutex);
        std::cout << "[File " << file_index << "] Processing: " << filename << std::endl;
    }

    // Leggi il file di campo
    FieldGrid grid = readFieldFile(field_file_path);
    
    if (grid.Nx == 0 || grid.Ny == 0) {
        std::lock_guard<std::mutex> lock(cout_mutex);
        std::cerr << "[File " << file_index << "] ERROR: Could not load grid from " << filename << std::endl;
        return;
    }

    {
        std::lock_guard<std::mutex> lock(cout_mutex);
        std::cout << "[File " << file_index << "] Grid loaded: Nx=" << grid.Nx << " Ny=" << grid.Ny << std::endl;
    }

    // Alloca spazio per le traiettorie
    std::vector<ProtonTrajectory> trajectories(params.n_protons);

    // Simula i protoni in parallelo
    unsigned int seed = std::random_device{}();
    
    #pragma omp parallel for schedule(dynamic)
    for (int p = 0; p < params.n_protons; ++p) {
        simulateProton(p, grid, params, trajectories, seed);
    }

    // Salva le traiettorie
    std::string output_file = output_dir + "/proton_trajectories.csv";
    saveTrajectories(output_file, trajectories);

    {
        std::lock_guard<std::mutex> lock(cout_mutex);
        std::cout << "[File " << file_index << "] Completed: " << params.n_protons << " protons simulated." << std::endl;
        std::cout << "[File " << file_index << "] Output saved to: " << output_file << std::endl;
    }
}

int main() {
    std::string base_field_dir = "export_field";
    std::string base_output_dir = "proton_trajectories_output";

    // Crea la cartella di output principale
    fs::create_directories(base_output_dir);

    // Parametri della simulazione
    SimulationParams params;

    // Raccogli tutti i file .fld
    std::vector<std::string> field_files;
    if (!fs::exists(base_field_dir)) {
        std::cerr << "Error: Directory " << base_field_dir << " does not exist." << std::endl;
        return 1;
    }

    for (const auto& entry : fs::directory_iterator(base_field_dir)) {
        if (entry.is_regular_file() && entry.path().extension() == ".fld") {
            field_files.push_back(entry.path().string());
        }
    }

    if (field_files.empty()) {
        std::cerr << "No .fld files found in " << base_field_dir << std::endl;
        return 1;
    }

    std::sort(field_files.begin(), field_files.end());

    std::cout << "Found " << field_files.size() << " field files to process." << std::endl;
    std::cout << "Starting parallel batch processing with " << omp_get_max_threads() << " threads..." << std::endl;

    // Processa i file in parallelo
    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < field_files.size(); ++i) {
        processFieldFile(field_files[i], base_output_dir, params, static_cast<int>(i + 1));
    }

    std::cout << "\n=== Batch processing completed! ===" << std::endl;
    std::cout << "All trajectories saved in: " << base_output_dir << std::endl;

    return 0;
}
