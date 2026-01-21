#include<iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <random>
#include <algorithm> // For std::min, std::max
#include <sstream>   // For std::stringstream
#include <omp.h>     // For OpenMP

int main(){

    //apro file di campo
    std::ifstream file("field.fld");
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file field.fld" << std::endl;
        return 1;
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
            // Converti le coordinate da micrometri a metri per coerenza con il resto della simulazione
            X_values.push_back(row_values[0] * 1e-6);
            Y_values.push_back(row_values[1] * 1e-6);
            Ex_values.push_back(row_values[3]); // V/m
            Ey_values.push_back(row_values[4]); // V/m
        }
    }


    file.close();

    // Output the number of values read for verification
    std::cout << "Read " << Ex_values.size() << " field values." << std::endl;

    if (Ex_values.empty() || Ey_values.empty()) {
        std::cerr << "No field data loaded. Check field.fld format." << std::endl;
        return 1;
    }

    // Debug: min/max grezzi
    auto minmax_x_raw = std::minmax_element(X_values.begin(), X_values.end());
    auto minmax_y_raw = std::minmax_element(Y_values.begin(), Y_values.end());
    std::cout << "Raw x range: [" << *minmax_x_raw.first << ", " << *minmax_x_raw.second << "] m\n";
    std::cout << "Raw y range: [" << *minmax_y_raw.first << ", " << *minmax_y_raw.second << "] m\n";

    // Ricostruisci la griglia 2D usando i valori unici di x e y (ordinati).
    std::vector<double> grid_x = X_values;
    std::vector<double> grid_y = Y_values;
    std::sort(grid_x.begin(), grid_x.end());
    grid_x.erase(std::unique(grid_x.begin(), grid_x.end()), grid_x.end());
    std::sort(grid_y.begin(), grid_y.end());
    grid_y.erase(std::unique(grid_y.begin(), grid_y.end()), grid_y.end());

    const size_t Nx = grid_x.size();
    const size_t Ny = grid_y.size();
    if (Nx == 0 || Ny == 0 || Ex_values.size() != Nx * Ny) {
        std::cerr << "Unable to deduce grid dimensions: Nx=" << Nx << " Ny=" << Ny
                  << " entries=" << Ex_values.size() << std::endl;
        return 1;
    }

    // Costruisci griglie Ex/Ey indicizzate (iy * Nx + ix)
    std::vector<double> Ex_grid(Nx * Ny, 0.0);
    std::vector<double> Ey_grid(Nx * Ny, 0.0);
    for (size_t k = 0; k < Ex_values.size(); ++k) {
        double x = X_values[k];
        double y = Y_values[k];
        auto itx = std::lower_bound(grid_x.begin(), grid_x.end(), x);
        auto ity = std::lower_bound(grid_y.begin(), grid_y.end(), y);
        if (itx == grid_x.end() || ity == grid_y.end()) {
            continue; // dovrebbe essere raro
        }
        size_t ix = static_cast<size_t>(std::distance(grid_x.begin(), itx));
        size_t iy = static_cast<size_t>(std::distance(grid_y.begin(), ity));
        Ex_grid[iy * Nx + ix] = Ex_values[k];
        Ey_grid[iy * Nx + ix] = Ey_values[k];
    }

    // Debug: stampa limiti della griglia e parametri iniziali
    std::cout << "Grid Nx=" << Nx << " Ny=" << Ny << "\n";
    std::cout << "Unique x count: " << grid_x.size() << " Unique y count: " << grid_y.size() << "\n";
    std::cout << "Grid x range: [" << grid_x.front() << ", " << grid_x.back() << "] m\n";
    std::cout << "Grid y range: [" << grid_y.front() << ", " << grid_y.back() << "] m\n";
    std::cout << "Start x = -5e-6 m, start y in [-20e-6, 20e-6] m" << std::endl;

    // Runge-Kutta 4 per protoni in campo elettrico
    double dt = 1e-14; // Time step in seconds
    double max_time = 1e-7; // Maximum simulation time in seconds
    double charge = 1.602e-19; // Charge of proton in Coulombs
    double initial_energy = 1000; // Initial energy in eV
    double mass = 1.672e-27; // Mass of proton in kg
    int n_protons = 50; // Numero di protoni da simulare

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-10, 10); // Angolo uniforme tra -10 e 10 gradi
    std::uniform_real_distribution<> dis_y(-40, 40.0); // y iniziale in micron
    std::ofstream traj_file("proton_trajectories.csv");
    traj_file << "proton_id,time_s,pos_x_m,pos_y_m,vel_x_m_s,vel_y_m_s\n";

    for (int p = 0; p < n_protons; ++p) {
        double pos_x = -9*1e-6; // Initial x position in meters
        double pos_y = dis_y(gen) * 1e-6;   // Initial y position in meters
        double vel_magnitude = std::sqrt((2 * initial_energy * 1.602e-19) / mass); // Initial velocity magnitude in m/s
        double angle_deg = dis(gen);
        double angle_rad = angle_deg * (M_PI / 180.0);
        double vel_x = vel_magnitude * std::cos(angle_rad);
        double vel_y = vel_magnitude * std::sin(angle_rad);
        double time = 0.0;

        while (time < max_time) {
            // Se esce dal dominio della griglia, interrompi (evita campi extrapolati non affidabili)
            if (pos_x < grid_x.front() || pos_x > grid_x.back() || pos_y < grid_y.front() || pos_y > grid_y.back()) {
                std::cout << "Proton " << p << " left the simulation domain at time " << time << " s." << std::endl;
                break;
            }

            // Trova il punto di griglia più vicino in x e y
            auto it_x = std::lower_bound(grid_x.begin(), grid_x.end(), pos_x);
            auto it_y = std::lower_bound(grid_y.begin(), grid_y.end(), pos_y);

            if (it_x == grid_x.end() || it_y == grid_y.end()) {
                std::cout << "Proton " << p << " left the simulation domain at time " << time << " s." << std::endl;
                break;
            }

            size_t ix = static_cast<size_t>(std::distance(grid_x.begin(), it_x));
            size_t iy = static_cast<size_t>(std::distance(grid_y.begin(), it_y));
            // Clamp per sicurezza
            ix = std::min(ix, Nx - 1);
            iy = std::min(iy, Ny - 1);
            size_t index = iy * Nx + ix;

            double Ex = Ex_grid[index];
            double Ey = Ey_grid[index];

            // Calculate acceleration
            double acc_x = (charge / mass) * Ex;
            double acc_y = (charge / mass) * Ey;

            // RK4 Integration (campo considerato costante su dt)
            double k1_vx = acc_x * dt;
            double k1_vy = acc_y * dt;
            double k1_px = vel_x * dt;
            double k1_py = vel_y * dt;

            double k2_vx = acc_x * dt;
            double k2_vy = acc_y * dt;
            double k2_px = (vel_x + k1_vx / 2.0) * dt;
            double k2_py = (vel_y + k1_vy / 2.0) * dt;

            double k3_vx = acc_x * dt;
            double k3_vy = acc_y * dt;
            double k3_px = (vel_x + k2_vx / 2.0) * dt;
            double k3_py = (vel_y + k2_vy / 2.0) * dt;

            double k4_vx = acc_x * dt;
            double k4_vy = acc_y * dt;
            double k4_px = (vel_x + k3_vx) * dt;
            double k4_py = (vel_y + k3_vy) * dt;

            vel_x += (k1_vx + 2.0 * k2_vx + 2.0 * k3_vx + k4_vx) / 6.0;
            vel_y += (k1_vy + 2.0 * k2_vy + 2.0 * k3_vy + k4_vy) / 6.0;
            pos_x += (k1_px + 2.0 * k2_px + 2.0 * k3_px + k4_px) / 6.0;
            pos_y += (k1_py + 2.0 * k2_py + 2.0 * k3_py + k4_py) / 6.0;
            time += dt;

            traj_file << p << "," << time << "," << pos_x << "," << pos_y << "," << vel_x << "," << vel_y << "\n";
        }
    }

    return 0;
}