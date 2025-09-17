#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <random>
#include <string>
#include <omp.h>  // Per OpenMP

// Definisce M_PI se non è già definito
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Costanti fisiche
const double c = 2.998e8;           // Velocità della luce [m/s]
const double eps0 = 8.854e-12;     // Permittività del vuoto [F/m]
const double mu0 = 4.0 * M_PI * 1e-7; // Permeabilità del vuoto [H/m]
const double e = 1.602e-19;        // Carica elementare [C]
const double mp = 1.673e-27;       // Massa del protone [kg]
const double me = 9.109e-31;       // Massa dell'elettrone [kg]

// Parametri della simulazione
const double L = 300e-6;           // Lunghezza totale [m] (300 µm)
const double H = 50e-6;            // Altezza totale [m] (50 µm)
const int Nx = 600;                // Punti griglia in x
const int Ny = 100;                // Punti griglia in y
const double dx = L / (Nx - 1);    // Spaziatura griglia x
const double dy = H / (Ny - 1);    // Spaziatura griglia y
const double dt = 0.5 * std::min(dx, dy) / c; // Passo temporale (stabilità CFL)

// Parametri laser
const double lambda0 = 2.0e-6;     // Lunghezza d'onda laser [m]
const double omega0 = 2 * M_PI * c / lambda0; // Frequenza angolare
const double E0_max = 1e9;         // Campo elettrico massimo [V/m]
const double pulse_duration = 50e-15; // Durata impulso [s]

// Parametri materiali
const double eps_sio2 = 3.9;       // Permittività relativa SiO2 (ossido di silicio)
const double eps_vacuum = 1.0;     // Permittività relativa vuoto

// Struttura per le particelle - ottimizzata per l'accesso vettorizzato
struct Particle {
    double x, y;        // Posizione [m]
    double vx, vy;      // Velocità [m/s]
    double q;           // Carica [C]
    double m;           // Massa [kg]
    bool active;        // Particella attiva
    
    // Costruttore di default necessario per resize()
    Particle() : x(0.0), y(0.0), vx(0.0), vy(0.0), q(0.0), m(0.0), active(false) {}
    
    // Costruttore parametrizzato
    Particle(double x0, double y0, double vx0, double vy0, double charge, double mass) 
        : x(x0), y(y0), vx(vx0), vy(vy0), q(charge), m(mass), active(true) {}
};

// Classe principale PIC con parallelizzazione
class PIC_DLA_Simulator {
private:
    // Campi elettromagnetici sulla griglia
    std::vector<std::vector<double>> Ex, Ey, Ez;    // Campo elettrico
    std::vector<std::vector<double>> Bx, By, Bz;    // Campo magnetico
    std::vector<std::vector<double>> Jx, Jy, Jz;    // Densità di corrente
    std::vector<std::vector<double>> rho;            // Densità di carica
    std::vector<std::vector<double>> eps_r;          // Permittività relativa
    
    // Particelle
    std::vector<Particle> particles;
    
    // Parametri temporali
    double time;
    int time_step;
    
    // Numero di thread per OpenMP
    int num_threads;
    
public:
    PIC_DLA_Simulator() : time(0.0), time_step(0) {
        // Ottieni il numero di core disponibili
        num_threads = omp_get_max_threads();
        std::cout << "Usando " << num_threads << " thread per la parallelizzazione" << std::endl;
        
        // Inizializza i campi
        initializeFields();
        // Setup geometria dielettrica
        setupDielectricGeometry();
        // Inizializza particelle
        initializeParticles();
    }
    
    void initializeFields() {
        // Ridimensiona tutti i campi
        Ex.resize(Nx, std::vector<double>(Ny, 0.0));
        Ey.resize(Nx, std::vector<double>(Ny, 0.0));
        Ez.resize(Nx, std::vector<double>(Ny, 0.0));
        
        Bx.resize(Nx, std::vector<double>(Ny, 0.0));
        By.resize(Nx, std::vector<double>(Ny, 0.0));
        Bz.resize(Nx, std::vector<double>(Ny, 0.0));
        
        Jx.resize(Nx, std::vector<double>(Ny, 0.0));
        Jy.resize(Nx, std::vector<double>(Ny, 0.0));
        Jz.resize(Nx, std::vector<double>(Ny, 0.0));
        
        rho.resize(Nx, std::vector<double>(Ny, 0.0));
        eps_r.resize(Nx, std::vector<double>(Ny, eps_vacuum));
    }
    
    void setupDielectricGeometry() {
        // Parallelizzazione della setup geometrica
        double period = 29e-6;         // Distanza centro-centro [m] (29 µm)
        double tooth_width = 14.5e-6;  // Larghezza dente [m] (50% del periodo)
        double tooth_height = 10e-6;   // Altezza dente [m]
        
        // Parallelizza il loop sulla griglia
        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                double x = i * dx;
                double y = j * dy;
                
                // Gap di accelerazione al centro
                double gap_center = H / 2.0;
                double gap_half_height = 5e-6;
                
                // Determina se siamo in un dente
                double x_in_period = fmod(x, period);
                bool in_tooth = (x_in_period < tooth_width);
                
                // Denti superiori e inferiori
                bool upper_tooth = (y > gap_center + gap_half_height) && 
                                  (y < gap_center + gap_half_height + tooth_height) && in_tooth;
                bool lower_tooth = (y < gap_center - gap_half_height) && 
                                  (y > gap_center - gap_half_height - tooth_height) && in_tooth;
                
                if (upper_tooth || lower_tooth) {
                    eps_r[i][j] = eps_sio2;
                } else {
                    eps_r[i][j] = eps_vacuum;
                }
            }
        }
    }
    
    void initializeParticles() {
        // Crea un fascio di protoni
        int n_particles = 10000;  // Aumentato per testare la parallelizzazione
        double beam_energy = 100e3 * e; // 100 keV in Joule
        double beam_x_start = 10e-6;    // Posizione iniziale x
        double beam_y_center = H / 2.0; // Centro del gap
        double beam_spread_y = 2e-6;    // Spread del fascio
        
        // Parallelizzazione della generazione particelle
        particles.resize(n_particles);
        
        #pragma omp parallel
        {
            // Ogni thread ha il proprio generatore random
            std::random_device rd;
            std::mt19937 gen(rd() + omp_get_thread_num()); // Seed diverso per ogni thread
            std::normal_distribution<double> dist_y(beam_y_center, beam_spread_y / 3.0);
            std::uniform_real_distribution<double> dist_x(beam_x_start - 1e-6, beam_x_start + 1e-6);
            
            #pragma omp for schedule(static)
            for (int i = 0; i < n_particles; i++) {
                double x0 = dist_x(gen);
                double y0 = dist_y(gen);
                
                // Velocità iniziale dall'energia cinetica
                double v0 = sqrt(2.0 * beam_energy / mp);
                
                particles[i] = Particle(x0, y0, v0, 0.0, e, mp);
            }
        }
        
        std::cout << "Inizializzate " << particles.size() << " particelle" << std::endl;
        std::cout << "Velocità iniziale: " << particles[0].vx / 1e6 << " Mm/s" << std::endl;
    }
    
    void updateLaserField() {
        // Impulso laser Gaussiano
        double t_center = 3.0 * pulse_duration;
        double gaussian = exp(-pow((time - t_center) / pulse_duration, 2));
        
        // Parallelizza il calcolo del campo laser
        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                double x = i * dx;
                double y = j * dy;
                
                // Solo nel gap di accelerazione
                double gap_center = H / 2.0;
                if (abs(y - gap_center) < 5e-6) {
                    double phase = omega0 * (time - x / c);
                    
                    // Campo elettrico nel gap (direzione y per accelerazione)
                    if (eps_r[i][j] == eps_vacuum) {
                        Ey[i][j] += E0_max * gaussian * sin(phase);
                    }
                }
            }
        }
    }
    
    void calculateCurrentDensity() {
        // Reset correnti - parallelizzato
        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                Jx[i][j] = 0.0;
                Jy[i][j] = 0.0;
                Jz[i][j] = 0.0;
                rho[i][j] = 0.0;
            }
        }
        
        // Deposizione delle particelle - con riduzione per evitare race conditions
        const int n_particles = particles.size();
        
        // Usa array temporanei per ogni thread per evitare race conditions
        std::vector<std::vector<std::vector<double>>> thread_Jx(num_threads, 
            std::vector<std::vector<double>>(Nx, std::vector<double>(Ny, 0.0)));
        std::vector<std::vector<std::vector<double>>> thread_Jy(num_threads, 
            std::vector<std::vector<double>>(Nx, std::vector<double>(Ny, 0.0)));
        std::vector<std::vector<std::vector<double>>> thread_rho(num_threads, 
            std::vector<std::vector<double>>(Nx, std::vector<double>(Ny, 0.0)));
        
        #pragma omp parallel
        {
            int thread_id = omp_get_thread_num();
            
            #pragma omp for schedule(static)
            for (int p = 0; p < n_particles; p++) {
                const auto& particle = particles[p];
                if (!particle.active) continue;
                
                int i = static_cast<int>(particle.x / dx);
                int j = static_cast<int>(particle.y / dy);
                
                if (i >= 0 && i < Nx && j >= 0 && j < Ny) {
                    double weight = 1.0 / (dx * dy); // Peso per densità
                    
                    thread_Jx[thread_id][i][j] += particle.q * particle.vx * weight;
                    thread_Jy[thread_id][i][j] += particle.q * particle.vy * weight;
                    thread_rho[thread_id][i][j] += particle.q * weight;
                }
            }
        }
        
        // Somma i contributi di tutti i thread
        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                for (int t = 0; t < num_threads; t++) {
                    Jx[i][j] += thread_Jx[t][i][j];
                    Jy[i][j] += thread_Jy[t][i][j];
                    rho[i][j] += thread_rho[t][i][j];
                }
            }
        }
    }
    
    void updateElectricField() {
        // Parallelizza l'aggiornamento del campo elettrico
        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 1; i < Nx - 1; i++) {
            for (int j = 1; j < Ny - 1; j++) {
                double eps_total = eps_r[i][j] * eps0;
                
                // ∂Ex/∂t = (1/ε)[(∂Bz/∂y) - μ₀Jx]
                double curl_B_x = (Bz[i][j+1] - Bz[i][j-1]) / (2.0 * dy);
                Ex[i][j] += dt * (curl_B_x - mu0 * Jx[i][j]) / eps_total;
                
                // ∂Ey/∂t = (1/ε)[-(∂Bz/∂x) - μ₀Jy]
                double curl_B_y = -(Bz[i+1][j] - Bz[i-1][j]) / (2.0 * dx);
                Ey[i][j] += dt * (curl_B_y - mu0 * Jy[i][j]) / eps_total;
            }
        }
    }
    
    void updateMagneticField() {
        // Parallelizza l'aggiornamento del campo magnetico
        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 1; i < Nx - 1; i++) {
            for (int j = 1; j < Ny - 1; j++) {
                double curl_E = (Ey[i+1][j] - Ey[i-1][j]) / (2.0 * dx) - 
                               (Ex[i][j+1] - Ex[i][j-1]) / (2.0 * dy);
                
                Bz[i][j] -= dt * curl_E;
            }
        }
    }
    
    void interpolateFieldsToParticles(const Particle& p, double& Ex_p, double& Ey_p, double& Bz_p) {
        // Interpolazione bilineare - questa funzione è chiamata per ogni particella
        int i = static_cast<int>(p.x / dx);
        int j = static_cast<int>(p.y / dy);
        
        if (i < 0 || i >= Nx - 1 || j < 0 || j >= Ny - 1) {
            Ex_p = Ey_p = Bz_p = 0.0;
            return;
        }
        
        double wx = (p.x - i * dx) / dx;
        double wy = (p.y - j * dy) / dy;
        
        Ex_p = (1-wx) * (1-wy) * Ex[i][j] + wx * (1-wy) * Ex[i+1][j] +
               (1-wx) * wy * Ex[i][j+1] + wx * wy * Ex[i+1][j+1];
               
        Ey_p = (1-wx) * (1-wy) * Ey[i][j] + wx * (1-wy) * Ey[i+1][j] +
               (1-wx) * wy * Ey[i][j+1] + wx * wy * Ey[i+1][j+1];
               
        Bz_p = (1-wx) * (1-wy) * Bz[i][j] + wx * (1-wy) * Bz[i+1][j] +
               (1-wx) * wy * Bz[i][j+1] + wx * wy * Bz[i+1][j+1];
    }
    
    void updateParticles() {
        // Parallelizza l'aggiornamento delle particelle
        const int n_particles = particles.size();
        
        #pragma omp parallel for schedule(dynamic, 100) // Dynamic scheduling per bilanciare il carico
        for (int p = 0; p < n_particles; p++) {
            auto& particle = particles[p];
            if (!particle.active) continue;
            
            // Interpolazione campi alla posizione della particella
            double Ex_p, Ey_p, Bz_p;
            interpolateFieldsToParticles(particle, Ex_p, Ey_p, Bz_p);
            
            // Forza di Lorentz: F = q(E + v×B)
            double Fx = particle.q * (Ex_p + particle.vy * Bz_p);
            double Fy = particle.q * (Ey_p - particle.vx * Bz_p);
            
            // Aggiornamento velocità (leap-frog)
            particle.vx += dt * Fx / particle.m;
            particle.vy += dt * Fy / particle.m;
            
            // Aggiornamento posizione
            particle.x += dt * particle.vx;
            particle.y += dt * particle.vy;
            
            // Controllo boundaries
            if (particle.x < 0 || particle.x > L || particle.y < 0 || particle.y > H) {
                particle.active = false;
                continue;
            }
            
            // Controllo collisione con dielettrico
            int i = static_cast<int>(particle.x / dx);
            int j = static_cast<int>(particle.y / dy);
            if (i >= 0 && i < Nx && j >= 0 && j < Ny) {
                if (eps_r[i][j] > eps_vacuum) {
                    particle.active = false; // Particella assorbita
                }
            }
        }
    }
    
    void applyBoundaryConditions() {
        // Parallelizza le condizioni al contorno
        #pragma omp parallel sections
        {
            #pragma omp section
            {
                // Pareti conduttrici sui bordi y
                #pragma omp parallel for schedule(static)
                for (int i = 0; i < Nx; i++) {
                    Ex[i][0] = Ex[i][Ny-1] = 0.0;
                    Ey[i][0] = Ey[i][Ny-1] = 0.0;
                }
            }
            
            #pragma omp section
            {
                // Condizioni periodiche o assorbenti sui bordi x
                #pragma omp parallel for schedule(static)
                for (int j = 0; j < Ny; j++) {
                    Ex[0][j] = Ex[1][j];
                    Ex[Nx-1][j] = Ex[Nx-2][j];
                    Ey[0][j] = Ey[1][j];
                    Ey[Nx-1][j] = Ey[Nx-2][j];
                }
            }
        }
    }
    
    void timeStep() {
        // Misura il tempo di esecuzione di ogni step per il profiling
        double start_time = omp_get_wtime();
        
        // 1. Aggiorna campo laser
        updateLaserField();
        
        // 2. Calcola densità di corrente dalle particelle
        calculateCurrentDensity();
        
        // 3. Aggiorna campi elettromagnetici (FDTD)
        updateMagneticField();
        updateElectricField();
        
        // 4. Applica condizioni al contorno
        applyBoundaryConditions();
        
        // 5. Aggiorna particelle
        updateParticles();
        
        // Avanza tempo
        time += dt;
        time_step++;
        
        // Profiling opzionale ogni 1000 step
        if (time_step % 1000 == 0) {
            double elapsed = omp_get_wtime() - start_time;
            std::cout << "Step " << time_step << " time: " << elapsed * 1000.0 << " ms" << std::endl;
        }
    }
    
    void saveResults(const std::string& filename_prefix) {
        // Salva campi elettromagnetici - parallelizza la scrittura se necessario
        std::ofstream field_file(filename_prefix + "_fields.csv");
        field_file << "x,y,Ex,Ey,Bz,eps_r\n";
        field_file << std::fixed << std::setprecision(6);
        
        // Per file di output grandi, si può parallelizzare la preparazione dei dati
        std::vector<std::string> field_lines(Nx * Ny);
        
        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                double x = i * dx * 1e6; // in µm
                double y = j * dy * 1e6; // in µm
                
                std::ostringstream oss;
                oss << std::fixed << std::setprecision(6);
                oss << x << "," << y << ","
                    << Ex[i][j] << "," << Ey[i][j] << ","
                    << Bz[i][j] << "," << eps_r[i][j] << "\n";
                
                field_lines[i * Ny + j] = oss.str();
            }
        }
        
        // Scrittura sequenziale (I/O non può essere parallelizzato efficacemente)
        for (const auto& line : field_lines) {
            field_file << line;
        }
        field_file.close();
        
        // Salva particelle con preparazione parallelizzata
        std::ofstream particle_file(filename_prefix + "_particles.csv");
        particle_file << "x_um,y_um,vx_mm_per_s,vy_mm_per_s,energy_keV,active\n";
        particle_file << std::fixed << std::setprecision(6);
        
        std::vector<std::string> particle_lines(particles.size());
        
        #pragma omp parallel for schedule(static)
        for (size_t i = 0; i < particles.size(); i++) {
            const auto& p = particles[i];
            double kinetic_energy = 0.5 * p.m * (p.vx * p.vx + p.vy * p.vy) / e / 1000.0; // keV
            
            std::ostringstream oss;
            oss << std::fixed << std::setprecision(6);
            oss << p.x * 1e6 << "," << p.y * 1e6 << ","
                << p.vx / 1e6 << "," << p.vy / 1e6 << ","
                << kinetic_energy << "," << (p.active ? 1 : 0) << "\n";
            
            particle_lines[i] = oss.str();
        }
        
        for (const auto& line : particle_lines) {
            particle_file << line;
        }
        particle_file.close();
        
        std::cout << "Risultati salvati: " << filename_prefix << "_*.csv" << std::endl;
    }
    
    int countActiveParticles() {
        // Parallelizza il conteggio delle particelle attive
        int active_count = 0;
        
        #pragma omp parallel for reduction(+:active_count) schedule(static)
        for (size_t i = 0; i < particles.size(); i++) {
            if (particles[i].active) {
                active_count++;
            }
        }
        
        return active_count;
    }
    
    void printStatus() {
        // Parallelizza il calcolo delle statistiche
        int active_particles = 0;
        double avg_energy = 0.0;
        double max_energy = 0.0;
        
        #pragma omp parallel for reduction(+:active_particles,avg_energy) reduction(max:max_energy) schedule(static)
        for (size_t i = 0; i < particles.size(); i++) {
            const auto& p = particles[i];
            if (p.active) {
                active_particles++;
                double energy = 0.5 * p.m * (p.vx * p.vx + p.vy * p.vy) / e / 1000.0; // keV
                avg_energy += energy;
                if (energy > max_energy) {
                    max_energy = energy;
                }
            }
        }
        
        if (active_particles > 0) {
            avg_energy /= active_particles;
        }
        
        std::cout << "Time: " << time * 1e12 << " ps, "
                  << "Active particles: " << active_particles << ", "
                  << "Avg energy: " << avg_energy << " keV, "
                  << "Max energy: " << max_energy << " keV" << std::endl;
    }
    
    void run(double total_time) {
        int total_steps = static_cast<int>(total_time / dt);
        int output_interval = total_steps / 10;
        
        std::cout << "Iniziando simulazione PIC DLA parallelizzata..." << std::endl;
        std::cout << "Parametri:" << std::endl;
        std::cout << "- Thread OpenMP: " << num_threads << std::endl;
        std::cout << "- Griglia: " << Nx << "x" << Ny << std::endl;
        std::cout << "- dt: " << dt * 1e15 << " fs" << std::endl;
        std::cout << "- Tempo totale: " << total_time * 1e12 << " ps" << std::endl;
        std::cout << "- Particelle iniziali: " << particles.size() << std::endl;
        
        double simulation_start = omp_get_wtime();
        
        for (int step = 0; step < total_steps; step++) {
            timeStep();
            
            // Controlla se ci sono ancora particelle attive
            int active_particles = countActiveParticles();
            
            if (step % output_interval == 0) {
                printStatus();
                saveResults("pic_dla_step_" + std::to_string(step));
            }
            
            // Se non ci sono più particelle attive, ferma la simulazione
            if (active_particles == 0) {
                std::cout << "\n*** SIMULAZIONE TERMINATA ANTICIPATAMENTE ***" << std::endl;
                std::cout << "Tutte le particelle sono diventate inattive al passo " << step << std::endl;
                std::cout << "Tempo simulato: " << time * 1e12 << " ps" << std::endl;
                std::cout << "Progresso completato: " << (100.0 * step / total_steps) << "%" << std::endl;
                break;
            }
        }
        
        double simulation_end = omp_get_wtime();
        double total_simulation_time = simulation_end - simulation_start;
        
        std::cout << "\nSimulazione completata!" << std::endl;
        std::cout << "Tempo totale di esecuzione: " << total_simulation_time << " secondi" << std::endl;
        std::cout << "Performance: " << (time_step / total_simulation_time) << " step/secondo" << std::endl;
        
        saveResults("pic_dla_final");
        printStatus();
    }
};

int main() {
    std::cout << "=== Simulazione PIC per Acceleratore DLA (Versione Parallelizzata) ===" << std::endl;
    
    // Imposta il numero di thread OpenMP se desiderato
    // omp_set_num_threads(8); // Decommentare per limitare il numero di thread
    
    try {
        PIC_DLA_Simulator simulator;
        
        // Simula per 200 ps
        double simulation_time = 200e-12; // 200 picosecondi
        simulator.run(simulation_time);
        
    } catch (const std::exception& e) {
        std::cerr << "Errore durante la simulazione: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}