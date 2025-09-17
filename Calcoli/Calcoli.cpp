#include <iostream>
#include <cmath>
#include <vector>   // Utile per salvare i dati per un plot
#include <iomanip>  // Per una stampa più pulita
#include <fstream>  // Per scrivere file CSV

int main() {
    // ---- INPUT ----
    const double initial_K_keV = 100.0;    // Energia Cinetica Iniziale [keV]
    const double G_GV_per_m = 1;         // Gradiente Accelerante [GV/m]
    const double lambda_um = 6.0;          // Lunghezza d'onda Laser [µm]
    const double L_total_um = 100.0;        // Lunghezza massima da simulare [µm]
    const double step_size_um = 0.0001;      // Dimensione del passo di integrazione [µm]

    // ---- COSTANTI FISICHE ----
    const double E0_MeV = 938.272;         // Energia a riposo del protone [MeV]

    // ---- CONVERSIONI E INIZIALIZZAZIONE ----
    const double G_MeV_per_um = G_GV_per_m * 1e-3; // Conversione corretta: 1 GV/m = 0.001 MeV/µm
    const double initial_K_MeV = initial_K_keV / 1000.0;

    double current_K_MeV = initial_K_MeV;
    double total_Delta_Phi = 0.0; // Sfasamento totale accumulato [gradi]

    // Calcoli iniziali
    double initial_E_MeV = initial_K_MeV + E0_MeV;
    const double initial_beta = sqrt(1.0 - pow(E0_MeV / initial_E_MeV, 2));
    const double Lambda_um = initial_beta * lambda_um;

    int num_steps = static_cast<int>(L_total_um / step_size_um);

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Inizio Simulazione..." << std::endl;
    
    // Apri file CSV per salvare i risultati
    std::ofstream csv_file("risultati_stadi.csv");
    csv_file << std::fixed << std::setprecision(6);
    csv_file << "Stadio,Z_inizio_um,Z_fine_um,Lunghezza_um,K_iniziale_keV,K_finale_keV,Beta_iniziale,Beta_finale,Lambda_iniziale_um,Lambda_finale_um,Numero_denti\n";

    // ---- CICLO DI INTEGRAZIONE NUMERICA ----
    for (int i = 0; i < num_steps; ++i) {
        double z_um = i * step_size_um;

        // 1. Calcola la beta corrente all'inizio di questo passo
        double current_E_MeV = current_K_MeV + E0_MeV;
        double current_beta = sqrt(1.0 - pow(E0_MeV / current_E_MeV, 2));

        // 2. Calcola lo sfasamento *in questo piccolo passo* e sommalo al totale
        double dPhi = 360.0 * (step_size_um / Lambda_um - step_size_um / (current_beta * lambda_um));
        total_Delta_Phi += dPhi;

        // Stampa i risultati ogni tanto per non affollare l'output
        if (i % 100 == 0) {
             std::cout << "z = " << z_um << " um, K = " << current_K_MeV * 1000 << " keV, Delta_Phi = " << -total_Delta_Phi << " gradi" << std::endl;
        }

        // 3. Controlla se abbiamo superato i 90 gradi
        if (std::abs(total_Delta_Phi) >= 90.0) {
            std::cout << "\nSTADIO 1:" << std::endl;
            std::cout << "Beta iniziale: " << initial_beta << ", Lambda struttura: " << Lambda_um << " um" << std::endl;
            std::cout << "--------------------------------------------------------" << std::endl;
            std::cout << "Sfasamento di 90 gradi raggiunto a z = " << z_um << " um" << std::endl;
            std::cout << "Energia cinetica finale: " << current_K_MeV * 1000 << " keV" << std::endl;
            std::cout << "Numero di denti: " << (z_um / Lambda_um) << std::endl;
            
            // Salva stadio 1 nel CSV
            double final_E_MeV = current_K_MeV + E0_MeV;
            double final_beta = sqrt(1.0 - pow(E0_MeV / final_E_MeV, 2));
            double final_Lambda_um = final_beta * lambda_um;
            csv_file << "1,0.0," << z_um << "," << z_um << "," << initial_K_MeV * 1000 << "," << current_K_MeV * 1000 
                     << "," << initial_beta << "," << final_beta << "," << Lambda_um << "," << final_Lambda_um 
                     << "," << (z_um / Lambda_um) << "\n";
            
            // --- STADI SUCCESSIVI (2-6) ---
            double stadio_z_um = z_um;
            double stadio_K_MeV = current_K_MeV;
            
            for (int stadio = 2; stadio <= 100; ++stadio) {
                // Calcola la nuova Lambda struttura per questo stadio
                double stadio_E_MeV_iniziale = stadio_K_MeV + E0_MeV;
                double stadio_beta_iniziale = sqrt(1.0 - pow(E0_MeV / stadio_E_MeV_iniziale, 2));
                double stadio_Lambda_um_iniziale = stadio_beta_iniziale * lambda_um;
                
                double stadio_total_Delta_Phi = 0.0;
                int stadio_steps = 0;
                int stadio_print_counter = 0;
                double stadio_z_start = stadio_z_um;
                
                std::cout << "\nSTADIO " << stadio << " (inizia a z = " << stadio_z_start << " um):" << std::endl;
                std::cout << "Beta iniziale stadio: " << stadio_beta_iniziale << ", Lambda struttura iniziale: " << stadio_Lambda_um_iniziale << " um" << std::endl;
                
                while (std::abs(stadio_total_Delta_Phi) < 90.0) {
                    // Aggiorna solo energia, ma usa sempre la Lambda iniziale dello stadio
                    double stadio_E_MeV = stadio_K_MeV + E0_MeV;
                    double stadio_beta = sqrt(1.0 - pow(E0_MeV / stadio_E_MeV, 2));
                    
                    // Usa la Lambda struttura COSTANTE per tutto lo stadio
                    double dPhi_stadio = 360.0 * (step_size_um / stadio_Lambda_um_iniziale - step_size_um / (stadio_beta * lambda_um));
                    stadio_total_Delta_Phi += dPhi_stadio;
                    stadio_z_um += step_size_um;
                    
                    double dK_MeV_stadio = G_MeV_per_um * step_size_um;
                    stadio_K_MeV += dK_MeV_stadio;
                    stadio_steps++;
                    
                    // Stampa ogni 10000 passi per monitorare l'avanzamento
                    if (stadio_steps % 10000 == 0) {
                        stadio_print_counter++;
                        std::cout << "[Stadio " << stadio << "] Step " << stadio_steps << " -> z = " << stadio_z_um << " um, K = " << stadio_K_MeV * 1000 << " keV, Delta_Phi = " << -stadio_total_Delta_Phi << " gradi" << std::endl;
                    }
                }
                
                // Calcola i parametri finali per lo stadio
                double stadio_E_MeV_finale = stadio_K_MeV + E0_MeV;
                double stadio_beta_finale = sqrt(1.0 - pow(E0_MeV / stadio_E_MeV_finale, 2));
                double stadio_Lambda_um_finale = stadio_beta_finale * lambda_um;
                double lunghezza_stadio = stadio_z_um - stadio_z_start;
                
                std::cout << "--------------------------------------------------------" << std::endl;
                std::cout << "Sfasamento di 90 gradi raggiunto a z = " << stadio_z_um << " um" << std::endl;
                std::cout << "Energia cinetica finale: " << stadio_K_MeV * 1000 << " keV" << std::endl;
                std::cout << "Beta finale stadio: " << stadio_beta_finale << ", Lambda struttura finale: " << stadio_Lambda_um_finale << " um" << std::endl;
                std::cout << "Lunghezza stadio: " << lunghezza_stadio << " um" << std::endl;
                std::cout << "Numero di denti: " << (lunghezza_stadio / stadio_Lambda_um_iniziale) << std::endl;
                
                // Salva stadio nel CSV
                double K_iniziale_stadio = (stadio_K_MeV - (G_MeV_per_um * lunghezza_stadio)) * 1000; // keV
                csv_file << stadio << "," << stadio_z_start << "," << stadio_z_um << "," << lunghezza_stadio 
                         << "," << K_iniziale_stadio << "," << stadio_K_MeV * 1000 
                         << "," << stadio_beta_iniziale << "," << stadio_beta_finale 
                         << "," << stadio_Lambda_um_iniziale << "," << stadio_Lambda_um_finale 
                         << "," << (lunghezza_stadio / stadio_Lambda_um_iniziale) << "\n";
            }
            
            break;
        }

        // 4. Aggiorna l'energia per il passo successivo
        double dK_MeV = G_MeV_per_um * step_size_um;
        current_K_MeV += dK_MeV;
    }

    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "Simulazione completata." << std::endl;

    // Chiudi il file CSV
    csv_file.close();

    return 0;
}
