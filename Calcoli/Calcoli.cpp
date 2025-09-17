#include <iostream>
#include <cmath>
#include <vector>   // Utile per salvare i dati per un plot
#include <iomanip>  // Per una stampa più pulita

int main() {
    // ---- INPUT ----
    const double initial_K_keV = 100.0;    // Energia Cinetica Iniziale [keV]
    const double G_GV_per_m = 1;         // Gradiente Accelerante [GV/m]
    const double lambda_um = 2.0;          // Lunghezza d'onda Laser [µm]
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
cd .
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Inizio Simulazione..." << std::endl;
    

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
            std::cout << "\nBeta iniziale: " << initial_beta << ", Lambda struttura: " << Lambda_um << " um" << std::endl;
            std::cout << "--------------------------------------------------------" << std::endl;
            std::cout << "\nSfasamento di 90 gradi raggiunto a z = " << z_um << " um" << std::endl;
            std::cout << "Energia cinetica finale: " << current_K_MeV * 1000 << " keV" << std::endl;
            std::cout << "Numero di denti: " << (z_um / Lambda_um) << std::endl;
            std::cout << "--------------------------------------------------------" << std::endl;
            break;
        }

        // 4. Aggiorna l'energia per il passo successivo
        double dK_MeV = G_MeV_per_um * step_size_um;
        current_K_MeV += dK_MeV;
    }

    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "Simulazione completata." << std::endl;

    return 0;
}
