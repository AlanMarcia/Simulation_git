/*
 * Configuration file for batch parallel simulation
 * Modifica questi valori per personalizzare la simulazione
 */

#ifndef SIMULATION_CONFIG_H
#define SIMULATION_CONFIG_H

// Parametri temporali
#define SIM_TIME_STEP          1e-12    // Time step in seconds
#define SIM_MAX_TIME           2e-8     // Maximum simulation time in seconds

// Parametri fisici
#define PROTON_CHARGE          1.602e-19  // Charge of proton in Coulombs
#define PROTON_MASS            1.672e-27  // Mass of proton in kg
#define INITIAL_ENERGY_EV      1000       // Initial kinetic energy in eV

// Parametri della simulazione
#define NUM_PROTONS            10         // Number of protons per simulation

// Posizioni iniziali
#define START_POS_X            -265e-6    // Initial x position in meters
#define START_POS_Y_MIN        -300e-6    // Minimum initial y position in meters
#define START_POS_Y_MAX        300e-6     // Maximum initial y position in meters

// Percorsi dei file
#define FIELD_INPUT_DIR        "export_field"              // Directory with field files
#define OUTPUT_BASE_DIR        "proton_trajectories_output" // Base output directory
#define TRAJECTORIES_FILENAME  "proton_trajectories.csv"    // CSV output filename

// OpenMP settings
#define OMP_SCHEDULE_TYPE      "dynamic"  // "static", "dynamic", or "guided"
#define OMP_SCHEDULE_CHUNK     1          // Chunk size for scheduling

// Debug settings
#define VERBOSE_OUTPUT         1          // 1 for detailed output, 0 for minimal
#define SAVE_FIELD_INFO        0          // 1 to save field grid info, 0 to skip

#endif // SIMULATION_CONFIG_H
