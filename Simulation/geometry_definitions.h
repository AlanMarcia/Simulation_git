#ifndef GEOMETRY_DEFINITIONS_H
#define GEOMETRY_DEFINITIONS_H

#include <string>
#include <vector>
#include <iostream> // For std::cerr in helper functions if needed directly

// Enum per i tipi di geometria
enum class GeometryType {
    PIANA,
    DENTI_SFASATI_PROFONDI,
    DENTI_UGUALI, 
    DENTI_SFASATI_PROFONDI_NM, // Aggiunto per la nuova geometria
    UNKNOWN
};

// Struttura per informazioni sul materiale in un punto x per la geometria "denti_sfasati_profondi"
struct PointMaterialInfo {
    bool is_silicon_tooth_region;
    double current_dig_depth_at_x;
    double current_tooth_height_at_x;
};

// Parametri di configurazione comuni della geometria
struct GeometryConfig {
    double h; // µm
    double L_total; // µm
    double x_free_space; // µm
    double x_structure_len; // µm
    double H_total_val; // µm
    double current_tolerance;
    double current_eps_material;
    double eps_vacuum_val; // Aggiunto per passarlo alla configurazione della permittività
};

// Parametri specifici per la geometria "piana"
struct PianaSpecificParams {
    double y_si_layer_thick; // µm
    double y_vacuum_gap_thick; // µm
};

// Parametri specifici per la geometria "denti_sfasati_profondi"
struct DentiSfasatiProfondiSpecificParams {
    double y_si_base_height; // µm
    double initial_y_teeth_height; // µm
    double y_teeth_height_decrement; // µm
    double x_teeth_width; // µm
    double initial_x_spacing_width; // µm
    double x_spacing_width_increment; // µm
    double initial_y_spacing_dig_depth; // µm
    double y_spacing_dig_depth_increment; // µm
    double y_vacuum_gap_thick; // µm (gap nominale sopra i denti)
};
// Parametri specifici per la geometria "denti_sfasati_profondi_nm"
struct DentiSfasatiProfondiNmSpecificParams {
    double y_si_base_height; // µm
    double initial_y_teeth_height; // µm
    double y_teeth_height_decrement; // µm
    double x_teeth_width; // µm
    double initial_x_spacing_width; // µm
    double x_spacing_width_increment; // µm
    double initial_y_spacing_dig_depth; // µm
    double y_spacing_dig_depth_increment; // µm
    double y_vacuum_gap_thick; // µm (gap nominale sopra i denti)
};


// Parametri specifici per la geometria "denti_uguali"
struct DentiUgualiSpecificParams {
    double y_si_base_height;      // µm
    double y_tooth_height;        // µm
    double x_tooth_width;         // µm
    double x_spacing_width;       // µm
    double y_vacuum_gap_thick;    // µm
};


// Funzioni di utilità
GeometryType stringToGeometryType(const std::string& str);
std::string geometryTypeToString(GeometryType type);

// Funzioni di inizializzazione dei parametri
void initializePianaGeometry(GeometryConfig& config, PianaSpecificParams& piana_params, double common_h, double eps_sio2, double eps_vac);
void initializeDentiSfasatiProfondiGeometry(GeometryConfig& config, DentiSfasatiProfondiSpecificParams& denti_params, double common_h, double eps_si, double eps_vac);
void initializeDentiSfasatiProfondiNmGeometry(GeometryConfig& config, DentiSfasatiProfondiSpecificParams& denti_params, double common_h, double eps_si, double eps_vac);
void initializeDentiUgualiGeometry(GeometryConfig& config, DentiUgualiSpecificParams& du_params, double common_h, double eps_si, double eps_vac); // Added

// Funzione helper per la geometria "denti_sfasati_profondi"
PointMaterialInfo get_denti_point_x_info(double x_target_rel, double total_struct_len,
                                         const DentiSfasatiProfondiSpecificParams& params);

// Funzioni di setup della mappa di permittività
void setupPianaPermittivity(std::vector<std::vector<double>>& eps_r,
                            const GeometryConfig& config,
                            const PianaSpecificParams& piana_params,
                            int Nx, int Ny);

void setupDentiSfasatiProfondiPermittivity(std::vector<std::vector<double>>& eps_r,
                                           const GeometryConfig& config,
                                           const DentiSfasatiProfondiSpecificParams& denti_params,
                                           int Nx, int Ny);
                                           
// Added for new geometry
void setupDentiSfasatiProfondiNmPermittivity(std::vector<std::vector<double>>& eps_r,
                                           const GeometryConfig& config,
                                           const DentiSfasatiProfondiSpecificParams& denti_params,
                                           int Nx, int Ny);

void setupDentiUgualiPermittivity(std::vector<std::vector<double>>& eps_r, // Added
                                  const GeometryConfig& config,
                                  const DentiUgualiSpecificParams& du_params,
                                  int Nx, int Ny);

// Funzione di setup delle condizioni al contorno
void setupBoundaryConditions(std::vector<std::vector<double>>& V,
                             std::vector<std::vector<bool>>& fixed_mask,
                             const std::vector<std::vector<double>>& eps_r, // Necessaria per applicare BC solo al silicio
                             const GeometryConfig& config,
                             double V_L, double V_R,
                             int Nx, int Ny);

// Funzione per salvare i parametri della geometria
void saveGeometryParams(const std::string& filename,
                        GeometryType type,
                        const GeometryConfig& config,
                        const PianaSpecificParams* piana_params, // Puntatore, può essere nullptr
                        const DentiSfasatiProfondiSpecificParams* denti_params, // Puntatore, può essere nullptr
                        const DentiUgualiSpecificParams* du_params); // Added, Puntatore, può essere nullptr

#endif // GEOMETRY_DEFINITIONS_H
