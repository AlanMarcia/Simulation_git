#include "geometry_definitions.h"
#include <fstream>
#include <iomanip>
#include <cmath> // Per std::max, std::abs
#include <algorithm> // Per std::min

GeometryType stringToGeometryType(const std::string& str) {
    if (str == "piana") return GeometryType::PIANA;
    if (str == "piana_rastremata") return GeometryType::PIANA_RASTREMATA;
    if (str == "piana_variabile") return GeometryType::PIANA_VARIABILE;
    if (str == "denti_sfasati_profondi") return GeometryType::DENTI_SFASATI_PROFONDI;
    if (str == "denti_uguali") return GeometryType::DENTI_UGUALI; // Added
    return GeometryType::UNKNOWN;
}

std::string geometryTypeToString(GeometryType type) {
    switch (type) {
        case GeometryType::PIANA: return "piana";
        case GeometryType::PIANA_RASTREMATA: return "piana_rastremata";
        case GeometryType::PIANA_VARIABILE: return "piana_variabile";
        case GeometryType::DENTI_SFASATI_PROFONDI: return "denti_sfasati_profondi";
        case GeometryType::DENTI_UGUALI: return "denti_uguali"; // Added
        default: return "unknown";
    }
}

void initializePianaGeometry(GeometryConfig& config, PianaSpecificParams& piana_params, double common_h, double eps_sio2, double eps_vac) {
    config.h = common_h;
    config.L_total = 320.0;
    config.x_free_space = 10.0;
    config.x_structure_len = 300.0;
    config.H_total_val = 30.0;
    config.current_tolerance = 1e-3;
    config.current_eps_material = eps_sio2;
    config.eps_vacuum_val = eps_vac;

    piana_params.y_si_layer_thick = 10.0;
    piana_params.y_vacuum_gap_thick = 10.0;
}

void initializePianaRastremataGeometry(GeometryConfig& config, PianaRastremataSpecificParams& piana_rastremata_params, double common_h, double eps_sio2, double eps_vac) {
    config.h = common_h;
    config.L_total = 320.0;
    config.x_free_space = 10.0;
    config.x_structure_len = 300.0;
    config.H_total_val = 45.0;
    config.current_tolerance = 1e-3;
    config.current_eps_material = eps_sio2;
    config.eps_vacuum_val = eps_vac;

    piana_rastremata_params.y_si_layer_thick_start = 20.0; // Spessore iniziale 20 µm
    piana_rastremata_params.y_si_layer_thick_end = 14.0;    // Spessore finale 14 µm
    piana_rastremata_params.y_vacuum_gap_thick = 5.0;     // Gap di vuoto
}

void initializePianaVariabileGeometry(GeometryConfig& config, PianaVariabileSpecificParams& piana_variabile_params, double common_h, double eps_sio2, double eps_vac) {
    config.h = common_h;
    config.L_total = 320.0;
    config.x_free_space = 10.0;
    config.x_structure_len = 300.0;
    config.H_total_val = 30.0;
    config.current_tolerance = 1e-3;
    config.current_eps_material = eps_sio2;
    config.eps_vacuum_val = eps_vac;

    piana_variabile_params.y_si_layer_thick_right = 10.0; // Spessore a destra 10 µm
    piana_variabile_params.y_si_layer_thick_left = 4.0;   // Spessore a sinistra 4 µm
    piana_variabile_params.y_vacuum_gap_thick = 10.0;     // Gap di vuoto
}

void initializeDentiSfasatiProfondiGeometry(GeometryConfig& config, DentiSfasatiProfondiSpecificParams& denti_params, double common_h, double eps_sio2, double eps_vac) {
    config.h = common_h;
    config.L_total = 320.0;
    config.x_free_space = 10.0;
    config.x_structure_len = 300.0;    config.H_total_val = 50.0; // Altezza maggiore per i denti
    config.current_tolerance = 1e-5;
    config.current_eps_material = eps_sio2;
    config.eps_vacuum_val = eps_vac;

    denti_params.y_si_base_height = 10.0;
    denti_params.initial_y_teeth_height = 10.0;
    denti_params.y_teeth_height_decrement = 1.0;
    denti_params.x_teeth_width = 10.0;
    denti_params.initial_x_spacing_width = 10.0;
    denti_params.x_spacing_width_increment = 2.0;
    denti_params.initial_y_spacing_dig_depth = 0.0;
    denti_params.y_spacing_dig_depth_increment = 1.0;
    denti_params.y_vacuum_gap_thick = 10.0;
}

void initializeDentiUgualiGeometry(GeometryConfig& config, DentiUgualiSpecificParams& du_params, double common_h, double eps_sio2, double eps_vac) {
    config.h = common_h;
    config.L_total = 320.0;
    config.x_free_space = 10.0;
    config.x_structure_len = 300.0;    config.H_total_val = 50.0; // Example H_total for denti uguali
    config.current_tolerance = 1e-5;
    config.current_eps_material = eps_sio2;
    config.eps_vacuum_val = eps_vac;

    du_params.y_si_base_height = 10.0;
    du_params.y_tooth_height = 5.0;
    du_params.x_tooth_width = 10.0;
    du_params.x_spacing_width = 10.0; // Spacing between teeth
    du_params.y_vacuum_gap_thick = 20.0;
}

PointMaterialInfo get_denti_point_x_info(double x_target_rel, double total_struct_len,
                                         const DentiSfasatiProfondiSpecificParams& params) {
    PointMaterialInfo info = {false, 0.0, 0.0};
    double current_x_pos = 0.0;
    double current_space_w = params.initial_x_spacing_width;
    double current_dig_d = params.initial_y_spacing_dig_depth;
    double current_tooth_h = params.initial_y_teeth_height;
    bool segment_is_tooth = true;

    if (x_target_rel < 0 || x_target_rel >= total_struct_len) {
        return info;
    }

    while (current_x_pos < total_struct_len) {
        if (segment_is_tooth) {
            double tooth_end_x = current_x_pos + params.x_teeth_width;
            if (x_target_rel >= current_x_pos && x_target_rel < tooth_end_x) {
                info.is_silicon_tooth_region = true;
                info.current_dig_depth_at_x = 0.0;
                info.current_tooth_height_at_x = current_tooth_h;
                return info;
            }
            current_x_pos = tooth_end_x;
            segment_is_tooth = false;
        } else {
            double space_end_x = current_x_pos + current_space_w;
            if (x_target_rel >= current_x_pos && x_target_rel < space_end_x) {
                info.is_silicon_tooth_region = false;
                info.current_dig_depth_at_x = current_dig_d;
                info.current_tooth_height_at_x = 0.0;
                return info;
            }
            current_x_pos = space_end_x;
            current_space_w += params.x_spacing_width_increment;
            current_dig_d += params.y_spacing_dig_depth_increment;
            if (current_dig_d > params.y_si_base_height) {
                current_dig_d = params.y_si_base_height;
            }
            current_tooth_h = std::max(0.0, current_tooth_h - params.y_teeth_height_decrement);
            segment_is_tooth = true;
        }
    }
    return info; // Dovrebbe essere raggiunto solo se x_target_rel è esattamente total_struct_len
}


void setupPianaPermittivity(std::vector<std::vector<double>>& eps_r,
                            const GeometryConfig& config,
                            const PianaSpecificParams& piana_params,
                            int Nx, int Ny) {
    const int idx_x_struct_start = static_cast<int>(config.x_free_space / config.h);
    const int idx_x_struct_end = static_cast<int>((config.x_free_space + config.x_structure_len) / config.h);
    const int idx_y_si_bot_end = static_cast<int>(piana_params.y_si_layer_thick / config.h);
    const int idx_y_si_top_start = static_cast<int>((piana_params.y_si_layer_thick + piana_params.y_vacuum_gap_thick) / config.h);

    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            eps_r[i][j] = config.eps_vacuum_val; // Default to vacuum
            if (i >= idx_x_struct_start && i <= idx_x_struct_end) {
                // Bottom Layer
                if (j >= 0 && j <= idx_y_si_bot_end) {
                    eps_r[i][j] = config.current_eps_material;
                }
                // Top Layer
                if (j >= idx_y_si_top_start && j < Ny) {
                    eps_r[i][j] = config.current_eps_material;
                }
            }
        }
    }
}

void setupPianaRastremataPermittivity(std::vector<std::vector<double>>& eps_r,
                                      const GeometryConfig& config,
                                      const PianaRastremataSpecificParams& piana_rastremata_params,
                                      int Nx, int Ny) {
    const int idx_x_struct_start = static_cast<int>(config.x_free_space / config.h);
    const int idx_x_struct_end = static_cast<int>((config.x_free_space + config.x_structure_len) / config.h);
    
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            eps_r[i][j] = config.eps_vacuum_val; // Default to vacuum
            
            if (i >= idx_x_struct_start && i <= idx_x_struct_end) {
                // Calcolo della posizione relativa x nella struttura (0 a 1)
                double x_position_rel = static_cast<double>(i - idx_x_struct_start) / static_cast<double>(idx_x_struct_end - idx_x_struct_start);
                
                // Interpolazione lineare dello spessore per ENTRAMBI i layer: da thick_start a thick_end
                double current_bottom_thickness = piana_rastremata_params.y_si_layer_thick_start + 
                    x_position_rel * (piana_rastremata_params.y_si_layer_thick_end - piana_rastremata_params.y_si_layer_thick_start);
                double current_top_thickness = current_bottom_thickness; // Stesso spessore sopra e sotto
                
                // Calcolo delle posizioni y: il canale si allarga automaticamente
                const int idx_y_bottom_end = static_cast<int>(current_bottom_thickness / config.h);
                const int idx_y_top_start = static_cast<int>((config.H_total_val - current_top_thickness) / config.h);
                
                // Bottom Layer (rastremata) - da y=0 a current_bottom_thickness
                if (j >= 0 && j <= idx_y_bottom_end) {
                    eps_r[i][j] = config.current_eps_material;
                }
                // Top Layer (rastremata) - da (H_total - current_top_thickness) a H_total
                if (j >= idx_y_top_start && j < Ny) {
                    eps_r[i][j] = config.current_eps_material;
                }
                // Il vacuum gap (canale) è automaticamente la regione tra idx_y_bottom_end e idx_y_top_start
            }
        }
    }
}

void setupPianaVariabilePermittivity(std::vector<std::vector<double>>& eps_r,
                                     const GeometryConfig& config,
                                     const PianaVariabileSpecificParams& piana_variabile_params,
                                     int Nx, int Ny) {
    const int idx_x_struct_start = static_cast<int>(config.x_free_space / config.h);
    const int idx_x_struct_end = static_cast<int>((config.x_free_space + config.x_structure_len) / config.h);
    
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            eps_r[i][j] = config.eps_vacuum_val; // Default to vacuum
            
            if (i >= idx_x_struct_start && i <= idx_x_struct_end) {
                // Calcolo della posizione relativa x nella struttura (0 = sinistra, 1 = destra)
                double x_position_rel = static_cast<double>(i - idx_x_struct_start) / static_cast<double>(idx_x_struct_end - idx_x_struct_start);
                
                // Interpolazione lineare dello spessore da sinistra (4μm) a destra (10μm)
                // x_position_rel = 0 -> sinistra (thick_left = 4μm)
                // x_position_rel = 1 -> destra (thick_right = 10μm)
                double current_bottom_thickness = piana_variabile_params.y_si_layer_thick_left + 
                    x_position_rel * (piana_variabile_params.y_si_layer_thick_right - piana_variabile_params.y_si_layer_thick_left);
                double current_top_thickness = current_bottom_thickness; // Stesso spessore sopra e sotto
                
                // Calcolo delle posizioni y: il canale si allarga automaticamente
                const int idx_y_bottom_end = static_cast<int>(current_bottom_thickness / config.h);
                const int idx_y_top_start = static_cast<int>((config.H_total_val - current_top_thickness) / config.h);
                
                // Bottom Layer (variabile) - da y=0 a current_bottom_thickness
                if (j >= 0 && j <= idx_y_bottom_end) {
                    eps_r[i][j] = config.current_eps_material;
                }
                // Top Layer (variabile) - da (H_total - current_top_thickness) a H_total
                if (j >= idx_y_top_start && j < Ny) {
                    eps_r[i][j] = config.current_eps_material;
                }
                // Il vacuum gap (canale) è automaticamente la regione tra idx_y_bottom_end e idx_y_top_start
            }
        }
    }
}

void setupDentiSfasatiProfondiPermittivity(std::vector<std::vector<double>>& eps_r,
                                           const GeometryConfig& config,
                                           const DentiSfasatiProfondiSpecificParams& denti_params,
                                           int Nx, int Ny) {
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            eps_r[i][j] = config.eps_vacuum_val; // Default to vacuum

            double x_abs = i * config.h;
            double y_abs = j * config.h;

            if (x_abs >= config.x_free_space && x_abs < (config.x_free_space + config.x_structure_len)) {
                double x_coord_in_structure = x_abs - config.x_free_space;
                PointMaterialInfo x_info = get_denti_point_x_info(x_coord_in_structure, config.x_structure_len, denti_params);

                double current_tooth_h = x_info.current_tooth_height_at_x;
                double y_bottom_base_top = denti_params.y_si_base_height;
                double y_actual_bottom_teeth_tip = y_bottom_base_top + (x_info.is_silicon_tooth_region ? current_tooth_h : 0.0);
                double y_actual_top_teeth_tip = config.H_total_val - y_bottom_base_top - (x_info.is_silicon_tooth_region ? current_tooth_h : 0.0);

                // Bottom Base region
                if (y_abs >= 0 && y_abs < y_bottom_base_top) {
                    if (x_info.is_silicon_tooth_region) {
                        eps_r[i][j] = config.current_eps_material;
                    } else {
                        if (y_abs < (y_bottom_base_top - x_info.current_dig_depth_at_x)) {
                            eps_r[i][j] = config.current_eps_material;
                        }
                    }
                }
                // Bottom Teeth region
                else if (x_info.is_silicon_tooth_region && current_tooth_h > 0 &&
                         y_abs >= y_bottom_base_top && y_abs < y_actual_bottom_teeth_tip) {
                    eps_r[i][j] = config.current_eps_material;
                }
                // Top Teeth Region
                else if (x_info.is_silicon_tooth_region && current_tooth_h > 0 &&
                         y_abs >= y_actual_top_teeth_tip && y_abs < (config.H_total_val - y_bottom_base_top)) {
                     eps_r[i][j] = config.current_eps_material;
                }
                // Top Base Region
                else if (y_abs >= (config.H_total_val - y_bottom_base_top) && y_abs < config.H_total_val) {
                    double y_top_base_bottom_surface = config.H_total_val - y_bottom_base_top;
                    if (x_info.is_silicon_tooth_region) {
                        eps_r[i][j] = config.current_eps_material;
                    } else {
                        if (y_abs >= (y_top_base_bottom_surface + x_info.current_dig_depth_at_x)) {
                            eps_r[i][j] = config.current_eps_material;
                        }
                    }
                }
            }
        }
    }
}

void setupDentiUgualiPermittivity(std::vector<std::vector<double>>& eps_r,
                                  const GeometryConfig& config,
                                  const DentiUgualiSpecificParams& du_params,
                                  int Nx, int Ny) {
    double period = du_params.x_tooth_width + du_params.x_spacing_width;

    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            eps_r[i][j] = config.eps_vacuum_val; // Default to vacuum

            double x_abs = i * config.h;
            double y_abs = j * config.h;

            if (x_abs >= config.x_free_space && x_abs < (config.x_free_space + config.x_structure_len)) {
                double x_rel_in_structure = x_abs - config.x_free_space;
                double x_in_period = fmod(x_rel_in_structure, period);

                bool is_in_bottom_tooth_x_range = (x_in_period < du_params.x_tooth_width);
                // For interdigitating, top teeth are shifted by tooth_width + spacing_width / 2, or simply check the other half of period
                // A simpler way for symmetric interdigitating: top teeth are where bottom teeth are not, within the period.
                // However, for standard "denti uguali" often the top teeth align or are offset by a full tooth/space.
                // Let's assume top teeth are offset by half a period for interdigitation, or aligned if not specified.
                // For simplicity here, let's make top teeth align with bottom teeth for "denti uguali"
                // If interdigitated, the logic for x_in_period for top teeth would be different.
                // For this example, let's assume top teeth are directly above bottom teeth spaces if not interdigitated,
                // or directly above bottom teeth if it's a simple stacked structure.
                // The typical "denti uguali" has teeth facing each other across the gap.

                // Bottom silicon layer: base + teeth
                double y_bottom_base_top = du_params.y_si_base_height;
                double y_bottom_tooth_tip = y_bottom_base_top + du_params.y_tooth_height;

                if (y_abs < y_bottom_base_top) { // In bottom base
                    eps_r[i][j] = config.current_eps_material;
                } else if (y_abs < y_bottom_tooth_tip && is_in_bottom_tooth_x_range) { // In bottom tooth
                    eps_r[i][j] = config.current_eps_material;
                }

                // Top silicon layer: base + teeth (inverted)
                double y_top_base_bottom = config.H_total_val - du_params.y_si_base_height;
                double y_top_tooth_root = y_top_base_bottom - du_params.y_tooth_height; // Tip of tooth points downwards

                // For interdigitating, top teeth are often shifted.
                // If top teeth are aligned with bottom spaces:
                bool is_in_top_tooth_x_range = (x_in_period >= du_params.x_tooth_width); // Tooth in the space of bottom layer
                // If top teeth are aligned with bottom teeth (simpler):
                // bool is_in_top_tooth_x_range = is_in_bottom_tooth_x_range;


                if (y_abs >= y_top_base_bottom) { // In top base
                    eps_r[i][j] = config.current_eps_material;
                } else if (y_abs >= y_top_tooth_root && is_in_bottom_tooth_x_range) { // In top tooth (assuming aligned with bottom for simplicity)
                                                                                // Change is_in_bottom_tooth_x_range to a different logic for interdigitated
                    eps_r[i][j] = config.current_eps_material;
                }
            }
        }
    }
}


void setupBoundaryConditions(std::vector<std::vector<double>>& V,
                             std::vector<std::vector<bool>>& fixed_mask,
                             const std::vector<std::vector<double>>& eps_r,
                             const GeometryConfig& config,
                             double V_L, double V_R,
                             int Nx, int Ny) {
    const int actual_idx_x_struct_start = static_cast<int>(config.x_free_space / config.h);
    // Ensure end index is within bounds and represents the last column of the structure
    const int actual_idx_x_struct_end = static_cast<int>((config.x_free_space + config.x_structure_len - config.h * 0.5) / config.h) ; 
                                       // Subtracting h*0.5 ensures we get the index for points up to x_structure_len

    for (int j = 0; j < Ny; ++j) {
        if (actual_idx_x_struct_start >= 0 && actual_idx_x_struct_start < Nx) {
            if (std::abs(eps_r[actual_idx_x_struct_start][j] - config.current_eps_material) < 1e-3) { // Check if it's material
                V[actual_idx_x_struct_start][j] = V_L;
                fixed_mask[actual_idx_x_struct_start][j] = true;
            }
        }
        if (actual_idx_x_struct_end >= 0 && actual_idx_x_struct_end < Nx) {
             if (std::abs(eps_r[actual_idx_x_struct_end][j] - config.current_eps_material) < 1e-3) { // Check if it's material
                V[actual_idx_x_struct_end][j] = V_R;
                fixed_mask[actual_idx_x_struct_end][j] = true;
            }
        }
    }
}

void saveGeometryParams(const std::string& filename,
                        GeometryType type,
                        const GeometryConfig& config,
                        const PianaSpecificParams* piana_params,
                        const PianaRastremataSpecificParams* piana_rastremata_params,
                        const PianaVariabileSpecificParams* piana_variabile_params,
                        const DentiSfasatiProfondiSpecificParams* denti_params,
                        const DentiUgualiSpecificParams* du_params) {
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }
    outfile << std::fixed << std::setprecision(10);
    outfile << "geometry_type," << geometryTypeToString(type) << std::endl;
    outfile << "h," << config.h << std::endl;
    outfile << "x_free_space," << config.x_free_space << std::endl;
    outfile << "x_structure_len," << config.x_structure_len << std::endl;
    outfile << "H_total," << config.H_total_val << std::endl;
    outfile << "eps_material," << config.current_eps_material << std::endl;
    outfile << "eps_vacuum," << config.eps_vacuum_val << std::endl;
    outfile << "tolerance," << config.current_tolerance << std::endl;


    if (type == GeometryType::PIANA && piana_params) {
        outfile << "y_si_layer_thick," << piana_params->y_si_layer_thick << std::endl;
        outfile << "y_vacuum_gap_thick," << piana_params->y_vacuum_gap_thick << std::endl;
    } else if (type == GeometryType::PIANA_RASTREMATA && piana_rastremata_params) {
        outfile << "y_si_layer_thick_start," << piana_rastremata_params->y_si_layer_thick_start << std::endl;
        outfile << "y_si_layer_thick_end," << piana_rastremata_params->y_si_layer_thick_end << std::endl;
        outfile << "y_vacuum_gap_thick," << piana_rastremata_params->y_vacuum_gap_thick << std::endl;
    } else if (type == GeometryType::PIANA_VARIABILE && piana_variabile_params) {
        outfile << "y_si_layer_thick_right," << piana_variabile_params->y_si_layer_thick_right << std::endl;
        outfile << "y_si_layer_thick_left," << piana_variabile_params->y_si_layer_thick_left << std::endl;
        outfile << "y_vacuum_gap_thick," << piana_variabile_params->y_vacuum_gap_thick << std::endl;
    } else if (type == GeometryType::DENTI_SFASATI_PROFONDI && denti_params) {
        outfile << "y_si_base_height," << denti_params->y_si_base_height << std::endl;
        outfile << "initial_y_teeth_height," << denti_params->initial_y_teeth_height << std::endl;
        outfile << "y_teeth_height_decrement," << denti_params->y_teeth_height_decrement << std::endl;
        outfile << "y_vacuum_gap_thick," << denti_params->y_vacuum_gap_thick << std::endl;
        outfile << "x_teeth_width," << denti_params->x_teeth_width << std::endl;
        outfile << "initial_x_spacing_width," << denti_params->initial_x_spacing_width << std::endl;
        outfile << "x_spacing_width_increment," << denti_params->x_spacing_width_increment << std::endl;
        outfile << "initial_y_spacing_dig_depth," << denti_params->initial_y_spacing_dig_depth << std::endl;
        outfile << "y_spacing_dig_depth_increment," << denti_params->y_spacing_dig_depth_increment << std::endl;
    } else if (type == GeometryType::DENTI_UGUALI && du_params) {
        outfile << "y_si_base_height," << du_params->y_si_base_height << std::endl;
        outfile << "y_tooth_height," << du_params->y_tooth_height << std::endl;
        outfile << "x_tooth_width," << du_params->x_tooth_width << std::endl;
        outfile << "x_spacing_width," << du_params->x_spacing_width << std::endl;
        outfile << "y_vacuum_gap_thick," << du_params->y_vacuum_gap_thick << std::endl;
    }
    outfile.close();
    std::cout << "Geometry parameters saved to " << filename << std::endl;
}
