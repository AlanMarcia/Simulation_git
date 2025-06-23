#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <map>
#include <cmath>
#include <numeric>

#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TApplication.h"
#include "TObjArray.h"
#include "TList.h"
#include "TPad.h"
#include "TLatex.h"
#include "TColor.h"
#include "TMultiGraph.h"
#include "TArrow.h"
#include "TLine.h"

// Function to load 1D coordinates from a CSV file (values expected in µm)
std::vector<double> load_coordinates(const std::string& filepath) {
    std::vector<double> coords;
    std::ifstream file(filepath);
    std::string line;

    if (!file.is_open()) {
        std::cerr << "Error: Could not open coordinates file: " << filepath << std::endl;
        return coords;
    }

    while (std::getline(file, line)) {
        if (line.empty() || line.find_first_not_of(" \t\n\v\f\r,") == std::string::npos) continue;
        
        std::stringstream ss(line);
        std::string token;
        
        while (std::getline(ss, token, ',')) {
            try {
                coords.push_back(std::stod(token));
            } catch (const std::exception& e) {
                std::cerr << "Warning: Could not convert '" << token << "' to double in coordinates file." << std::endl;
            }
        }
    }
    file.close();
    return coords;
}

// Function to load 2D data from a CSV file, returns a vector of vectors
std::vector<std::vector<double>> load_2d_csv(const std::string& filename, int& nRows, int& nCols) {
    std::vector<std::vector<double>> data;
    std::ifstream file(filename);
    std::string line;

    if (!file.is_open()) {
        std::cerr << "Error: Could not open 2D CSV file: " << filename << std::endl;
        nRows = 0; nCols = 0;
        return data;
    }

    while (std::getline(file, line)) {
        if (line.empty() || line.find_first_not_of(" \t\n\v\f\r,") == std::string::npos) continue;
        
        std::vector<double> row;
        std::stringstream ss(line);
        std::string token;
        
        while (std::getline(ss, token, ',')) {
            try {
                row.push_back(std::stod(token));
            } catch (const std::exception& e) {
                std::cerr << "Warning: Could not convert '" << token << "' to double." << std::endl;
                row.push_back(0.0);
            }
        }
        
        if (!row.empty()) {
            data.push_back(row);
        }
    }
    file.close();
    
    nRows = data.size();
    nCols = (nRows > 0) ? data[0].size() : 0;
    
    return data;
}

// Function to load geometry parameters from CSV
std::map<std::string, double> load_geometry_params(const std::string& filepath) {
    std::map<std::string, double> params;
    std::ifstream file(filepath);
    std::string line;

    if (!file.is_open()) {
        std::cerr << "Warning: Could not open geometry parameters file: " << filepath << std::endl;
        return params;
    }

    while (std::getline(file, line)) {
        if (line.empty()) continue;
        
        std::stringstream ss(line);
        std::string key, value_str;
        
        if (std::getline(ss, key, ',') && std::getline(ss, value_str)) {
            try {
                params[key] = std::stod(value_str);
            } catch (const std::exception& e) {
                std::cerr << "Warning: Could not convert parameter '" << key << "' value to double." << std::endl;
            }
        }
    }
    file.close();
    return params;
}

// Function to create contour plots with structure outlines
void draw_detailed_outlines(TCanvas* canvas, TH2D* histogram, const std::vector<std::vector<double>>& eps_r_data, 
                           const std::vector<double>& x_coords, const std::vector<double>& y_coords, 
                           double threshold, int color, int linestyle) {
    if (threshold < 0) return; // Skip if no valid threshold
    
    // Create a histogram for contour generation
    int Nx = x_coords.size();
    int Ny = y_coords.size();
    
    if (Nx == 0 || Ny == 0 || eps_r_data.empty()) return;
    
    // Use unique name to avoid memory leak warnings
    static int contour_counter = 0;
    std::string hist_name = "h_eps_contour_" + std::to_string(contour_counter++);
    TH2D* h_eps_contour = new TH2D(hist_name.c_str(), "Permittivity for Contour",
                                   Nx-1, x_coords.data(), Ny-1, y_coords.data());
    
    // Fill histogram with permittivity data
    for (int i = 0; i < Nx-1; ++i) {
        for (int j = 0; j < Ny-1; ++j) {
            if (j < static_cast<int>(eps_r_data.size()) && i < static_cast<int>(eps_r_data[j].size())) {
                h_eps_contour->SetBinContent(i + 1, j + 1, eps_r_data[j][i]);
            }
        }
    }
    
    h_eps_contour->SetContour(1, &threshold);
    h_eps_contour->SetLineColor(color);
    h_eps_contour->SetLineStyle(linestyle);
    h_eps_contour->SetLineWidth(1);
    h_eps_contour->Draw("CONT3 SAME");
    
    canvas->Update();
}

void plot_results_main(const std::string& folder_path_str = "") {
    gROOT->SetBatch(kFALSE); // Run with GUI to display plots
    
    std::string output_folder_name;
    if (!folder_path_str.empty()) {
        output_folder_name = folder_path_str;
        if (output_folder_name.back() == '/' || output_folder_name.back() == '\\') {
            output_folder_name.pop_back();
        }
        std::cout << "Using specified folder: " << output_folder_name << std::endl;
    } else {
        output_folder_name = "geometria_Denti_sfasati_profondi_5um"; // Default folder
        std::cout << "No folder specified, using default: " << output_folder_name << std::endl;
    }

    if (!gSystem->OpenDirectory(output_folder_name.c_str())) {
        std::cerr << "Error: The specified folder '" << output_folder_name << "' does not exist or is not a directory." << std::endl;
        std::cerr << "Please ensure the C++ simulation has run and created this folder with CSV files." << std::endl;
        return;
    }

    // Define file paths
    std::string potential_file = output_folder_name + "/potential.csv";
    std::string ex_file = output_folder_name + "/electric_field_x.csv";
    std::string ey_file = output_folder_name + "/electric_field_y.csv";
    std::string eps_r_file = output_folder_name + "/permittivity.csv";
    std::string x_coords_file = output_folder_name + "/x_coordinates.csv";
    std::string y_coords_file = output_folder_name + "/y_coordinates.csv";
    std::string geometry_params_file = output_folder_name + "/geometry_params.csv";

    // Load data
    int nRows_V, nCols_V, nRows_Ex, nCols_Ex, nRows_Ey, nCols_Ey, nRows_eps, nCols_eps;
    
    std::vector<std::vector<double>> V_data = load_2d_csv(potential_file, nRows_V, nCols_V);
    std::vector<std::vector<double>> Ex_data = load_2d_csv(ex_file, nRows_Ex, nCols_Ex);
    std::vector<std::vector<double>> Ey_data = load_2d_csv(ey_file, nRows_Ey, nCols_Ey);
    std::vector<std::vector<double>> eps_r_data = load_2d_csv(eps_r_file, nRows_eps, nCols_eps);
    std::vector<double> x_coords = load_coordinates(x_coords_file);
    std::vector<double> y_coords = load_coordinates(y_coords_file);
    std::map<std::string, double> geo_params = load_geometry_params(geometry_params_file);

    if (V_data.empty() || Ex_data.empty() || Ey_data.empty() || eps_r_data.empty() || 
        x_coords.empty() || y_coords.empty()) {
        std::cerr << "One or more data files could not be loaded. Aborting plot." << std::endl;
        return;
    }

    int Nx = x_coords.size();
    int Ny = y_coords.size();

    std::cout << "Data loaded successfully. Nx=" << Nx << ", Ny=" << Ny << std::endl;

    // Calculate E-field magnitude
    std::vector<std::vector<double>> E_mag_data(nRows_Ex, std::vector<double>(nCols_Ex));
    for (int i = 0; i < nRows_Ex; ++i) {
        for (int j = 0; j < nCols_Ex; ++j) {
            double ex = Ex_data[i][j];
            double ey = Ey_data[i][j];
            E_mag_data[i][j] = std::sqrt(ex * ex + ey * ey);
        }
    }

    // Calculate outline threshold dynamically
    double outline_threshold = -1.0;
    if (!eps_r_data.empty()) {
        double eps_min = eps_r_data[0][0];
        double eps_max = eps_r_data[0][0];
        
        for (const auto& row : eps_r_data) {
            for (double val : row) {
                if (val < eps_min) eps_min = val;
                if (val > eps_max) eps_max = val;
            }
        }
        
        if (eps_max > eps_min + 0.5) {
            outline_threshold = (eps_min + eps_max) / 2.0;
            std::cout << "Dynamically calculated outline threshold: " << outline_threshold 
                     << " (from eps_r min/max: " << eps_min << ", " << eps_max << ")" << std::endl;
        } else {
            std::cout << "Warning: Permittivity values are too close. Could not reliably determine outline threshold." << std::endl;
        }
    }

    // Determine y-center for profile plot
    int y_center_gap_idx = -1;
    double y_center_gap_abs = -1.0;
    
    if (!geo_params.empty()) {
        // Try toothed geometry logic first
        auto it_base = geo_params.find("y_si_base_height");
        auto it_teeth = geo_params.find("y_teeth_height");
        auto it_gap = geo_params.find("y_vacuum_gap_thick");
        
        if (it_base != geo_params.end() && it_teeth != geo_params.end() && it_gap != geo_params.end()) {
            y_center_gap_abs = it_base->second + it_teeth->second + (it_gap->second / 2.0);
            std::cout << "Using toothed geometry logic for profile plot y-center." << std::endl;
        } else {
            // Try flat geometry logic
            auto it_layer = geo_params.find("y_si_layer_thick");
            if (it_layer != geo_params.end() && it_gap != geo_params.end()) {
                y_center_gap_abs = it_layer->second + (it_gap->second / 2.0);
                std::cout << "Using flat geometry logic for profile plot y-center." << std::endl;
            }
        }
        
        if (y_center_gap_abs > 0 && !y_coords.empty()) {
            // Find closest y-index
            double min_diff = std::abs(y_coords[0] - y_center_gap_abs);
            y_center_gap_idx = 0;
            for (int i = 1; i < static_cast<int>(y_coords.size()); ++i) {
                double diff = std::abs(y_coords[i] - y_center_gap_abs);
                if (diff < min_diff) {
                    min_diff = diff;
                    y_center_gap_idx = i;
                }
            }
            std::cout << "Calculated y-center for profile plot: " << y_center_gap_abs 
                     << " µm (index: " << y_center_gap_idx << ")" << std::endl;
        }
    }

    // Plot 1: Electric Potential
    TCanvas* c_potential = new TCanvas("c_potential", "Electric Potential", 1000, 750);
    
    TH2D* h_potential = new TH2D("h_potential", "Electric Potential (V);x (#mum);y (#mum)",
                                Nx-1, x_coords.data(), Ny-1, y_coords.data());
    
    for (int i = 0; i < Nx-1; ++i) {
        for (int j = 0; j < Ny-1; ++j) {
            if (j < nRows_V && i < nCols_V) {
                h_potential->SetBinContent(i + 1, j + 1, V_data[j][i]);
            }
        }
    }
      h_potential->Draw("COLZ");
    if (outline_threshold > 0) {
        draw_detailed_outlines(c_potential, h_potential, eps_r_data, x_coords, y_coords, 
                              outline_threshold, kWhite, 2);
    }
    
    c_potential->Update();
    c_potential->Show();
    c_potential->WaitPrimitive(); // Wait for user interaction before continuing
    c_potential->SaveAs((output_folder_name + "/potential_plot_root.png").c_str());
    c_potential->SaveAs((output_folder_name + "/potential_plot_root.pdf").c_str());
    std::cout << "Potential plot saved." << std::endl;

    // Plot 2: Electric Field Magnitude
    TCanvas* c_efield = new TCanvas("c_efield", "Electric Field Magnitude", 1000, 750);
    
    TH2D* h_efield = new TH2D("h_efield", "Electric Field Magnitude |E| (V/#mum);x (#mum);y (#mum)",
                             Nx-1, x_coords.data(), Ny-1, y_coords.data());
    
    for (int i = 0; i < Nx-1; ++i) {
        for (int j = 0; j < Ny-1; ++j) {
            if (j < static_cast<int>(E_mag_data.size()) && i < static_cast<int>(E_mag_data[j].size())) {
                h_efield->SetBinContent(i + 1, j + 1, E_mag_data[j][i]);
            }
        }
    }
      h_efield->Draw("COLZ");
    if (outline_threshold > 0) {
        draw_detailed_outlines(c_efield, h_efield, eps_r_data, x_coords, y_coords, 
                              outline_threshold, kWhite, 2);
    }
    
    c_efield->Update();
    c_efield->Show();
    c_efield->WaitPrimitive(); // Wait for user interaction before continuing
    c_efield->SaveAs((output_folder_name + "/efield_magnitude_plot_root.png").c_str());
    c_efield->SaveAs((output_folder_name + "/efield_magnitude_plot_root.pdf").c_str());
    std::cout << "E-field magnitude plot saved." << std::endl;

    // Plot 3: Permittivity Map
    TCanvas* c_permittivity = new TCanvas("c_permittivity", "Permittivity Map", 1000, 750);
    
    TH2D* h_permittivity = new TH2D("h_permittivity", "Relative Permittivity Map (#epsilon_{r});x (#mum);y (#mum)",
                                   Nx-1, x_coords.data(), Ny-1, y_coords.data());
    
    for (int i = 0; i < Nx-1; ++i) {
        for (int j = 0; j < Ny-1; ++j) {
            if (j < nRows_eps && i < nCols_eps) {
                h_permittivity->SetBinContent(i + 1, j + 1, eps_r_data[j][i]);
            }
        }
    }
      h_permittivity->Draw("COLZ");
    if (outline_threshold > 0) {
        draw_detailed_outlines(c_permittivity, h_permittivity, eps_r_data, x_coords, y_coords, 
                              outline_threshold, kBlack, 2);
    }
    
    c_permittivity->Update();
    c_permittivity->Show();
    c_permittivity->WaitPrimitive(); // Wait for user interaction before continuing
    c_permittivity->SaveAs((output_folder_name + "/permittivity_map_plot_root.png").c_str());
    c_permittivity->SaveAs((output_folder_name + "/permittivity_map_plot_root.pdf").c_str());
    std::cout << "Permittivity map plot saved." << std::endl;

    // Plot 4: Electric Field Vector Plot (Quiver equivalent)
    TCanvas* c_vectors = new TCanvas("c_vectors", "Electric Field Vectors", 1200, 900);
    
    // Create background field magnitude plot
    TH2D* h_background = new TH2D("h_background", "Electric Field Vectors in Vacuum;x (#mum);y (#mum)",
                                 Nx-1, x_coords.data(), Ny-1, y_coords.data());
    
    for (int i = 0; i < Nx-1; ++i) {
        for (int j = 0; j < Ny-1; ++j) {
            if (j < static_cast<int>(E_mag_data.size()) && i < static_cast<int>(E_mag_data[j].size())) {
                h_background->SetBinContent(i + 1, j + 1, E_mag_data[j][i]);
            }
        }
    }
    
    h_background->Draw("COLZ");
    
    // Add field vectors (simplified - every 5th point)
    int skip_rate = 5;
    double arrow_scale = 0.1; // Adjust as needed
    
    for (int i = 0; i < Nx; i += skip_rate) {
        for (int j = 0; j < Ny; j += skip_rate) {
            if (i < static_cast<int>(x_coords.size()) && j < static_cast<int>(y_coords.size()) &&
                j < static_cast<int>(Ex_data.size()) && i < static_cast<int>(Ex_data[j].size())) {
                
                // Check if in vacuum (skip if in material)
                if (outline_threshold > 0 && j < static_cast<int>(eps_r_data.size()) && 
                    i < static_cast<int>(eps_r_data[j].size()) && 
                    eps_r_data[j][i] >= outline_threshold) {
                    continue; // Skip material regions
                }
                
                double x_pos = x_coords[i];
                double y_pos = y_coords[j];
                double ex = Ex_data[j][i];
                double ey = Ey_data[j][i];
                
                double magnitude = std::sqrt(ex * ex + ey * ey);
                if (magnitude > 1e-10) { // Avoid zero-length arrows
                    double scale = arrow_scale / magnitude; // Normalize and scale
                    double dx = ex * scale;
                    double dy = ey * scale;
                    
                    TArrow* arrow = new TArrow(x_pos, y_pos, x_pos + dx, y_pos + dy, 0.01, ">");
                    arrow->SetLineColor(kBlack);
                    arrow->SetLineWidth(1);
                    arrow->Draw();
                }
            }
        }
    }
      if (outline_threshold > 0) {
        draw_detailed_outlines(c_vectors, h_background, eps_r_data, x_coords, y_coords, 
                              outline_threshold, kBlack, 2);
    }
    
    c_vectors->Update();
    c_vectors->Show();
    c_vectors->WaitPrimitive(); // Wait for user interaction before continuing
    c_vectors->SaveAs((output_folder_name + "/efield_quiver_vacuum_plot_root.png").c_str());
    c_vectors->SaveAs((output_folder_name + "/efield_quiver_vacuum_plot_root.pdf").c_str());
    std::cout << "E-field vector plot saved." << std::endl;

    // Plot 5: Profile plot at the center of the vacuum gap
    if (y_center_gap_idx >= 0 && y_center_gap_idx < static_cast<int>(y_coords.size())) {
        TCanvas* c_profile = new TCanvas("c_profile", "Center Gap Profile", 1200, 800);
        
        // Extract profile data
        std::vector<double> V_profile, E_profile;
        for (int i = 0; i < Nx; ++i) {
            if (y_center_gap_idx < static_cast<int>(V_data.size()) && 
                i < static_cast<int>(V_data[y_center_gap_idx].size())) {
                V_profile.push_back(V_data[y_center_gap_idx][i]);
            }
            if (y_center_gap_idx < static_cast<int>(E_mag_data.size()) && 
                i < static_cast<int>(E_mag_data[y_center_gap_idx].size())) {
                E_profile.push_back(E_mag_data[y_center_gap_idx][i]);
            }
        }
        
        if (!V_profile.empty() && !E_profile.empty()) {
            // Create graphs
            TGraph* g_potential = new TGraph(x_coords.size(), x_coords.data(), V_profile.data());
            TGraph* g_efield = new TGraph(x_coords.size(), x_coords.data(), E_profile.data());
            
            // Set up multi-graph for dual y-axis effect
            TMultiGraph* mg = new TMultiGraph();
            
            g_potential->SetLineColor(kBlue);
            g_potential->SetLineWidth(2);
            g_potential->SetTitle("Potential (V)");
            
            g_efield->SetLineColor(kRed);
            g_efield->SetLineWidth(2);
            g_efield->SetLineStyle(2);
            g_efield->SetTitle("|E| (V/#mum)");
            
            mg->Add(g_potential, "L");
            mg->Add(g_efield, "L");
            
            mg->SetTitle(Form("Profile at y = %.2f #mum (Center of Vacuum Gap);x (#mum);Values", 
                             y_coords[y_center_gap_idx]));
            mg->Draw("A");
            
            TLegend* leg_profile = new TLegend(0.7, 0.75, 0.88, 0.88);
            leg_profile->AddEntry(g_potential, "Potential (V)", "l");
            leg_profile->AddEntry(g_efield, "|E| (V/#mum)", "l");
            leg_profile->Draw();
              c_profile->SetGrid();
            c_profile->Update();
            c_profile->Show();
            c_profile->WaitPrimitive(); // Wait for user interaction before continuing
            c_profile->SaveAs((output_folder_name + "/center_gap_profile_plot_root.png").c_str());
            c_profile->SaveAs((output_folder_name + "/center_gap_profile_plot_root.pdf").c_str());
            std::cout << "Profile plot saved." << std::endl;
        }
          // delete c_profile;
    }

    // Cleanup
    delete c_potential;
    delete c_efield;
    delete c_permittivity;
    delete c_vectors;
    
    std::cout << "ROOT plotting finished. All plots saved to " << output_folder_name << std::endl;
    std::cout << "Click on each plot window to proceed to the next plot." << std::endl;
}

// Main guard for script execution
#ifndef __CINT__
int main(int argc, char **argv) {
    std::string folder_path;
    if (argc > 1) {
        folder_path = argv[1];
    }
    plot_results_main(folder_path);
    return 0;
}
#endif

// Wrapper function for ROOT interpreter
void plot_results_root(const std::string& folder_path = "") {
    plot_results_main(folder_path);
}
