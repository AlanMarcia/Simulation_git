// filepath: c:\Users\admin\Desktop\Simulation_git\Simulation\plot_trajectories_root.C
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm> // For std::min/max
#include <map>       // For std::map
#include <cmath>     // For std::sqrt, std::stod, std::stoi
#include <numeric>   // For std::accumulate

#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TApplication.h" // To keep plots open
#include "TObjArray.h"
#include "TList.h"
#include "TPad.h"
#include "TLatex.h"     // For text on plots
#include "TColor.h"     // For custom colors if needed
#include "TMultiGraph.h" // For combining multiple graphs

// Physical constants
const double M_PROTON_KG = 1.67262192e-27;  // Mass of proton in kg
const double E_CHARGE_C = 1.60217663e-19;   // Elementary charge in Coulombs

// Function to load 1D coordinates from a CSV file (values expected in µm)
std::vector<double> load_coordinates(const std::string& filepath) {
    std::vector<double> coords;
    std::ifstream file(filepath);
    std::string line;

    if (!file.is_open()) {
        std::cerr << "Error: Coordinates file not found at " << filepath << std::endl;
        return coords; // Return empty vector
    }

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string value_str;
        while(std::getline(ss, value_str, ',')) {
            try {
                coords.push_back(std::stod(value_str));
            } catch (const std::invalid_argument& ia) {
                std::cerr << "Warning: Could not convert value to double in " << filepath << ": " << value_str << " - " << ia.what() << std::endl;
            } catch (const std::out_of_range& oor) {
                std::cerr << "Warning: Value out of range in " << filepath << ": " << value_str << " - " << oor.what() << std::endl;
            }
        }
    }
    file.close();
    return coords;
}

// Function to load 2D data from a CSV file, returns a vector of vectors
// Assumes CSV is saved Ny rows, Nx columns
std::vector<std::vector<double>> load_2d_csv(const std::string& filename, int& nRows, int& nCols) {
    std::vector<std::vector<double>> data;
    std::ifstream file(filename);
    std::string line;
    nRows = 0;
    nCols = 0;

    if (!file.is_open()) {
        std::cerr << "Error: File " << filename << " not found." << std::endl;
        return data; // Return empty vector
    }

    while (std::getline(file, line)) {
        if (line.empty() || line.find_first_not_of(" \\t\\n\\v\\f\\r,") == std::string::npos) { // Skip empty or comma-only lines
            continue;
        }
        nRows++;
        std::vector<double> row_data;
        std::stringstream ss(line);
        std::string value_str;
        int currentCol = 0;
        while(std::getline(ss, value_str, ',')) {
            try {
                row_data.push_back(std::stod(value_str));
                currentCol++;
            } catch (const std::invalid_argument& ia) {
                std::cerr << "Warning: Could not convert value to double in " << filename << " at row " << nRows << ": " << value_str << " - " << ia.what() << std::endl;
            } catch (const std::out_of_range& oor) {
                std::cerr << "Warning: Value out of range in " << filename << " at row " << nRows << ": " << value_str << " - " << oor.what() << std::endl;
            }
        }
        if (nRows == 1) {
            nCols = currentCol;
        } else if (currentCol != nCols && currentCol != 0) {
             std::cerr << "Warning: Inconsistent number of columns in " << filename << " at row " << nRows << ". Expected " << nCols << ", got " << currentCol << std::endl;
        }
        if (!row_data.empty()) {
            data.push_back(row_data);
        }
    }
    file.close();
    nRows = data.size(); // Actual number of non-empty rows read
    if (nRows > 0) {
        nCols = data[0].size(); 
    } else {
        nCols = 0;
    }
    return data;
}

// Structure to hold trajectory data for a single proton
struct ProtonTrajectory {
    int id;
    std::vector<double> time_s;
    std::vector<double> x_m;
    std::vector<double> y_m;
    std::vector<double> vx_m_per_s;
    std::vector<double> vy_m_per_s;
};

// Function to load all proton trajectories from the single CSV file
std::vector<ProtonTrajectory> load_all_trajectories(const std::string& filepath) {
    std::vector<ProtonTrajectory> trajectories_vec;
    std::map<int, ProtonTrajectory> trajectory_map; 

    std::ifstream file(filepath);
    std::string line;

    if (!file.is_open()) {
        std::cerr << "Error: Trajectory file not found: " << filepath << std::endl;
        return trajectories_vec;
    }

    std::getline(file, line); // Skip header line

    int line_num = 1;
    while (std::getline(file, line)) {
        line_num++;
        if (line.empty() || line.find_first_not_of(" \\t\\n\\v\\f\\r,") == std::string::npos) continue;

        std::stringstream ss(line);
        std::string value_str;
        std::vector<std::string> values;

        while(std::getline(ss, value_str, ',')) {
            values.push_back(value_str);
        }

        if (values.size() >= 6) { 
            try {
                int proton_id = std::stoi(values[0]);
                double time_val = std::stod(values[1]);
                double x_val = std::stod(values[2]);
                double y_val = std::stod(values[3]);
                double vx_val = std::stod(values[4]);
                double vy_val = std::stod(values[5]);

                trajectory_map[proton_id].id = proton_id;
                trajectory_map[proton_id].time_s.push_back(time_val);
                trajectory_map[proton_id].x_m.push_back(x_val);
                trajectory_map[proton_id].y_m.push_back(y_val);
                trajectory_map[proton_id].vx_m_per_s.push_back(vx_val);
                trajectory_map[proton_id].vy_m_per_s.push_back(vy_val);

            } catch (const std::exception& e) {
                std::cerr << "Warning: Error parsing line " << line_num << " in " << filepath << ": \"" << line << "\" - " << e.what() << std::endl;
            }
        } else {
             std::cerr << "Warning: Skipping malformed line " << line_num << " in " << filepath << ": \"" << line << "\" (expected at least 6 columns, got " << values.size() << ")" << std::endl;
        }
    }
    file.close();

    for (auto const& pair_item : trajectory_map) {
        trajectories_vec.push_back(pair_item.second);
    }
    return trajectories_vec;
}


void plot_trajectories_main(const std::string& folder_path_str = "") {
    // gROOT->SetBatch(kTRUE); // Run in batch mode (no GUI)

    std::string input_base_folder;
    if (!folder_path_str.empty()) {
        input_base_folder = folder_path_str;
        if (input_base_folder.back() == '/' || input_base_folder.back() == '\\') {
            input_base_folder.pop_back();
        }
        std::cout << "Using specified folder for trajectory data: " << input_base_folder << std::endl;
    } else {
        input_base_folder = "geometria_Denti_uguali_5um"; // Default
        std::cout << "No folder specified, using default for trajectory data: " << input_base_folder << std::endl;
    }

    if (!gSystem->OpenDirectory(input_base_folder.c_str())) {
         std::cerr << "Error: The specified folder '" << input_base_folder << "' does not exist or is not a directory." << std::endl;
         std::cerr << "Please ensure the C++ simulation has run and created this folder with CSV files." << std::endl;
         return;
    }

    std::string all_trajectories_file = input_base_folder + "/all_proton_trajectories.csv";
    std::string output_plot_file_png = input_base_folder + "/proton_trajectories_plot_root.png";
    std::string output_plot_file_pdf = input_base_folder + "/proton_trajectories_plot_root.pdf";
    std::string output_hist_plot_file_png = input_base_folder + "/proton_final_energy_histogram_root.png";
    std::string output_hist_plot_file_pdf = input_base_folder + "/proton_final_energy_histogram_root.pdf";
    
    std::string output_accel_plot_file_png = input_base_folder + "/proton_acceleration_profile_root.png";
    std::string output_accel_plot_file_pdf = input_base_folder + "/proton_acceleration_profile_root.pdf";
    std::string output_velo_plot_file_png = input_base_folder + "/proton_velocity_profile_root.png";
    std::string output_velo_plot_file_pdf = input_base_folder + "/proton_velocity_profile_root.pdf";

    std::string eps_r_file = input_base_folder + "/permittivity.csv";
    std::string x_coords_file = input_base_folder + "/x_coordinates.csv";
    std::string y_coords_file = input_base_folder + "/y_coordinates.csv";

    int Ny_eps, Nx_eps;
    std::vector<std::vector<double>> eps_r_data = load_2d_csv(eps_r_file, Ny_eps, Nx_eps);
    std::vector<double> x_coords = load_coordinates(x_coords_file); 
    std::vector<double> y_coords = load_coordinates(y_coords_file); 

    if (eps_r_data.empty() || x_coords.empty() || y_coords.empty()) {
        std::cerr << "Could not load necessary data (permittivity, coordinates). Exiting." << std::endl;
        return;
    }
    std::cout << "Loaded x_coords size: " << x_coords.size() << ", y_coords size: " << y_coords.size() << std::endl;
    std::cout << "Loaded eps_r_data dimensions: Ny_eps=" << Ny_eps << ", Nx_eps=" << Nx_eps << std::endl;


    double L_total_sim_um = x_coords.empty() ? 0 : x_coords.back(); 
    double H_total_sim_um = y_coords.empty() ? 0 : y_coords.back(); 
    double grid_spacing_h_um = 0.5; 

    // --- Trajectory Plot ---
    TCanvas *c_trajectories = new TCanvas("c_trajectories", "Proton Trajectories", 1200, 800);
    gStyle->SetOptStat(0); 

    TH2F *h_frame = new TH2F("h_frame", Form("Proton Trajectories (plotted in #mum);X (#mum);Y (#mum)"),
                             100, x_coords.front(), x_coords.back(), // Use actual min/max from coords
                             100, y_coords.front(), y_coords.back());
    h_frame->Draw(); 

    double outline_threshold = -1.0;
    if (!eps_r_data.empty()) {
        double eps_min_val = eps_r_data[0][0];
        double eps_max_val = eps_r_data[0][0];
        for (const auto& row : eps_r_data) {
            for (double val : row) {
                if (val < eps_min_val) eps_min_val = val;
                if (val > eps_max_val) eps_max_val = val;
            }
        }

        if (eps_min_val != eps_max_val && eps_max_val > eps_min_val + 0.5) {
            outline_threshold = (eps_min_val + eps_max_val) / 2.0;
            printf("Dynamically calculated outline threshold for trajectories plot: %.2f (from eps_r min/max: %.2f, %.2f)\\n",
                   outline_threshold, eps_min_val, eps_max_val);
        } else if (eps_min_val == eps_max_val) {
            printf("Warning (trajectories plot): Permittivity data contains only one unique value (%.2f). No outlines will be drawn.\\n", eps_min_val);
        } else {
            printf("Warning (trajectories plot): Permittivity values (%.2f-%.2f) are too close. Could not reliably determine outline threshold.\\n", eps_min_val, eps_max_val);
        }
    } else {
        printf("Warning (trajectories plot): Permittivity data (eps_r_data) is empty. Cannot draw outlines.\\n");
    }

    std::vector<TGraph*> contour_graphs_vec;
    TH2D *h_eps_r_for_contour = nullptr;

    if (outline_threshold != -1.0 && !x_coords.empty() && !y_coords.empty() && 
        Nx_eps == static_cast<int>(x_coords.size()) && Ny_eps == static_cast<int>(y_coords.size()) && Ny_eps == static_cast<int>(eps_r_data.size()) && Nx_eps == static_cast<int>(eps_r_data[0].size())) {
        
        h_eps_r_for_contour = new TH2D("h_eps_r_for_contour", "Permittivity Map for Contour",
                                 Nx_eps, x_coords.data(), Ny_eps, y_coords.data()); // Use variable binning

        for (int i = 0; i < Nx_eps; ++i) { // Loop over x bins
            for (int j = 0; j < Ny_eps; ++j) { // Loop over y bins
                // eps_r_data is [row][col] which is [y_idx][x_idx]
                // TH2D SetBinContent is (binx, biny, content) where binx, biny are 1-indexed
                h_eps_r_for_contour->SetBinContent(i + 1, j + 1, eps_r_data[j][i]);
            }
        }
        h_eps_r_for_contour->SetContour(1, &outline_threshold);
        h_eps_r_for_contour->Draw("CONT Z LIST SAME"); 
        gPad->Update(); 

        TObjArray *contours = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
        if (contours) {
            std::cout << "Found " << contours->GetEntries() << " contour levels." << std::endl;
            for (int i = 0; i < contours->GetEntries(); ++i) {
                TList *list = (TList*)contours->At(i);
                std::cout << "  Level " << i << " has " << list->GetSize() << " contour segments." << std::endl;
                for (int j = 0; j < list->GetSize(); ++j) {
                    TGraph *gr = dynamic_cast<TGraph*>(list->At(j)->Clone()); 
                    if (gr) {
                        gr->SetLineColor(kBlue);
                        gr->SetLineStyle(2); 
                        gr->SetLineWidth(2); // Made slightly thicker
                        gr->Draw("L SAME");
                        contour_graphs_vec.push_back(gr); 
                    }
                }
            }
        } else {
             std::cout << "No contours generated by CONT Z LIST." << std::endl;
        }
    } else {
        std::cout << "Skipping structure outline: threshold not determined, or data mismatch." << std::endl;
        if (Nx_eps != static_cast<int>(x_coords.size())) std::cout << "Nx_eps (" << Nx_eps << ") != x_coords.size() (" << x_coords.size() << ")" << std::endl;
        if (Ny_eps != static_cast<int>(y_coords.size())) std::cout << "Ny_eps (" << Ny_eps << ") != y_coords.size() (" << y_coords.size() << ")" << std::endl;
        if (!eps_r_data.empty()){
            if (Ny_eps != static_cast<int>(eps_r_data.size())) std::cout << "Ny_eps (" << Ny_eps << ") != eps_r_data.size() (" << eps_r_data.size() << ")" << std::endl;
            if (Nx_eps != static_cast<int>(eps_r_data[0].size())) std::cout << "Nx_eps (" << Nx_eps << ") != eps_r_data[0].size() (" << eps_r_data[0].size() << ")" << std::endl;
        }
    }


    std::vector<ProtonTrajectory> all_trajectories = load_all_trajectories(all_trajectories_file);
    std::cout << "Loaded " << all_trajectories.size() << " proton trajectories from file." << std::endl;

    std::vector<double> final_energies_eV;
    int successful_protons_count = 0;
    
    int num_protons_in_file = all_trajectories.size();
    int num_trajectories_to_plot = std::min(num_protons_in_file, 1000);
    double convert_to_um = 1e6;

    TGraph* traj_graph_for_legend = nullptr;
    std::vector<TGraph*> drawn_traj_graphs; // To manage memory

    for (int i = 0; i < num_trajectories_to_plot; ++i) {
        const auto& traj = all_trajectories[i];
        if (traj.x_m.empty()) continue;

        std::vector<double> x_vals_um, y_vals_um;
        for (size_t k = 0; k < traj.x_m.size(); ++k) {
            x_vals_um.push_back(traj.x_m[k] * convert_to_um);
            y_vals_um.push_back(traj.y_m[k] * convert_to_um);
        }

        if (!x_vals_um.empty()) {
            TGraph *gTraj = new TGraph(x_vals_um.size(), x_vals_um.data(), y_vals_um.data());
            gTraj->SetLineColorAlpha(kRed, 0.7f); // Added alpha
            gTraj->SetLineWidth(1);
            gTraj->SetLineStyle(1);
            gTraj->Draw("L SAME"); 
            if (i == 0) traj_graph_for_legend = gTraj; 
            drawn_traj_graphs.push_back(gTraj);
        }
    }
    
    h_frame->SetTitle(Form("Proton Trajectories (%d of %d shown, plotted in #mum);X (#mum);Y (#mum)", num_trajectories_to_plot, num_protons_in_file));

    TLegend *legend = new TLegend(0.65, 0.75, 0.88, 0.88); 
    if (!contour_graphs_vec.empty() && contour_graphs_vec[0]) {
         legend->AddEntry(contour_graphs_vec[0], "Structure Outline (from #epsilon_{r})", "l");
    } else {
        TGraph* dummy_outline = new TGraph(); // Create a dummy for the legend entry
        dummy_outline->SetLineColor(kBlue);
        dummy_outline->SetLineStyle(2);
        dummy_outline->SetLineWidth(2);
        legend->AddEntry(dummy_outline, "Structure Outline (not drawn)", "l");
        // delete dummy_outline; // Add to a management list if needed, or handle here
    }
    if (traj_graph_for_legend) {
        legend->AddEntry(traj_graph_for_legend, "Proton Trajectories", "l");
    }
    legend->Draw();

    c_trajectories->SetGrid();
    c_trajectories->Update();
    c_trajectories->SaveAs(output_plot_file_png.c_str());
    c_trajectories->SaveAs(output_plot_file_pdf.c_str());
    std::cout << "Trajectory plot saved to " << output_plot_file_png << " and " << output_plot_file_pdf << std::endl;


    // --- Energy Histogram ---
    double x_success_threshold_m = (L_total_sim_um * 1e-6) - (grid_spacing_h_um * 1e-6 / 2.0);

    for (const auto& traj : all_trajectories) {
        if (traj.x_m.empty()) continue;

        double last_x_m = traj.x_m.back();
        if (last_x_m >= x_success_threshold_m) {
            successful_protons_count++;
            double vx_mps = traj.vx_m_per_s.back();
            double vy_mps = traj.vy_m_per_s.back();
            double v_sq_mps = vx_mps * vx_mps + vy_mps * vy_mps;
            double ke_joules = 0.5 * M_PROTON_KG * v_sq_mps;
            double ke_eV = ke_joules / E_CHARGE_C;
            final_energies_eV.push_back(ke_eV);
        }
    }

    TH1F *h_energy = nullptr; // Declare outside if
    TCanvas *c_hist_energy = nullptr;

    if (!final_energies_eV.empty()) {
        std::cout << "Number of protons considered successful for energy histogram: " << successful_protons_count << std::endl;
        c_hist_energy = new TCanvas("c_hist_energy", "Final Kinetic Energy Histogram", 1000, 600);
        
        double min_ke_val = final_energies_eV[0], max_ke_val = final_energies_eV[0];
        double sum_ke = 0;
        for (double ke : final_energies_eV) {
            if (ke < min_ke_val) min_ke_val = ke;
            if (ke > max_ke_val) max_ke_val = ke;
            sum_ke += ke;
        }
        double mean_ke_calc = sum_ke / final_energies_eV.size();
        double std_dev_ke_calc = 0;
        for (double ke : final_energies_eV) {
            std_dev_ke_calc += (ke - mean_ke_calc) * (ke - mean_ke_calc);
        }
        std_dev_ke_calc = final_energies_eV.size() > 1 ? std::sqrt(std_dev_ke_calc / (final_energies_eV.size() -1)) : 0;


        h_energy = new TH1F("h_energy", "Histogram of Final Kinetic Energies;Final Kinetic Energy (eV);Number of Protons",
                                  50, min_ke_val - (max_ke_val-min_ke_val)*0.05, max_ke_val + (max_ke_val-min_ke_val)*0.1); 
        for (double ke : final_energies_eV) {
            h_energy->Fill(ke);
        }
        h_energy->SetFillColor(TColor::GetColor(152, 251, 152)); // mediumseagreen
        h_energy->SetLineColor(kBlack);
        h_energy->SetFillStyle(1001); 
        h_energy->Draw();

        c_hist_energy->SetGrid();
        
        double efficiency = (num_protons_in_file > 0) ? (double)successful_protons_count / num_protons_in_file : 0;
        TLatex *stats_text = new TLatex();
        stats_text->SetNDC();
        stats_text->SetTextSize(0.025); 
        stats_text->SetTextAlign(33); 
        
        // Create a pseudo-box with TPad or draw lines for bbox
        TPad *text_pad = new TPad("text_pad", "text_pad", 1, 1, 2, 2); // x1,y1,x2,y2 NDC
        text_pad->SetFillStyle(4000); // Transparent
        text_pad->SetFillColorAlpha(TColor::GetColor("#FFDEAD"), 0.5); // wheat with alpha
        text_pad->SetBorderSize(1); // Border for the box
        text_pad->Draw();
        text_pad->cd(); // Draw text inside this pad

        stats_text->DrawLatex(0.95, 0.90, Form("Mean: %.2f eV", mean_ke_calc));
        stats_text->DrawLatex(0.95, 0.80, Form("Std Dev: %.2f eV", std_dev_ke_calc));
        stats_text->DrawLatex(0.95, 0.70, Form("Min: %.2f eV", min_ke_val));
        stats_text->DrawLatex(0.95, 0.60, Form("Max: %.2f eV", max_ke_val));
        stats_text->DrawLatex(0.95, 0.50, Form("Count: %zu", final_energies_eV.size()));
        stats_text->DrawLatex(0.95, 0.40, Form("Efficiency: %.2f%%", efficiency * 100.0));
        
        c_hist_energy->cd(); // Switch back to main canvas before saving
        c_hist_energy->Update();
        c_hist_energy->SaveAs(output_hist_plot_file_png.c_str());
        c_hist_energy->SaveAs(output_hist_plot_file_pdf.c_str());
        std::cout << "Energy histogram saved to " << output_hist_plot_file_png << " and " << output_hist_plot_file_pdf << std::endl;

    } else {
        std::cout << "No successful protons found to generate an energy histogram." << std::endl;
    }    // --- Acceleration and Velocity Profiles vs Position (Enhanced) ---
    std::vector<double> all_acceleration_x_pos, all_acceleration_mag;
    std::vector<double> all_velocity_x_pos, all_velocity_mag;

    // Define spatial binning parameters (matching Python approach)
    double x_min_m = x_coords.empty() ? 0 : x_coords.front() * 1e-6; // Convert from um to m
    double x_max_m = x_coords.empty() ? 1e-6 : x_coords.back() * 1e-6;
    const int num_bins = 100; // Number of spatial bins
    double x_bin_size_m = (x_max_m - x_min_m) / num_bins;
    
    for (const auto& traj : all_trajectories) {
        if (traj.time_s.empty() || traj.x_m.empty()) continue;

        for (size_t k = 0; k < traj.time_s.size(); ++k) {
            double x_curr = traj.x_m[k];
            double vx_curr = traj.vx_m_per_s[k];
            double vy_curr = traj.vy_m_per_s[k];
            double v_mag_curr = std::sqrt(vx_curr * vx_curr + vy_curr * vy_curr);

            // Bin the x position
            double x_binned = std::round(x_curr / x_bin_size_m) * x_bin_size_m;

            x_pos_to_vx_list[x_binned].push_back(vx_curr);
            x_pos_to_vy_list[x_binned].push_back(vy_curr);
            x_pos_to_v_mag_list[x_binned].push_back(v_mag_curr);

            if (k > 0) { // Acceleration requires a previous point
                double t_prev = traj.time_s[k-1];
                double t_curr = traj.time_s[k];
                double dt = t_curr - t_prev;

                if (dt > 1e-12) { // Avoid division by zero or very small dt
                    double vx_prev = traj.vx_m_per_s[k-1];
                    double vy_prev = traj.vy_m_per_s[k-1];
                    
                    double ax = (vx_curr - vx_prev) / dt;
                    double ay = (vy_curr - vy_prev) / dt;
                    double a_mag = std::sqrt(ax * ax + ay * ay);

                    x_pos_to_ax_list[x_binned].push_back(ax);
                    x_pos_to_ay_list[x_binned].push_back(ay);
                    x_pos_to_a_mag_list[x_binned].push_back(a_mag);
                }
            }
        }
    }

    // Prepare data for TGraphs (Acceleration vs X position)
    std::vector<double> common_x_points_accel;
    std::vector<double> avg_ax_vec, avg_ay_vec, avg_a_mag_vec;

    for (auto const& [x_pos, ax_values] : x_pos_to_ax_list) {
        if (ax_values.empty()) continue;

        auto it_ay = x_pos_to_ay_list.find(x_pos);
        auto it_a_mag = x_pos_to_a_mag_list.find(x_pos);

        if (it_ay != x_pos_to_ay_list.end() && !it_ay->second.empty() &&
            it_a_mag != x_pos_to_a_mag_list.end() && !it_a_mag->second.empty()) {
            
            common_x_points_accel.push_back(x_pos * 1e6); // Convert back to micrometers for plotting
            avg_ax_vec.push_back(std::accumulate(ax_values.begin(), ax_values.end(), 0.0) / ax_values.size());
            avg_ay_vec.push_back(std::accumulate(it_ay->second.begin(), it_ay->second.end(), 0.0) / it_ay->second.size());
            avg_a_mag_vec.push_back(std::accumulate(it_a_mag->second.begin(), it_a_mag->second.end(), 0.0) / it_a_mag->second.size());
        }
    }

    // Prepare data for TGraphs (Velocity vs X position)
    std::vector<double> common_x_points_velo;
    std::vector<double> avg_vx_vec, avg_vy_vec, avg_v_mag_vec;

    for (auto const& [x_pos, vx_values] : x_pos_to_vx_list) {
        if (vx_values.empty()) continue;

        auto it_vy = x_pos_to_vy_list.find(x_pos);
        auto it_v_mag = x_pos_to_v_mag_list.find(x_pos);

        if (it_vy != x_pos_to_vy_list.end() && !it_vy->second.empty() &&
            it_v_mag != x_pos_to_v_mag_list.end() && !it_v_mag->second.empty()) {

            common_x_points_velo.push_back(x_pos * 1e6); // Convert back to micrometers for plotting
            avg_vx_vec.push_back(std::accumulate(vx_values.begin(), vx_values.end(), 0.0) / vx_values.size());
            avg_vy_vec.push_back(std::accumulate(it_vy->second.begin(), it_vy->second.end(), 0.0) / it_vy->second.size());
            avg_v_mag_vec.push_back(std::accumulate(it_v_mag->second.begin(), it_v_mag->second.end(), 0.0) / it_v_mag->second.size());
        }
    }    TCanvas *c_accel = nullptr;
    TMultiGraph *mg_accel = nullptr;
    TLegend *leg_accel = nullptr;

    if (!common_x_points_accel.empty()) {
        c_accel = new TCanvas("c_accel", "Average Acceleration vs Position", 1200, 800);
        c_accel->SetGrid();
        mg_accel = new TMultiGraph("mg_accel", "Average Proton Acceleration vs Position;X Position (#mum);Acceleration (m/s^{2})");

        TGraph *g_avg_ax = new TGraph(common_x_points_accel.size(), common_x_points_accel.data(), avg_ax_vec.data());
        g_avg_ax->SetLineColor(kBlue);
        g_avg_ax->SetLineWidth(2);
        mg_accel->Add(g_avg_ax, "L");

        TGraph *g_avg_ay = new TGraph(common_x_points_accel.size(), common_x_points_accel.data(), avg_ay_vec.data());
        g_avg_ay->SetLineColor(kGreen+2);
        g_avg_ay->SetLineWidth(2);
        mg_accel->Add(g_avg_ay, "L");

        TGraph *g_avg_a_mag = new TGraph(common_x_points_accel.size(), common_x_points_accel.data(), avg_a_mag_vec.data());
        g_avg_a_mag->SetLineColor(kRed);
        g_avg_a_mag->SetLineWidth(2);
        mg_accel->Add(g_avg_a_mag, "L");
        
        mg_accel->Draw("A"); // Draw axis for multigraph

        leg_accel = new TLegend(0.7, 0.75, 0.88, 0.88);
        leg_accel->AddEntry(g_avg_ax, "Avg. A_{x}", "l");
        leg_accel->AddEntry(g_avg_ay, "Avg. A_{y}", "l");
        leg_accel->AddEntry(g_avg_a_mag, "Avg. |A|", "l");
        leg_accel->Draw();

        c_accel->Update();
        c_accel->SaveAs(output_accel_plot_file_png.c_str());
        c_accel->SaveAs(output_accel_plot_file_pdf.c_str());
        std::cout << "Acceleration vs position plot saved to " << output_accel_plot_file_png << " and " << output_accel_plot_file_pdf << std::endl;
    } else {
        std::cout << "Not enough data to generate acceleration vs position profiles." << std::endl;
    }

    TCanvas *c_velo = nullptr;
    TMultiGraph *mg_velo = nullptr;
    TLegend *leg_velo = nullptr;

    if (!common_x_points_velo.empty()) {
        c_velo = new TCanvas("c_velo", "Average Velocity vs Position", 1200, 800);
        c_velo->SetGrid();
        mg_velo = new TMultiGraph("mg_velo", "Average Proton Velocity vs Position;X Position (#mum);Velocity (m/s)");

        TGraph *g_avg_vx = new TGraph(common_x_points_velo.size(), common_x_points_velo.data(), avg_vx_vec.data());
        g_avg_vx->SetLineColor(kBlue);
        g_avg_vx->SetLineWidth(2);
        mg_velo->Add(g_avg_vx, "L");

        TGraph *g_avg_vy = new TGraph(common_x_points_velo.size(), common_x_points_velo.data(), avg_vy_vec.data());
        g_avg_vy->SetLineColor(kGreen+2);
        g_avg_vy->SetLineWidth(2);
        mg_velo->Add(g_avg_vy, "L");

        TGraph *g_avg_v_mag = new TGraph(common_x_points_velo.size(), common_x_points_velo.data(), avg_v_mag_vec.data());
        g_avg_v_mag->SetLineColor(kRed);
        g_avg_v_mag->SetLineWidth(2);
        mg_velo->Add(g_avg_v_mag, "L");

        mg_velo->Draw("A");

        leg_velo = new TLegend(0.7, 0.75, 0.88, 0.88);
        leg_velo->AddEntry(g_avg_vx, "Avg. V_{x}", "l");
        leg_velo->AddEntry(g_avg_vy, "Avg. V_{y}", "l");
        leg_velo->AddEntry(g_avg_v_mag, "Avg. |V|", "l");
        leg_velo->Draw();

        c_velo->Update();
        c_velo->SaveAs(output_velo_plot_file_png.c_str());
        c_velo->SaveAs(output_velo_plot_file_pdf.c_str());
        std::cout << "Velocity vs position plot saved to " << output_velo_plot_file_png << " and " << output_velo_plot_file_pdf << std::endl;
    } else {
        std::cout << "Not enough data to generate velocity vs position profiles." << std::endl;
    }

    // Cleanup
    // delete c_trajectories;
    // delete h_frame;
    // if (h_eps_r_for_contour) delete h_eps_r_for_contour;
    // for (TGraph* gr : contour_graphs_vec) delete gr;
    // for (TGraph* gr : drawn_traj_graphs) delete gr; // traj_graph_for_legend is one of these
    // // legend is owned by the canvas, no need to delete explicitly if canvas is deleted.
    
    // if (c_hist_energy) delete c_hist_energy; // h_energy is owned by c_hist_energy
    // // if (h_energy) delete h_energy; // Only if not drawn on a canvas that is deleted

    // if (c_accel) delete c_accel; // Deletes legend leg_accel too
    // if (mg_accel) delete mg_accel; // Deletes graphs g_avg_ax, g_avg_ay, g_avg_a_mag

    // if (c_velo) delete c_velo; // Deletes legend leg_velo too
    // if (mg_velo) delete mg_velo; // Deletes graphs g_avg_vx, g_avg_vy, g_avg_v_mag

    // std::cout << "ROOT script processing finished." << std::endl;
    // // gROOT->SetBatch(kFALSE); // Turn off batch mode if it was set

    // If you want to keep the application alive to see plots interactively when not in batch mode:
    // TApplication theApp("App", nullptr, nullptr);
    // theApp.Run();
}

// Main guard for script execution (e.g. when compiling with ACLiC or standalone)
#ifndef __CINT__ // or __CLING__ for newer ROOT versions
int main(int argc, char **argv) {
    std::string folder_path = "";
    if (argc > 1) {
        folder_path = argv[1];
    }
    plot_trajectories_main(folder_path);
    return 0;
}
#endif
