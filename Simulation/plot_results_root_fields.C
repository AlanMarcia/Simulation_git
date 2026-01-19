// ROOT macro to plot electric field results from simulation
// Usage: root -l 'plot_results_root_fields.C("folder_name")'
// Or in ROOT: .x plot_results_root_fields.C("folder_name")

#include <TCanvas.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TStyle.h>
#include <TColor.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TArrow.h>
#include <TMath.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>

using namespace std;

// Structure to hold 2D field data
struct FieldData {
    int nx, ny;
    double xmin, xmax, ymin, ymax;
    double* x_coords;
    double* y_coords;
    double** potential;
    double** Ex;
    double** Ey;
    double** eps_r;
    double** E_mag;
    
    FieldData() : nx(0), ny(0), xmin(0), xmax(0), ymin(0), ymax(0),
                  x_coords(nullptr), y_coords(nullptr), potential(nullptr),
                  Ex(nullptr), Ey(nullptr), eps_r(nullptr), E_mag(nullptr) {}
};

// Load 1D CSV data - ROOT-compatible version
double* load_1d_csv(const string& filename, int& size) {
    // First pass: count elements
    ifstream file1(filename);
    if (!file1.is_open()) {
        cerr << "Error: Could not open file " << filename << endl;
        size = 0;
        return nullptr;
    }
    
    size = 0;
    string line;
    while (getline(file1, line)) {
        stringstream ss(line);
        string value;
        while (getline(ss, value, ',')) size++;
    }
    file1.close();
    
    if (size == 0) return nullptr;
    
    // Allocate array
    double* data = new double[size];
    
    // Second pass: load data
    ifstream file2(filename);
    int idx = 0;
    while (getline(file2, line)) {
        stringstream ss(line);
        string value;
        while (getline(ss, value, ',') && idx < size) {
            data[idx++] = stod(value);
        }
    }
    file2.close();
    
    return data;
}

// Load 2D CSV data - ROOT-compatible version
double** load_2d_csv(const string& filename, int& rows, int& cols) {
    // First pass: count dimensions
    ifstream file1(filename);
    if (!file1.is_open()) {
        cerr << "Error: Could not open file " << filename << endl;
        rows = 0;
        cols = 0;
        return nullptr;
    }
    
    string line;
    rows = 0;
    cols = 0;
    while (getline(file1, line)) {
        if (cols == 0) {
            stringstream ss(line);
            string value;
            while (getline(ss, value, ',')) cols++;
        }
        rows++;
    }
    file1.close();
    
    if (rows == 0 || cols == 0) return nullptr;
    
    // Allocate 2D array
    double** data = new double*[rows];
    for (int i = 0; i < rows; i++) {
        data[i] = new double[cols];
    }
    
    // Second pass: load data
    ifstream file2(filename);
    int row = 0;
    while (getline(file2, line) && row < rows) {
        stringstream ss(line);
        string value;
        int col = 0;
        while (getline(ss, value, ',') && col < cols) {
            data[row][col] = stod(value);
            col++;
        }
        row++;
    }
    file2.close();
    
    return data;
}

// Load geometry parameters
map<string, double> load_geometry_params(const string& filename) {
    map<string, double> params;
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Warning: Could not open geometry parameters file " << filename << endl;
        return params;
    }
    
    string line;
    while (getline(file, line)) {
        stringstream ss(line);
        string key, value;
        if (getline(ss, key, ',') && getline(ss, value, ',')) {
            try {
                params[key] = stod(value);
            } catch (...) {
                // Skip non-numeric values
            }
        }
    }
    file.close();
    return params;
}

// Create contour plot with outlines
TH2D* create_contour_plot(const string& name, const string& title,
                          const FieldData& data,
                          double** field,
                          bool mask_dielectric = false,
                          double outline_threshold = -1) {
    
    TH2D* hist = new TH2D(name.c_str(), title.c_str(),
                          data.nx, data.xmin, data.xmax,
                          data.ny, data.ymin, data.ymax);
    
    // Fill histogram
    for (int iy = 0; iy < data.ny; iy++) {
        for (int ix = 0; ix < data.nx; ix++) {
            double value = field[iy][ix];
            
            // Mask dielectric regions if requested
            if (mask_dielectric && outline_threshold > 0) {
                if (data.eps_r[iy][ix] >= outline_threshold) {
                    value = 0; // Set to zero in dielectric
                }
            }
            
            hist->SetBinContent(ix + 1, iy + 1, value);
        }
    }
    
    return hist;
}

// Draw geometry outlines on current pad
void draw_outlines(const FieldData& data, double threshold, int color, int style) {
    if (threshold <= 0) return;
    
    // Simple outline drawing - find boundaries
    for (int iy = 1; iy < data.ny - 1; iy++) {
        for (int ix = 1; ix < data.nx - 1; ix++) {
            bool is_boundary = false;
            double current_eps = data.eps_r[iy][ix];
            
            // Check if this is a boundary point
            if ((current_eps < threshold && data.eps_r[iy][ix+1] >= threshold) ||
                (current_eps >= threshold && data.eps_r[iy][ix+1] < threshold) ||
                (current_eps < threshold && data.eps_r[iy+1][ix] >= threshold) ||
                (current_eps >= threshold && data.eps_r[iy+1][ix] < threshold)) {
                is_boundary = true;
            }
            
            if (is_boundary) {
                // Draw a small marker at boundary
                TMarker* m = new TMarker(data.x_coords[ix], data.y_coords[iy], 1);
                m->SetMarkerColor(color);
                m->SetMarkerSize(0.1);
                m->Draw();
            }
        }
    }
}

void plot_results_root_fields(const string& folder_name = "geometria_Denti_sfasati_profondi_5um") {
    
    cout << "Loading data from folder: " << folder_name << endl;
    
    // Set ROOT style
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kViridis);
    gStyle->SetNumberContours(50);
    
    // Load data
    FieldData data;
    
    string base_path = folder_name + "/";
    int nx_size, ny_size;
    data.x_coords = load_1d_csv(base_path + "x_coordinates.csv", nx_size);
    data.y_coords = load_1d_csv(base_path + "y_coordinates.csv", ny_size);
    
    int nx_temp, ny_temp;
    data.potential = load_2d_csv(base_path + "potential.csv", ny_temp, nx_temp);
    data.Ex = load_2d_csv(base_path + "electric_field_x.csv", ny_temp, nx_temp);
    data.Ey = load_2d_csv(base_path + "electric_field_y.csv", ny_temp, nx_temp);
    data.eps_r = load_2d_csv(base_path + "permittivity.csv", ny_temp, nx_temp);
    
    if (!data.x_coords || !data.y_coords || !data.potential) {
        cerr << "Error: Could not load data files!" << endl;
        return;
    }
    
    data.nx = nx_size;
    data.ny = ny_size;
    data.xmin = data.x_coords[0];
    data.xmax = data.x_coords[data.nx - 1];
    data.ymin = data.y_coords[0];
    data.ymax = data.y_coords[data.ny - 1];
    
    cout << "Data loaded: " << data.nx << " x " << data.ny << " grid" << endl;
    
    // Calculate E magnitude
    data.E_mag = new double*[data.ny];
    for (int iy = 0; iy < data.ny; iy++) {
        data.E_mag[iy] = new double[data.nx];
        for (int ix = 0; ix < data.nx; ix++) {
            data.E_mag[iy][ix] = sqrt(data.Ex[iy][ix] * data.Ex[iy][ix] + 
                                      data.Ey[iy][ix] * data.Ey[iy][ix]);
        }
    }
    
    // Determine thresholds for material boundaries
    double outline_threshold_silicon = -1;
    double outline_threshold_aluminum = -1;
    
    set<double> unique_eps;
    for (int iy = 0; iy < data.ny; iy++) {
        for (int ix = 0; ix < data.nx; ix++) {
            unique_eps.insert(data.eps_r[iy][ix]);
        }
    }
    
    vector<double> eps_values(unique_eps.begin(), unique_eps.end());
    sort(eps_values.begin(), eps_values.end());
    
    if (eps_values.size() >= 3) {
        double val_vacuum = eps_values[0];
        double val_silicon = eps_values[1];
        double val_aluminum = eps_values[eps_values.size() - 1];
        outline_threshold_silicon = (val_vacuum + val_silicon) / 2.0;
        outline_threshold_aluminum = (val_silicon + val_aluminum) / 2.0;
        cout << "Detected 3 materials: vacuum=" << val_vacuum 
             << ", silicon=" << val_silicon 
             << ", aluminum=" << val_aluminum << endl;
    } else if (eps_values.size() >= 2) {
        outline_threshold_silicon = (eps_values[0] + eps_values.back()) / 2.0;
        cout << "Detected 2 materials, threshold=" << outline_threshold_silicon << endl;
    }
    
    // Load geometry parameters
    auto geo_params = load_geometry_params(base_path + "geometry_params.csv");
    
    // Calculate center of vacuum gap for profile plot
    int y_center_idx = data.ny / 2; // Default to middle
    if (geo_params.count("y_si_base_height") && geo_params.count("y_teeth_height") && 
        geo_params.count("y_vacuum_gap_thick")) {
        double y_center = geo_params["y_vacuum_padding_bottom"] + 
                         geo_params["y_si_base_height"] + 
                         geo_params["y_teeth_height"] + 
                         geo_params["y_vacuum_gap_thick"] / 2.0;
        // Find closest index
        double min_dist = 1e10;
        for (int i = 0; i < data.ny; i++) {
            double dist = fabs(data.y_coords[i] - y_center);
            if (dist < min_dist) {
                min_dist = dist;
                y_center_idx = i;
            }
        }
        cout << "Profile plot at y = " << data.y_coords[y_center_idx] << " um (index " << y_center_idx << ")" << endl;
    }
    
    // ============= PLOT 1: Electric Potential =============
    TCanvas* c1 = new TCanvas("c1", "Electric Potential", 1000, 800);
    c1->SetRightMargin(0.15);
    
    TH2D* h_potential = create_contour_plot("h_potential", "Electric Potential;x (#mum);y (#mum)",
                                           data, data.potential);
    h_potential->Draw("COLZ");
    h_potential->GetZaxis()->SetTitle("Potential (V)");
    
    c1->SaveAs((folder_name + "/potential_plot_root.png").c_str());
    c1->SaveAs((folder_name + "/potential_plot_root.pdf").c_str());
    
    // ============= PLOT 2: Electric Field Magnitude (Vacuum Only) =============
    TCanvas* c2 = new TCanvas("c2", "Electric Field Magnitude", 1000, 800);
    c2->SetRightMargin(0.15);
    
    // Use hot colormap for E-field
    gStyle->SetPalette(kTemperatureMap);
    
    TH2D* h_emag = create_contour_plot("h_emag", "Electric Field Magnitude |E| in Vacuum;x (#mum);y (#mum)",
                                      data, data.E_mag, true, outline_threshold_silicon);
    h_emag->Draw("COLZ");
    h_emag->GetZaxis()->SetTitle("Electric Field Magnitude (V/#mum)");
    h_emag->GetZaxis()->SetTitleOffset(1.3);
    
    c2->SaveAs((folder_name + "/efield_magnitude_plot_root.png").c_str());
    c2->SaveAs((folder_name + "/efield_magnitude_plot_root.pdf").c_str());
    
    // ============= PLOT 3: Permittivity Map =============
    gStyle->SetPalette(kCool);
    TCanvas* c3 = new TCanvas("c3", "Permittivity Map", 1000, 800);
    c3->SetRightMargin(0.15);
    
    TH2D* h_eps = create_contour_plot("h_eps", "Relative Permittivity Map;x (#mum);y (#mum)",
                                     data, data.eps_r);
    h_eps->Draw("COLZ");
    h_eps->GetZaxis()->SetTitle("Relative Permittivity #varepsilon_{r}");
    
    c3->SaveAs((folder_name + "/permittivity_map_plot_root.png").c_str());
    c3->SaveAs((folder_name + "/permittivity_map_plot_root.pdf").c_str());
    
    // ============= PLOT 4: Electric Field Quiver (Vacuum Only) =============
    TCanvas* c4 = new TCanvas("c4", "Electric Field Vectors", 1200, 900);
    c4->SetRightMargin(0.15);
    
    // Create base histogram for axes
    TH2D* h_base = new TH2D("h_base", "Electric Field Vectors in Vacuum;x (#mum);y (#mum)",
                           100, data.xmin, data.xmax, 100, data.ymin, data.ymax);
    h_base->Draw();
    
    // Draw arrows (downsampled)
    int skip = 5;
    double scale_factor = (data.xmax - data.xmin) / data.nx * 2.0; // Adjust arrow scaling
    
    for (int iy = 0; iy < data.ny; iy += skip) {
        for (int ix = 0; ix < data.nx; ix += skip) {
            // Only draw in vacuum
            if (data.eps_r[iy][ix] < outline_threshold_silicon) {
                double x = data.x_coords[ix];
                double y = data.y_coords[iy];
                double ex = data.Ex[iy][ix];
                double ey = data.Ey[iy][ix];
                double mag = data.E_mag[iy][ix];
                
                if (mag > 1e-10) {
                    // Normalize and scale
                    double dx = ex / mag * scale_factor;
                    double dy = ey / mag * scale_factor;
                    
                    TArrow* arrow = new TArrow(x, y, x + dx, y + dy, 0.01, "|>");
                    arrow->SetLineColor(kBlue);
                    arrow->SetFillColor(kBlue);
                    arrow->SetLineWidth(1);
                    arrow->Draw();
                }
            }
        }
    }
    
    c4->SaveAs((folder_name + "/efield_quiver_vacuum_plot_root.png").c_str());
    c4->SaveAs((folder_name + "/efield_quiver_vacuum_plot_root.pdf").c_str());
    
    // ============= PLOT 5: Profile at center of vacuum gap =============
    TCanvas* c5 = new TCanvas("c5", "Field Profile", 1400, 600);
    
    // Extract profiles
    double* V_profile = new double[data.nx];
    double* Ex_profile = new double[data.nx];
    double* Ey_profile = new double[data.nx];
    double* Emag_profile = new double[data.nx];
    
    for (int ix = 0; ix < data.nx; ix++) {
        V_profile[ix] = data.potential[y_center_idx][ix];
        Ex_profile[ix] = data.Ex[y_center_idx][ix];
        Ey_profile[ix] = data.Ey[y_center_idx][ix];
        Emag_profile[ix] = data.E_mag[y_center_idx][ix];
    }
    
    TGraph* g_V = new TGraph(data.nx, data.x_coords, V_profile);
    TGraph* g_Emag = new TGraph(data.nx, data.x_coords, Emag_profile);
    
    g_V->SetLineColor(kBlue);
    g_V->SetLineWidth(2);
    g_V->SetTitle(Form("Profile at y = %.2f #mum;x (#mum);Potential (V)", data.y_coords[y_center_idx]));
    
    g_Emag->SetLineColor(kRed);
    g_Emag->SetLineWidth(2);
    g_Emag->SetLineStyle(2);
    
    g_V->Draw("AL");
    
    // Add second y-axis for E-field (overlay)
    TGraph* g_Emag_overlay = (TGraph*)g_Emag->Clone();
    // Scale E-field to fit on same plot
    double V_max = TMath::MaxElement(data.nx, V_profile);
    double E_max = TMath::MaxElement(data.nx, Emag_profile);
    if (E_max > 0) {
        for (int i = 0; i < data.nx; i++) {
            g_Emag_overlay->SetPoint(i, data.x_coords[i], Emag_profile[i] * V_max / E_max * 0.8);
        }
    }
    g_Emag_overlay->Draw("L SAME");
    
    TLegend* leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg->AddEntry(g_V, "Potential (V)", "l");
    leg->AddEntry(g_Emag_overlay, "|E| (V/#mum) scaled", "l");
    leg->Draw();
    
    c5->SaveAs((folder_name + "/center_gap_profile_plot_root.png").c_str());
    c5->SaveAs((folder_name + "/center_gap_profile_plot_root.pdf").c_str());
    
    // ============= PLOT 6: E-field components profile =============
    TCanvas* c6 = new TCanvas("c6", "E-field Components", 1400, 1200);
    c6->Divide(1, 3);
    
    TGraph* g_Ex = new TGraph(data.nx, data.x_coords, Ex_profile);
    TGraph* g_Ey = new TGraph(data.nx, data.x_coords, Ey_profile);
    TGraph* g_Emag_comp = new TGraph(data.nx, data.x_coords, Emag_profile);
    
    c6->cd(1);
    g_Ex->SetLineColor(kBlue+2);
    g_Ex->SetLineWidth(2);
    g_Ex->SetTitle("Ex Component;x (#mum);Ex (V/#mum)");
    g_Ex->Draw("AL");
    gPad->SetGrid();
    
    c6->cd(2);
    g_Ey->SetLineColor(kGreen+2);
    g_Ey->SetLineWidth(2);
    g_Ey->SetTitle("Ey Component;x (#mum);Ey (V/#mum)");
    g_Ey->Draw("AL");
    gPad->SetGrid();
    
    c6->cd(3);
    g_Emag_comp->SetLineColor(kRed+2);
    g_Emag_comp->SetLineWidth(2);
    g_Emag_comp->SetTitle("|E| Magnitude;x (#mum);|E| (V/#mum)");
    g_Emag_comp->Draw("AL");
    gPad->SetGrid();
    
    // Add statistics text
    double ex_max = TMath::MaxElement(data.nx, Ex_profile);
    double ey_max = TMath::MaxElement(data.nx, Ey_profile);
    double emag_max = TMath::MaxElement(data.nx, Emag_profile);
    double emag_mean = TMath::Mean(data.nx, Emag_profile);
    
    TLatex* text = new TLatex();
    text->SetNDC();
    text->SetTextSize(0.03);
    text->DrawLatex(0.15, 0.85, Form("Max |Ex|: %.3f V/#mum", ex_max));
    text->DrawLatex(0.15, 0.80, Form("Max |Ey|: %.3f V/#mum", ey_max));
    text->DrawLatex(0.15, 0.75, Form("Max |E|: %.3f V/#mum", emag_max));
    text->DrawLatex(0.15, 0.70, Form("Mean |E|: %.3f V/#mum", emag_mean));
    
    c6->SaveAs((folder_name + "/efield_components_profile_plot_root.png").c_str());
    c6->SaveAs((folder_name + "/efield_components_profile_plot_root.pdf").c_str());
    
    cout << "\n=== All plots saved to " << folder_name << " ===" << endl;
    cout << "Files created:" << endl;
    cout << "  - potential_plot_root.png/pdf" << endl;
    cout << "  - efield_magnitude_plot_root.png/pdf" << endl;
    cout << "  - permittivity_map_plot_root.png/pdf" << endl;
    cout << "  - efield_quiver_vacuum_plot_root.png/pdf" << endl;
    cout << "  - center_gap_profile_plot_root.png/pdf" << endl;
    cout << "  - efield_components_profile_plot_root.png/pdf" << endl;
}
