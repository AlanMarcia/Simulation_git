import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os
import glob
import numpy as np # For NaN comparison if needed, and math.
import sys # Import sys module for command-line arguments

# Physical constants
M_PROTON_KG = 1.67262192e-27  # Mass of proton in kg
E_CHARGE_C = 1.60217663e-19   # Elementary charge in Coulombs

# def load_geometry_params(filepath): # No longer needed for plotting geometry
#     """Loads geometry parameters from a CSV file (values expected in µm)."""
#     params = {}
#     try:
#         df = pd.read_csv(filepath, header=None, index_col=0)
#         params = df[1].to_dict()
#     except FileNotFoundError:
#         print(f"Error: Geometry parameters file not found at {filepath}")
#         return None
#     except Exception as e:
#         print(f"Error reading geometry parameters: {e}")
#         return None
#     # Ensure numeric types
#     for key in params:
#         try:
#             params[key] = float(params[key])
#         except ValueError:
#             print(f"Warning: Could not convert param {key} to float.")
#             pass 
#     return params

def load_coordinates(filepath):
    """Loads 1D coordinates from a CSV file (values expected in µm)."""
    try:
        df = pd.read_csv(filepath, header=None)
        return df[0].tolist()
    except FileNotFoundError:
        print(f"Error: Coordinates file not found at {filepath}")
        return None
    except Exception as e:
        print(f"Error reading coordinates file {filepath}: {e}")
        return None

def load_2d_csv(filename):
    """Loads 2D data from a CSV file, returns (Ny, Nx) array."""
    if not os.path.exists(filename):
        print(f"Error: File {filename} not found.")
        return None
    try:
        data = np.loadtxt(filename, delimiter=',')
        return data # Assumes CSV is saved Ny rows, Nx columns
    except Exception as e:
        print(f"Error loading 2D CSV file {filename}: {e}")
        return None

def main(folder_path=None): # Add folder_path argument
    if folder_path:
        input_base_folder = folder_path
        print(f"Using specified folder for trajectory data: {input_base_folder}")
    else:
        input_base_folder = "geometria_Denti_uguali_5um" # Default base folder for input files
        print(f"No folder specified, using default for trajectory data: {input_base_folder}")

    if not os.path.isdir(input_base_folder):
        print(f"Error: The specified folder '{input_base_folder}' does not exist or is not a directory.")
        print("Please ensure the C++ simulation has run and created this folder with CSV files.")
        return

    # trajectories_folder = os.path.join(input_base_folder, "proton_trajectories") # Removed
    all_trajectories_file = os.path.join(input_base_folder, "all_proton_trajectories.csv") # New single file
    output_plot_file = os.path.join(input_base_folder, "proton_trajectories_plot.png") # Changed to .png
    output_hist_plot_file = os.path.join(input_base_folder, "proton_final_energy_histogram.png") # Changed to .png

    # geom_params_file = os.path.join(input_base_folder, "geometry_params.csv") # File still exists, but not used for plotting shapes
    eps_r_file = os.path.join(input_base_folder, "permittivity.csv")
    x_coords_file = os.path.join(input_base_folder, "x_coordinates.csv")
    y_coords_file = os.path.join(input_base_folder, "y_coordinates.csv")

    # geom = load_geometry_params(geom_params_file) # geom dimensions are in µm - No longer used for plotting shapes
    eps_r_data = load_2d_csv(eps_r_file) # Expected shape (Ny, Nx)
    x_coords = load_coordinates(x_coords_file)   # x_coords are in µm
    y_coords = load_coordinates(y_coords_file)   # y_coords are in µm

    if eps_r_data is None or x_coords is None or y_coords is None:
        print("Could not load necessary data (permittivity, coordinates). Exiting.")
        return

    L_total_sim_um = x_coords[-1] if x_coords else 0 # µm
    H_total_sim_um = y_coords[-1] if y_coords else 0 # µm
    
    # Default grid spacing for tolerance, as geom_params is not loaded for 'h'
    # This value should ideally match the 'h' used in the simulation for accurate success check.
    # If geometry_params.csv is available and 'h' is needed, it could be loaded separately for just this value.
    grid_spacing_h_um = 0.5 # µm, default for tolerance, adjust if needed

    fig, ax = plt.subplots(figsize=(12, 8))

    # Plot geometry outlines using permittivity map
    X_mesh, Y_mesh = np.meshgrid(x_coords, y_coords)
    eps_vac_assumed = 1.0
    eps_si_assumed = 11.7 # Typical value for silicon
    outline_threshold = (eps_vac_assumed + eps_si_assumed) / 2.0
    
    ax.contour(X_mesh, Y_mesh, eps_r_data, levels=[outline_threshold], colors='blue', linewidths=0.8, linestyles='--')
    
    # Add a proxy artist for the legend entry for structure outline
    structure_outline_proxy = plt.Line2D([0], [0], linestyle="--", color="blue", label='Structure Outline (from εr)')
    # ax.legend(handles=[structure_outline_proxy], loc='upper right') # Initial legend setup

    final_energies_eV = []
    successful_protons_count = 0

    # trajectory_files = glob.glob(os.path.join(trajectories_folder, "proton_*_trajectory.csv")) # Removed
    
    # Load the single trajectory file
    try:
        df_all_trajectories = pd.read_csv(all_trajectories_file)
        if df_all_trajectories.empty:
            print(f"Trajectory file {all_trajectories_file} is empty.")
            df_all_trajectories = None # Set to None to skip plotting/analysis
    except FileNotFoundError:
        print(f"Trajectory file not found: {all_trajectories_file}")
        df_all_trajectories = None
    except Exception as e:
        print(f"Error reading trajectory file {all_trajectories_file}: {e}")
        df_all_trajectories = None

    
    plotted_traj_legend_added = False
    num_protons_in_file = 0
    
    if df_all_trajectories is not None:
        grouped_trajectories = df_all_trajectories.groupby('proton_id')
        proton_ids = list(grouped_trajectories.groups.keys())
        num_protons_in_file = len(proton_ids)
        
        num_trajectories_to_plot = min(num_protons_in_file, 1000)

        # Plot a subset of trajectories
        for i, proton_id in enumerate(proton_ids):
            if i >= num_trajectories_to_plot:
                break # Stop plotting if limit is reached
            
            df_traj = grouped_trajectories.get_group(proton_id)
            if not df_traj.empty:
                # Convert m to µm for plotting
                line, = ax.plot(df_traj['x_m'] * 1e6, df_traj['y_m'] * 1e6, linestyle='-', linewidth=0.5, alpha=0.7, color='red')
                if not plotted_traj_legend_added:
                    line.set_label('Proton Trajectories') # Add label to one line for the legend
                    plotted_traj_legend_added = True
        
        # Process all trajectories for energy histogram
        for proton_id in proton_ids:
            df_traj = grouped_trajectories.get_group(proton_id)
            if not df_traj.empty:
                last_point = df_traj.iloc[-1]
                # last_point['x_m'] is in meters. Compare with L_total_sim_um (converted to m)
                if last_point['x_m'] >= (L_total_sim_um * 1e-6) - (grid_spacing_h_um * 1e-6 / 2.0) :
                    successful_protons_count += 1
                    # Velocities are already in m/s from the CSV
                    vx_mps = last_point['vx_m_per_s']
                    vy_mps = last_point['vy_m_per_s']
                                        
                    v_sq_mps = vx_mps**2 + vy_mps**2
                    
                    ke_joules = 0.5 * M_PROTON_KG * v_sq_mps
                    ke_eV = ke_joules / E_CHARGE_C
                    final_energies_eV.append(ke_eV)
    else:
        print(f"No trajectory data loaded from {all_trajectories_file}. Skipping trajectory plotting and energy analysis.")
        num_trajectories_to_plot = 0


    ax.set_xlabel("X (µm)") # Axis labels remain in µm for visualization
    ax.set_ylabel("Y (µm)") # Axis labels remain in µm for visualization
    ax.set_title(f"Proton Trajectories ({num_trajectories_to_plot} of {num_protons_in_file} shown, plotted in µm)")
    
    # Set plot limits based on the simulation domain (in µm)
    ax.set_xlim(0, L_total_sim_um)
    ax.set_ylim(0, H_total_sim_um)
    
    # Update legend to include both structure and trajectories if any were plotted
    handles, labels = ax.get_legend_handles_labels()
    custom_handles = [structure_outline_proxy]
    if any(label == 'Proton Trajectories' for label in labels):
         # Find the trajectory handle to ensure it's included
        traj_handle = next((h for h, l in zip(handles, labels) if l == 'Proton Trajectories'), None)
        if traj_handle:
            custom_handles.append(traj_handle)
    
   # ax.legend(handles=custom_handles, loc='upper right')
    ax.set_aspect('equal', adjustable='box') # Important for correct geometric representation
    plt.grid(True, linestyle=':', alpha=0.7)
    
    try:
        plt.savefig(output_plot_file, dpi=300) # Corrected to savefig, adjusted dpi
        print(f"Plot saved to {output_plot_file}")
    except Exception as e:
        print(f"Error saving plot: {e}")
    
    plt.show() # Show trajectory plot
    plt.close(fig) # Close the trajectory plot figure

    # Plotting the histogram of final energies
    if final_energies_eV:
        print(f"Number of protons considered successful for energy histogram: {successful_protons_count}")
        fig_hist, ax_hist = plt.subplots(figsize=(10, 6))
        ax_hist.hist(final_energies_eV, bins=50, edgecolor='black', alpha=0.75, color='mediumseagreen') # Added color for consistency
        ax_hist.set_xlabel("Final Kinetic Energy (eV)")
        ax_hist.set_ylabel("Number of Protons")
        ax_hist.set_title("Histogram of Final Kinetic Energies of Successful Protons")
        ax_hist.grid(True, linestyle=':', alpha=0.7)

        # Calculate and display statistics
        mean_ke = np.mean(final_energies_eV)
        std_ke = np.std(final_energies_eV)
        min_ke = np.min(final_energies_eV)
        max_ke = np.max(final_energies_eV)
        num_protons = len(final_energies_eV)
        Efficiency = successful_protons_count / num_protons_in_file if num_protons_in_file > 0 else 0
        stats_text_ke = f'Mean: {mean_ke:.2f} eV\nStd Dev: {std_ke:.2f} eV\nMin: {min_ke:.2f} eV\nMax: {max_ke:.2f} eV \nCount: {num_protons} Efficiency: {Efficiency:.2%}'
        
        # Position the text box in the upper right corner of the plot
        ax_hist.text(0.95, 0.95, stats_text_ke, transform=ax_hist.transAxes, fontsize=9,
                     verticalalignment='top', horizontalalignment='right',
                     bbox=dict(boxstyle='round,pad=0.5', fc='wheat', alpha=0.5))
        
        plt.tight_layout() # Adjust layout to prevent text box from being clipped
        try:
            plt.savefig(output_hist_plot_file, dpi=300) 
            print(f"Energy histogram saved to {output_hist_plot_file}")
        except Exception as e:
            print(f"Error saving energy histogram: {e}")
        plt.show() # Show histogram plot
        plt.close(fig_hist)
    else:
        print("No successful protons found to generate an energy histogram.")

if __name__ == "__main__":
    cli_folder_path = None
    if len(sys.argv) > 1:
        cli_folder_path = sys.argv[1]
    main(folder_path=cli_folder_path)
