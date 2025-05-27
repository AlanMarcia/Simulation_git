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
    output_accel_plot_file = os.path.join(input_base_folder, "proton_acceleration_profile.png") # New plot file
    output_vel_plot_file = os.path.join(input_base_folder, "proton_velocity_profile.png") # New plot file for velocity

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
    
    # Dynamically calculate outline_threshold from eps_r_data
    outline_threshold = None
    if eps_r_data is not None:
        unique_eps_values = np.unique(eps_r_data)
        unique_eps_values.sort() # Ensure sorted
        if len(unique_eps_values) >= 2:
            val_low = unique_eps_values[0]
            val_high = unique_eps_values[-1]
            if val_high > val_low + 0.5: # Heuristic: difference must be at least 0.5
                outline_threshold = (val_low + val_high) / 2.0
                print(f"Dynamically calculated outline threshold for trajectories plot: {outline_threshold:.2f} (from eps_r min/max: {val_low:.2f}, {val_high:.2f})")
            else:
                print(f"Warning (trajectories plot): Unique permittivity values ({unique_eps_values}) are too close. Could not reliably determine outline threshold.")
        elif len(unique_eps_values) == 1:
            print(f"Warning (trajectories plot): Permittivity data contains only one unique value ({unique_eps_values[0]:.2f}). No outlines will be drawn.")
        else:
            print("Warning (trajectories plot): Could not process permittivity data for outlines.")
    else:
        print("Warning (trajectories plot): Permittivity data (eps_r_data) is None. Cannot draw outlines.")

    # eps_vac_assumed = 1.0 # Old method
    # eps_si_assumed = 11.7 # Typical value for silicon # Old method
    # outline_threshold = (eps_vac_assumed + eps_si_assumed) / 2.0 # Old method
    
    if outline_threshold is not None:
        ax.contour(X_mesh, Y_mesh, eps_r_data, levels=[outline_threshold], colors='blue', linewidths=0.8, linestyles='--')
        # Add a proxy artist for the legend entry for structure outline
        structure_outline_proxy = plt.Line2D([0], [0], linestyle="--", color="blue", label='Structure Outline (from εr)')
    else:
        # Create a dummy proxy if no outline is drawn, so legend logic doesn't break
        structure_outline_proxy = plt.Line2D([0], [0], linestyle="none", label='Structure Outline (not drawn)')
        print("Skipping structure outline on trajectory plot as threshold could not be determined.")
    
    # ax.legend(handles=[structure_outline_proxy], loc='upper right') # Initial legend setup

    final_energies_eV = []
    successful_protons_count = 0
    all_acceleration_data = [] # To store {'x_m': x_position, 'a_mag_mps2': acceleration_magnitude}
    all_velocity_data = [] # To store {'x_m': x_position, 'v_mag_mps': velocity_magnitude}

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
        
        # Process all trajectories for energy histogram and acceleration
        for proton_id in proton_ids:
            df_traj = grouped_trajectories.get_group(proton_id)
            if not df_traj.empty:
                # Energy calculation (existing code)
                last_point = df_traj.iloc[-1]
                if last_point['x_m'] >= (L_total_sim_um * 1e-6) - (grid_spacing_h_um * 1e-6 / 2.0) :
                    successful_protons_count += 1
                    # Velocities are already in m/s from the CSV
                    vx_mps = last_point['vx_m_per_s']
                    vy_mps = last_point['vy_m_per_s']
                                        
                    v_sq_mps = vx_mps**2 + vy_mps**2
                    
                    ke_joules = 0.5 * M_PROTON_KG * v_sq_mps
                    ke_eV = ke_joules / E_CHARGE_C
                    final_energies_eV.append(ke_eV)

                # Acceleration calculation
                if len(df_traj) >= 2:
                    times = df_traj['time_s'].values
                    x_positions = df_traj['x_m'].values
                    vx_mps = df_traj['vx_m_per_s'].values
                    vy_mps = df_traj['vy_m_per_s'].values

                    for j in range(1, len(df_traj)):
                        dt = times[j] - times[j-1]
                        if dt > 0: # Avoid division by zero
                            accel_x = (vx_mps[j] - vx_mps[j-1]) / dt
                            accel_y = (vy_mps[j] - vy_mps[j-1]) / dt
                            accel_mag = np.sqrt(accel_x**2 + accel_y**2)
                            # Use x-position of the point where velocity is v[j]
                            current_x_m = x_positions[j] 
                            all_acceleration_data.append({'x_m': current_x_m, 'a_mag_mps2': accel_mag})
                
                # Velocity data extraction
                for idx, row in df_traj.iterrows():
                    v_mag = np.sqrt(row['vx_m_per_s']**2 + row['vy_m_per_s']**2)
                    all_velocity_data.append({'x_m': row['x_m'], 'v_mag_mps': v_mag})

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
    custom_handles = []
    if outline_threshold is not None: # Only add structure outline to legend if it was plotted
        custom_handles.append(structure_outline_proxy)
    
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

    # Plotting acceleration profile
    if all_acceleration_data:
        df_accel = pd.DataFrame(all_acceleration_data)
        df_accel.replace([np.inf, -np.inf], np.nan, inplace=True) # Handle potential infinities
        df_accel.dropna(subset=['a_mag_mps2'], inplace=True) # Remove NaNs

        if not df_accel.empty:
            num_bins = 100 # Number of bins along x-axis
            # Ensure x_bins cover the full range of x_coords for context, even if no accel data there
            min_x_coord_m = x_coords[0] * 1e-6 if x_coords else df_accel['x_m'].min()
            max_x_coord_m = x_coords[-1] * 1e-6 if x_coords else df_accel['x_m'].max()

            # Filter df_accel to be within the simulation's x_coords range before binning
            df_accel_filtered = df_accel[(df_accel['x_m'] >= min_x_coord_m) & (df_accel['x_m'] <= max_x_coord_m)].copy() # Use .copy() to avoid SettingWithCopyWarning
            
            if not df_accel_filtered.empty:
                x_bins = np.linspace(min_x_coord_m, max_x_coord_m, num_bins + 1)
                # Use .loc for assignment to ensure it works on the DataFrame/copy
                df_accel_filtered.loc[:, 'x_bin_group'] = pd.cut(df_accel_filtered['x_m'], bins=x_bins, include_lowest=True, right=True)
                
                binned_accel_stats = df_accel_filtered.groupby('x_bin_group', observed=False)['a_mag_mps2'].agg(['mean', 'std']).reset_index()
                binned_accel_stats['x_bin_center_m'] = binned_accel_stats['x_bin_group'].apply(lambda x: x.mid if isinstance(x, pd.Interval) else np.nan)
                binned_accel_stats.dropna(subset=['x_bin_center_m'], inplace=True) # Remove rows where center couldn't be calculated
                binned_accel_stats['std'].fillna(0, inplace=True) # Fill NaN std (e.g. for bins with 1 point)

                if not binned_accel_stats.empty:
                    fig_accel, ax_accel = plt.subplots(figsize=(12, 7))
                    
                    # Ensure x_bin_center_m is float before multiplication
                    x_plot_um = binned_accel_stats['x_bin_center_m'].astype(float) * 1e6 # Convert x to µm for plotting
                    mean_a_plot = binned_accel_stats['mean']
                    std_a_plot = binned_accel_stats['std']

                    ax_accel.plot(x_plot_um, mean_a_plot, color='dodgerblue', linestyle='-', linewidth=1.5, label='Mean Acceleration Magnitude')
                    ax_accel.fill_between(x_plot_um, mean_a_plot - std_a_plot, mean_a_plot + std_a_plot, color='lightskyblue', alpha=0.4, label='Std Dev (Fluctuation)')
                    
                    ax_accel.set_xlabel("X-position (µm)")
                    ax_accel.set_ylabel("Acceleration Magnitude (m/s²)")
                    ax_accel.set_title("Mean Proton Acceleration Magnitude vs. X-position")
                    ax_accel.legend(loc='upper right')
                    ax_accel.grid(True, linestyle=':', alpha=0.7)
                    ax_accel.set_xlim(0, L_total_sim_um) # Use simulation domain limits
                    # Optionally set y-limits if needed, e.g., ax_accel.set_ylim(bottom=0)
                    
                    plt.tight_layout()
                    try:
                        plt.savefig(output_accel_plot_file, dpi=300)
                        print(f"Acceleration profile plot saved to {output_accel_plot_file}")
                    except Exception as e:
                        print(f"Error saving acceleration profile plot: {e}")
                    plt.show()
                    plt.close(fig_accel)
                else:
                    print("No valid binned acceleration data to plot.")
            else:
                print("No acceleration data points fall within the simulation's x-coordinate range.")
        else:
            print("No valid acceleration data to plot after cleaning.")
    else:
        print("No acceleration data calculated to generate a profile plot.")

    # Plotting velocity profile
    if all_velocity_data:
        df_vel = pd.DataFrame(all_velocity_data)
        df_vel.replace([np.inf, -np.inf], np.nan, inplace=True) # Handle potential infinities
        df_vel.dropna(subset=['v_mag_mps'], inplace=True) # Remove NaNs

        if not df_vel.empty:
            num_bins_vel = 100 # Number of bins along x-axis for velocity
            min_x_coord_m_vel = x_coords[0] * 1e-6 if x_coords else df_vel['x_m'].min()
            max_x_coord_m_vel = x_coords[-1] * 1e-6 if x_coords else df_vel['x_m'].max()

            df_vel_filtered = df_vel[(df_vel['x_m'] >= min_x_coord_m_vel) & (df_vel['x_m'] <= max_x_coord_m_vel)].copy()

            if not df_vel_filtered.empty:
                x_bins_vel = np.linspace(min_x_coord_m_vel, max_x_coord_m_vel, num_bins_vel + 1)
                df_vel_filtered.loc[:, 'x_bin_group_vel'] = pd.cut(df_vel_filtered['x_m'], bins=x_bins_vel, include_lowest=True, right=True)
                
                binned_vel_stats = df_vel_filtered.groupby('x_bin_group_vel', observed=False)['v_mag_mps'].agg(['mean', 'std']).reset_index()
                binned_vel_stats['x_bin_center_m_vel'] = binned_vel_stats['x_bin_group_vel'].apply(lambda x: x.mid if isinstance(x, pd.Interval) else np.nan)
                binned_vel_stats.dropna(subset=['x_bin_center_m_vel'], inplace=True)
                binned_vel_stats['std'].fillna(0, inplace=True)

                if not binned_vel_stats.empty:
                    fig_vel, ax_vel = plt.subplots(figsize=(12, 7))
                    
                    x_plot_um_vel = binned_vel_stats['x_bin_center_m_vel'].astype(float) * 1e6
                    mean_v_plot = binned_vel_stats['mean']
                    std_v_plot = binned_vel_stats['std']

                    ax_vel.plot(x_plot_um_vel, mean_v_plot, color='green', linestyle='-', linewidth=1.5, label='Mean Velocity Magnitude')
                    ax_vel.fill_between(x_plot_um_vel, mean_v_plot - std_v_plot, mean_v_plot + std_v_plot, color='lightgreen', alpha=0.4, label='Std Dev (Fluctuation)')
                    
                    ax_vel.set_xlabel("X-position (µm)")
                    ax_vel.set_ylabel("Velocity Magnitude (m/s)")
                    ax_vel.set_title("Mean Proton Velocity Magnitude vs. X-position")
                    ax_vel.legend(loc='upper right')
                    ax_vel.grid(True, linestyle=':', alpha=0.7)
                    ax_vel.set_xlim(0, L_total_sim_um)
                    ax_vel.set_ylim(bottom=0) # Velocity magnitude is non-negative
                    
                    plt.tight_layout()
                    try:
                        plt.savefig(output_vel_plot_file, dpi=300)
                        print(f"Velocity profile plot saved to {output_vel_plot_file}")
                    except Exception as e:
                        print(f"Error saving velocity profile plot: {e}")
                    plt.show()
                    plt.close(fig_vel)
                else:
                    print("No valid binned velocity data to plot.")
            else:
                print("No velocity data points fall within the simulation's x-coordinate range.")
        else:
            print("No valid velocity data to plot after cleaning.")
    else:
        print("No velocity data calculated to generate a profile plot.")


if __name__ == "__main__":
    cli_folder_path = None
    if len(sys.argv) > 1:
        cli_folder_path = sys.argv[1]
    main(folder_path=cli_folder_path)
