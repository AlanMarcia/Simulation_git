import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os
import glob
import numpy as np # For NaN comparison if needed, and math.
import sys # Import sys module for command-line arguments
from functools import lru_cache  # For caching expensive operations
import gc  # Garbage collection for managing memory

# Set matplotlib to use a faster non-interactive backend
# Use 'Agg' backend for better performance when generating plots to files
plt.switch_backend('Agg')

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
        # Use NumPy directly for better performance
        # This avoids the overhead of creating a pandas DataFrame
        coords = np.loadtxt(filepath, delimiter=',')
        # If needed, convert to float32 to save memory
        if np.all(np.abs(coords) < 1e6):  # Check if values are within a reasonable range
            coords = coords.astype(np.float32)
        return coords.tolist()
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
        # Use optimized loading parameters
        data = np.loadtxt(filename, delimiter=',', dtype=np.float64)
        # Memory optimization: Check if data can be downcast to float32 without significant loss
        if np.all(np.abs(data) < 1e6):  # Check if values are within a reasonable range
            data = data.astype(np.float32)  # Use float32 instead of float64 to save memory
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
    output_angle_hist_file = os.path.join(input_base_folder, "proton_final_angle_histogram.png") # New histogram for angles
    output_accel_plot_file = os.path.join(input_base_folder, "proton_acceleration_profile.png") # New plot file
    output_vel_plot_file = os.path.join(input_base_folder, "proton_velocity_profile.png") # New plot file for velocity

    # geom_params_file = os.path.join(input_base_folder, "geometry_params.csv") # File still exists, but not used for plotting shapes
    eps_r_file = os.path.join(input_base_folder, "permittivity.csv")
    x_coords_file = os.path.join(input_base_folder, "x_coordinates.csv")
    y_coords_file = os.path.join(input_base_folder, "y_coordinates.csv")

    # geom = load_geometry_params(geom_params_file) # geom dimensions are in µm - No longer used for plotting shapes
    eps_r_data = load_2d_csv(eps_r_file) # Expected shape (Ny, Nx)
    x_coords = load_coordinates(x_coords_file)   # x_coords are in μm
    y_coords = load_coordinates(y_coords_file)   # y_coords are in μm

    if eps_r_data is None or x_coords is None or y_coords is None:
        print("Could not load necessary data (permittivity, coordinates). Exiting.")
        return

    L_total_sim_um = x_coords[-1] if x_coords else 0 # μm
    H_total_sim_um = y_coords[-1] if y_coords else 0 # μm
    
    # Default grid spacing for tolerance, as geom_params is not loaded for 'h'
    # This value should ideally match the 'h' used in the simulation for accurate success check.
    # If geometry_params.csv is available and 'h' is needed, it could be loaded separately for just this value.
    grid_spacing_h_um = 0.5 # μm, default for tolerance, adjust if needed

    fig, ax = plt.subplots(figsize=(12, 8))

    # Plot geometry outlines using permittivity map
    X_mesh, Y_mesh = np.meshgrid(x_coords, y_coords)
      # Dynamically calculate outline_threshold from eps_r_data - optimized
    outline_threshold = None
    if eps_r_data is not None:
        # Use more efficient approach to find min and max values
        eps_min = np.min(eps_r_data)
        eps_max = np.max(eps_r_data)
        
        # Only compute unique values if needed for reporting
        if eps_min != eps_max and eps_max > eps_min + 0.5:
            outline_threshold = (eps_min + eps_max) / 2.0
            print(f"Dynamically calculated outline threshold for trajectories plot: {outline_threshold:.2f} (from eps_r min/max: {eps_min:.2f}, {eps_max:.2f})")
        elif eps_min == eps_max:
            print(f"Warning (trajectories plot): Permittivity data contains only one unique value ({eps_min:.2f}). No outlines will be drawn.")
        else:
            print(f"Warning (trajectories plot): Permittivity values ({eps_min:.2f}-{eps_max:.2f}) are too close. Could not reliably determine outline threshold.")
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
    final_angles_deg = []  # To store final trajectory angles for successful protons

    # trajectory_files = glob.glob(os.path.join(trajectories_folder, "proton_*_trajectory.csv")) # Removed
      # Load the single trajectory file
    try:
        # Use chunksize for better memory efficiency when dealing with very large files
        # Only read columns we need
        required_columns = ['proton_id', 'time_s', 'x_m', 'y_m', 'vx_m_per_s', 'vy_m_per_s']
        df_all_trajectories = pd.read_csv(all_trajectories_file, usecols=required_columns)
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
        
        num_trajectories_to_plot = min(num_protons_in_file, 1000)        # Plot a subset of trajectories - optimized
        # Pre-calculate conversion factor from m to μm
        convert_to_um = 1e6
        
        # Use batch processing for trajectory plotting
        batch_size = 20  # Process trajectories in batches for better performance
        legend_added = False
        
        for batch_start in range(0, min(num_trajectories_to_plot, len(proton_ids)), batch_size):
            batch_end = min(batch_start + batch_size, num_trajectories_to_plot, len(proton_ids))
            batch_ids = proton_ids[batch_start:batch_end]
            
            for i, proton_id in enumerate(batch_ids):
                df_traj = grouped_trajectories.get_group(proton_id)
                if not df_traj.empty:
                    # Use numpy arrays directly for better performance
                    # Note: x_m and y_m are already in meters from proton_simulator.cpp
                    x_vals = df_traj['x_m'].values * convert_to_um  # Convert from m to μm for plotting
                    y_vals = df_traj['y_m'].values * convert_to_um  # Convert from m to μm for plotting
                    
                    line, = ax.plot(x_vals, y_vals, linestyle='-', linewidth=0.5, 
                                   alpha=0.7, color='red')
                    
                    if not plotted_traj_legend_added:
                        line.set_label('Proton Trajectories') # Add label to one line for the legend
                        plotted_traj_legend_added = True
            
            # Force garbage collection after each batch
            gc.collect()
          # Process all trajectories for energy histogram and acceleration - optimized
        # Pre-allocate lists for better performance
        max_protons = len(proton_ids)
        acceleration_data_capacity = max_protons * 100  # Estimate average number of points per trajectory
        velocity_data_capacity = max_protons * 100
        
        all_acceleration_data = []
        all_velocity_data = []
        
        # Convert threshold once - L_total_sim_um is in μm, convert to meters for comparison
        x_success_threshold = (L_total_sim_um * 1e-6) - (grid_spacing_h_um * 1e-6 / 2.0)
        
        for proton_id in proton_ids:
            df_traj = grouped_trajectories.get_group(proton_id)
            if not df_traj.empty:
                # Energy calculation - use vectorized operations and avoid iloc for better performance
                # Note: x_m and velocities are already in SI units (m, m/s) from proton_simulator.cpp
                last_point = df_traj.iloc[-1]
                if last_point['x_m'] >= x_success_threshold:
                    successful_protons_count += 1
                    # Vectorized calculation for energy - all units are already SI
                    vx_mps = last_point['vx_m_per_s']  # Already in m/s
                    vy_mps = last_point['vy_m_per_s']  # Already in m/s
                    v_sq_mps = vx_mps**2 + vy_mps**2
                    ke_joules = 0.5 * M_PROTON_KG * v_sq_mps  # Correct SI calculation
                    ke_eV = ke_joules / E_CHARGE_C
                    final_energies_eV.append(ke_eV)
                    
                    # Calculate final trajectory angle (in degrees)
                    final_angle_rad = np.arctan2(vy_mps, vx_mps)
                    final_angle_deg = np.degrees(final_angle_rad)
                    final_angles_deg.append(final_angle_deg)

                # Acceleration calculation - vectorized approach
                if len(df_traj) >= 2:
                    # Extract all necessary arrays at once - all already in SI units
                    times = df_traj['time_s'].values        # seconds
                    x_positions = df_traj['x_m'].values     # meters (converted from μm by proton_simulator.cpp)
                    vx_mps = df_traj['vx_m_per_s'].values   # m/s
                    vy_mps = df_traj['vy_m_per_s'].values   # m/s
                    
                    # Calculate time differences
                    dt_values = np.diff(times)
                    # Filter out zero time differences to avoid division by zero
                    valid_indices = np.where(dt_values > 0)[0]
                    
                    if len(valid_indices) > 0:
                        # Calculate acceleration vectors using numpy operations
                        accel_x_values = np.diff(vx_mps)[valid_indices] / dt_values[valid_indices]
                        accel_y_values = np.diff(vy_mps)[valid_indices] / dt_values[valid_indices]
                        accel_mag_values = np.sqrt(accel_x_values**2 + accel_y_values**2)
                        
                        # Corresponding x-positions (at j+1)
                        x_pos_values = x_positions[valid_indices + 1]
                        
                        # Extend all_acceleration_data efficiently
                        for i in range(len(valid_indices)):
                            all_acceleration_data.append({
                                'x_m': x_pos_values[i],
                                'a_mag_mps2': accel_mag_values[i]
                            })
                
                # Velocity data extraction - vectorized approach
                # Calculate velocity magnitudes for all points at once
                v_mag_values = np.sqrt(df_traj['vx_m_per_s'].values**2 + df_traj['vy_m_per_s'].values**2)
                
                # Create and extend all_velocity_data efficiently
                for i, x_m in enumerate(df_traj['x_m'].values):
                    all_velocity_data.append({
                        'x_m': x_m, 
                        'v_mag_mps': v_mag_values[i]
                    })

    else:
        print(f"No trajectory data loaded from {all_trajectories_file}. Skipping trajectory plotting and energy analysis.")
        num_trajectories_to_plot = 0

    # --- Print Summary Statistics ---
    print("\n" + "="*60)
    print("PROTON TRAJECTORY ANALYSIS SUMMARY")
    print("="*60)
    print(f"Total protons simulated: {num_protons_in_file}")
    print(f"Successful protons (reached end): {successful_protons_count}")
    if num_protons_in_file > 0:
        efficiency = (successful_protons_count / num_protons_in_file) * 100
        print(f"Transmission efficiency: {efficiency:.2f}%")
    
    if final_energies_eV:
        print(f"\nFinal Energy Statistics (successful protons):")
        print(f"  Mean energy: {np.mean(final_energies_eV):.2f} eV")
        print(f"  Std deviation: {np.std(final_energies_eV):.2f} eV")
        print(f"  Min energy: {np.min(final_energies_eV):.2f} eV")
        print(f"  Max energy: {np.max(final_energies_eV):.2f} eV")
        print(f"  Median energy: {np.median(final_energies_eV):.2f} eV")
    
    if final_angles_deg:
        print(f"\nFinal Trajectory Angle Statistics (successful protons):")
        print(f"  Mean angle: {np.mean(final_angles_deg):.2f}°")
        print(f"  Std deviation: {np.std(final_angles_deg):.2f}°")
        print(f"  Min angle: {np.min(final_angles_deg):.2f}°")
        print(f"  Max angle: {np.max(final_angles_deg):.2f}°")
    print("="*60 + "\n")


    ax.set_xlabel("X (μm)") # Axis labels remain in μm for visualization
    ax.set_ylabel("Y (μm)") # Axis labels remain in μm for visualization
    ax.set_title(f"Proton Trajectories ({num_trajectories_to_plot} of {num_protons_in_file} shown, plotted in μm)")
    
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
    
    # Plotting the histogram of final trajectory angles
    if final_angles_deg:
        print(f"Generating angle distribution histogram for {len(final_angles_deg)} successful protons...")
        fig_angle, ax_angle = plt.subplots(figsize=(10, 6))
        ax_angle.hist(final_angles_deg, bins=50, edgecolor='black', alpha=0.75, color='coral')
        ax_angle.set_xlabel("Final Trajectory Angle (degrees)")
        ax_angle.set_ylabel("Number of Protons")
        ax_angle.set_title("Histogram of Final Trajectory Angles of Successful Protons")
        ax_angle.grid(True, linestyle=':', alpha=0.7)
        
        # Calculate and display angle statistics
        mean_angle = np.mean(final_angles_deg)
        std_angle = np.std(final_angles_deg)
        min_angle = np.min(final_angles_deg)
        max_angle = np.max(final_angles_deg)
        median_angle = np.median(final_angles_deg)
        
        stats_text_angle = f'Mean: {mean_angle:.2f}°\nStd Dev: {std_angle:.2f}°\nMedian: {median_angle:.2f}°\nMin: {min_angle:.2f}°\nMax: {max_angle:.2f}°'
        
        ax_angle.text(0.95, 0.95, stats_text_angle, transform=ax_angle.transAxes, fontsize=9,
                     verticalalignment='top', horizontalalignment='right',
                     bbox=dict(boxstyle='round,pad=0.5', fc='lightyellow', alpha=0.5))
        
        # Add a vertical line at 0 degrees for reference
        ax_angle.axvline(x=0, color='red', linestyle='--', linewidth=1.5, alpha=0.7, label='0° (horizontal)')
        ax_angle.legend(loc='upper left')
        
        plt.tight_layout()
        try:
            plt.savefig(output_angle_hist_file, dpi=300)
            print(f"Angle distribution histogram saved to {output_angle_hist_file}")
        except Exception as e:
            print(f"Error saving angle histogram: {e}")
        plt.show()
        plt.close(fig_angle)
    else:
        print("No successful protons found to generate an angle histogram.")
        # Plotting acceleration profile - optimized
    if all_acceleration_data:
        # Convert to numpy arrays directly for faster processing
        # Note: x_m values are in meters, a_mag_mps2 are in m/s²
        x_vals = np.array([item['x_m'] for item in all_acceleration_data])  # meters
        a_vals = np.array([item['a_mag_mps2'] for item in all_acceleration_data])  # m/s²
        
        # Handle infinities and NaNs - vectorized operation
        valid_indices = np.isfinite(a_vals)
        x_vals = x_vals[valid_indices]
        a_vals = a_vals[valid_indices]
        
        if len(x_vals) > 0:
            num_bins = 100 # Number of bins along x-axis
            # Ensure x_bins cover the full range of x_coords for context
            # x_coords are in μm, convert to meters
            min_x_coord_m = x_coords[0] * 1e-6 if x_coords else np.min(x_vals)
            max_x_coord_m = x_coords[-1] * 1e-6 if x_coords else np.max(x_vals)

            # Filter data to be within the simulation's x_coords range before binning
            valid_range_indices = (x_vals >= min_x_coord_m) & (x_vals <= max_x_coord_m)
            x_vals_filtered = x_vals[valid_range_indices]
            a_vals_filtered = a_vals[valid_range_indices]
            
            if len(x_vals_filtered) > 0:
                # Use numpy's histogram function for binning - much faster than pandas cut+groupby
                x_bins = np.linspace(min_x_coord_m, max_x_coord_m, num_bins + 1)
                bin_indices = np.digitize(x_vals_filtered, x_bins) - 1
                
                # Ensure bin indices are within valid range
                valid_bin_indices = (bin_indices >= 0) & (bin_indices < len(x_bins) - 1)
                x_vals_filtered = x_vals_filtered[valid_bin_indices]
                a_vals_filtered = a_vals_filtered[valid_bin_indices]
                bin_indices = bin_indices[valid_bin_indices]
                
                # Calculate bin centers
                bin_centers = (x_bins[:-1] + x_bins[1:]) / 2
                
                # Calculate statistics per bin
                mean_vals = np.zeros(num_bins)
                std_vals = np.zeros(num_bins)
                count_vals = np.zeros(num_bins, dtype=int)
                
                for i in range(num_bins):
                    bin_data = a_vals_filtered[bin_indices == i]
                    if len(bin_data) > 0:
                        mean_vals[i] = np.mean(bin_data)
                        std_vals[i] = np.std(bin_data) if len(bin_data) > 1 else 0
                        count_vals[i] = len(bin_data)
                  # Filter out empty bins
                valid_bins = count_vals > 0
                bin_centers = bin_centers[valid_bins]
                mean_vals = mean_vals[valid_bins]
                std_vals = std_vals[valid_bins]
                
                if len(bin_centers) > 0:
                    fig_accel, ax_accel = plt.subplots(figsize=(12, 7))
                    
                    # Convert bin centers to µm for plotting
                    x_plot_um = bin_centers * 1e6
                    
                    # Use optimized plotting approach
                    ax_accel.plot(x_plot_um, mean_vals, color='dodgerblue', linestyle='-', linewidth=1.5, label='Mean Acceleration Magnitude')
                    ax_accel.fill_between(x_plot_um, mean_vals - std_vals, mean_vals + std_vals, color='lightskyblue', alpha=0.4, label='Std Dev (Fluctuation)')
                    
                    ax_accel.set_xlabel("X-position (μm)")
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
        print("No acceleration data calculated to generate a profile plot.")    # Plotting velocity profile - optimized using numpy
    if all_velocity_data:
        # Convert to numpy arrays directly for faster processing
        # Note: x_m values are in meters, v_mag_mps are in m/s
        x_vals_vel = np.array([item['x_m'] for item in all_velocity_data])  # meters
        v_vals = np.array([item['v_mag_mps'] for item in all_velocity_data])  # m/s
        
        # Handle infinities and NaNs in one step
        valid_indices_vel = np.isfinite(v_vals)
        x_vals_vel = x_vals_vel[valid_indices_vel]
        v_vals = v_vals[valid_indices_vel]
        
        if len(x_vals_vel) > 0:
            num_bins_vel = 100 # Number of bins along x-axis for velocity
            # x_coords are in μm, convert to meters for comparison
            min_x_coord_m_vel = x_coords[0] * 1e-6 if x_coords else np.min(x_vals_vel)
            max_x_coord_m_vel = x_coords[-1] * 1e-6 if x_coords else np.max(x_vals_vel)

            # Filter data to be within simulation bounds
            valid_range_indices_vel = (x_vals_vel >= min_x_coord_m_vel) & (x_vals_vel <= max_x_coord_m_vel)
            x_vals_vel_filtered = x_vals_vel[valid_range_indices_vel]
            v_vals_filtered = v_vals[valid_range_indices_vel]
            
            if len(x_vals_vel_filtered) > 0:
                # Use numpy's histogram function for binning
                x_bins_vel = np.linspace(min_x_coord_m_vel, max_x_coord_m_vel, num_bins_vel + 1)
                bin_indices_vel = np.digitize(x_vals_vel_filtered, x_bins_vel) - 1
                
                # Ensure bin indices are within valid range
                valid_bin_indices_vel = (bin_indices_vel >= 0) & (bin_indices_vel < len(x_bins_vel) - 1)
                x_vals_vel_filtered = x_vals_vel_filtered[valid_bin_indices_vel]
                v_vals_filtered = v_vals_filtered[valid_bin_indices_vel]
                bin_indices_vel = bin_indices_vel[valid_bin_indices_vel]
                
                # Calculate bin centers
                bin_centers_vel = (x_bins_vel[:-1] + x_bins_vel[1:]) / 2
                
                # Calculate statistics per bin - optimized approach
                mean_vals_vel = np.zeros(num_bins_vel)
                std_vals_vel = np.zeros(num_bins_vel)
                count_vals_vel = np.zeros(num_bins_vel, dtype=int)
                
                for i in range(num_bins_vel):
                    bin_data = v_vals_filtered[bin_indices_vel == i]
                    if len(bin_data) > 0:
                        mean_vals_vel[i] = np.mean(bin_data)
                        std_vals_vel[i] = np.std(bin_data) if len(bin_data) > 1 else 0
                        count_vals_vel[i] = len(bin_data)
                
                # Filter out empty bins
                valid_bins_vel = count_vals_vel > 0
                bin_centers_vel = bin_centers_vel[valid_bins_vel]
                mean_vals_vel = mean_vals_vel[valid_bins_vel]
                std_vals_vel = std_vals_vel[valid_bins_vel]

                if len(bin_centers_vel) > 0:
                    fig_vel, ax_vel = plt.subplots(figsize=(12, 7))
                    
                    x_plot_um_vel = bin_centers_vel * 1e6  # Convert from meters to μm for plotting
                    
                    # Efficient plotting
                    ax_vel.plot(x_plot_um_vel, mean_vals_vel, color='green', linestyle='-', linewidth=1.5, label='Mean Velocity Magnitude')
                    ax_vel.fill_between(x_plot_um_vel, mean_vals_vel - std_vals_vel, mean_vals_vel + std_vals_vel, color='lightgreen', alpha=0.4, label='Std Dev (Fluctuation)')
                    
                    ax_vel.set_xlabel("X-position (μm)")
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
        print("No velocity data calculated to generate a profile plot.")    # --- Create and Save Animation of First 10 Protons (optimized) ---
    try:
        print("Preparing animation of first 10 protons...")
        NUM_PARTICLES = 10
        FPS = 5
        OUTPUT_VIDEO = os.path.join(input_base_folder, "particle_motion.mp4")

        # Check for animation module availability before processing data
        animation_module_available = False
        try:
            from matplotlib.animation import FuncAnimation
            animation_module_available = True
        except ImportError:
            print("Warning: matplotlib.animation not available. Skipping animation creation.")
            
        # Only proceed if animation module is available
        if animation_module_available and df_all_trajectories is not None and not df_all_trajectories.empty:
            # Filter data first to reduce memory usage
            particle_ids = df_all_trajectories['proton_id'].unique()[:NUM_PARTICLES]
            
            # Use query for more efficient filtering
            query_str = ' | '.join([f'proton_id == {pid}' for pid in particle_ids])
            anim_df = df_all_trajectories.query(query_str).copy()
            
            # Get time values and limits only once
            times = np.sort(anim_df['time_s'].unique())
            xlim = (anim_df['x_m'].min(), anim_df['x_m'].max())
            ylim = (anim_df['y_m'].min(), anim_df['y_m'].max())

            # Set up plot once
            fig_anim, ax_anim = plt.subplots(figsize=(10, 6), dpi=100)  # Lower dpi for animation performance
            scat = ax_anim.scatter([], [], s=20)
            ax_anim.set_xlim(xlim)
            ax_anim.set_ylim(ylim)
            ax_anim.set_xlabel("X (m)")
            ax_anim.set_ylabel("Y (m)")
            ax_anim.set_title(f"Motion of {NUM_PARTICLES} Protons")

            # Add structure contour if available - do this once
            if eps_r_data is not None and outline_threshold is not None:
                # Convert coordinates efficiently using broadcasting
                # x_coords and y_coords are in μm, convert to meters for mesh
                x_coords_m = np.array(x_coords)[:, np.newaxis] * 1e-6  # Add dimension for broadcasting
                y_coords_m = np.array(y_coords)[np.newaxis, :] * 1e-6  # Add dimension for broadcasting
                X_mesh_anim, Y_mesh_anim = np.meshgrid(x_coords_m.flatten(), y_coords_m.flatten())
                
                # Use a simpler contour for animation (fewer levels, simpler lines)
                ax_anim.contour(X_mesh_anim, Y_mesh_anim, eps_r_data, 
                               levels=[outline_threshold], colors='blue', 
                               linewidths=0.8, linestyles='--')            # Pre-compute particle positions for all time steps
            # This is a major optimization to avoid doing lookups in every animation frame
            print("Pre-computing particle positions for animation...")
            particle_positions = {}
            
            # Create a dictionary for faster lookups
            grouped_dict = {pid: anim_df[anim_df['proton_id'] == pid] for pid in particle_ids}
            
            # Create a time-sorted version of positions for each particle
            for t in times:
                xs, ys = [], []
                for pid in particle_ids:
                    group_df = grouped_dict[pid]
                    # Use binary search for finding the right time index (much faster)
                    # Find the index of the largest time value <= t
                    idx = np.searchsorted(group_df['time_s'].values, t, side='right') - 1
                    if idx >= 0:  # If a valid position exists
                        xs.append(group_df['x_m'].values[idx])
                        ys.append(group_df['y_m'].values[idx])
                    else:
                        # If no valid position yet, use starting position or skip
                        if len(group_df) > 0:
                            xs.append(group_df['x_m'].values[0])
                            ys.append(group_df['y_m'].values[0])
                
                # Store the positions for this time step
                particle_positions[t] = (xs, ys)
            
            print(f"Position data prepared for {len(times)} time steps")
            
            def update(frame):
                t = times[frame]
                # Get pre-computed positions
                xs, ys = particle_positions[t]
                scat.set_offsets(list(zip(xs, ys)))
                ax_anim.set_title(f"Time: {t:.2e} s")
                return scat,            # Use a more optimized approach for animation
            ani = None
            if len(times) > 1:
                from matplotlib.animation import FuncAnimation
                
                # Subsample frames if there are too many (for performance)
                max_frames = 200  # Limit number of frames for performance
                frame_indices = np.linspace(0, len(times)-1, min(max_frames, len(times)), dtype=int)
                selected_times = times[frame_indices]
                
                # Create animation with subsampled frames
                ani = FuncAnimation(
                    fig_anim, update, frames=frame_indices, 
                    interval=1000 / FPS, blit=True)
                
                print("Saving video...")
                # Use a more optimized writer
                try:
                    # Try using ffmpeg writer with optimized settings
                    from matplotlib.animation import FFMpegWriter
                    writer = FFMpegWriter(fps=FPS, metadata=dict(artist='Physics Simulation'),
                                         bitrate=1800)  # Lower bitrate for smaller file size
                    ani.save(OUTPUT_VIDEO, writer=writer)
                except (ImportError, ValueError):
                    # Fall back to default writer if FFMpegWriter not available
                    ani.save(OUTPUT_VIDEO, fps=FPS, extra_args=['-vcodec', 'libx264', 
                                                            '-pix_fmt', 'yuv420p',
                                                            '-profile:v', 'baseline'])
                print(f"Video saved as {OUTPUT_VIDEO}")
            else:
                print("Not enough time steps for animation.")
                
            # Explicitly clean up to free memory
            plt.close(fig_anim)
            del ani  # Explicitly delete the animation object
            if 'particle_positions' in locals():
                del particle_positions  # Clean up the pre-computed positions
    except Exception as e:
        print(f"Error during animation creation: {e}")

if __name__ == "__main__":
    cli_folder_path = None
    if len(sys.argv) > 1:
        cli_folder_path = sys.argv[1]
    main(folder_path=cli_folder_path)
