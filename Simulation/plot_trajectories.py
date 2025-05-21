import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os
import glob
import numpy as np # For NaN comparison if needed, and math.

# Physical constants
M_PROTON_KG = 1.67262192e-27  # Mass of proton in kg
E_CHARGE_C = 1.60217663e-19   # Elementary charge in Coulombs

def load_geometry_params(filepath):
    """Loads geometry parameters from a CSV file (values expected in µm)."""
    params = {}
    try:
        df = pd.read_csv(filepath, header=None, index_col=0)
        params = df[1].to_dict()
    except FileNotFoundError:
        print(f"Error: Geometry parameters file not found at {filepath}")
        return None
    except Exception as e:
        print(f"Error reading geometry parameters: {e}")
        return None
    # Ensure numeric types
    for key in params:
        try:
            params[key] = float(params[key])
        except ValueError:
            print(f"Warning: Could not convert param {key} to float.")
            pass # Keep as string if not convertible, though all expected are float
    return params

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

def main():
    input_base_folder = "geometria_piana"
    trajectories_folder = os.path.join(input_base_folder, "proton_trajectories")
    output_plot_file = os.path.join(input_base_folder, "proton_trajectories_plot.png")
    output_hist_plot_file = os.path.join(input_base_folder, "proton_final_energy_histogram.png")

    geom_params_file = os.path.join(input_base_folder, "geometry_params.csv")
    x_coords_file = os.path.join(input_base_folder, "x_coordinates.csv")
    y_coords_file = os.path.join(input_base_folder, "y_coordinates.csv")

    geom = load_geometry_params(geom_params_file) # geom dimensions are in µm
    x_coords = load_coordinates(x_coords_file)   # x_coords are in µm
    y_coords = load_coordinates(y_coords_file)   # y_coords are in µm

    if geom is None or x_coords is None or y_coords is None:
        print("Could not load necessary data. Exiting.")
        return

    L_total_sim_um = x_coords[-1] if x_coords else 0 # µm
    H_total_sim_um = y_coords[-1] if y_coords else geom.get('H_total', 0) # µm
    grid_spacing_h_um = geom.get('h', 0.1) # µm, for tolerance

    fig, ax = plt.subplots(figsize=(12, 8))

    # Plot geometry boundaries (using µm values directly)
    x_fs_um = geom.get('x_free_space', 0)
    x_sl_um = geom.get('x_structure_len', 0)
    y_slt_um = geom.get('y_si_layer_thick', 0)
    y_vgt_um = geom.get('y_vacuum_gap_thick', 0)
    
    # Bottom Silicon Layer
    ax.add_patch(patches.Rectangle((x_fs_um, 0), x_sl_um, y_slt_um, facecolor='lightblue', edgecolor='blue', label='Bottom Si Layer'))
    
    H_structure_total_um = geom.get('H_total', H_total_sim_um) 
    top_si_height_um = H_structure_total_um - (y_slt_um + y_vgt_um)
    if top_si_height_um > 0:
         ax.add_patch(patches.Rectangle((x_fs_um, y_slt_um + y_vgt_um), x_sl_um, top_si_height_um, facecolor='lightblue', edgecolor='blue', label='Top Si Layer'))
    else:
        print("Warning: Top Si layer has zero or negative height based on parameters.")

    final_energies_eV = []
    successful_protons_count = 0

    trajectory_files = glob.glob(os.path.join(trajectories_folder, "proton_*_trajectory.csv"))
    
    if not trajectory_files:
        print(f"No trajectory files found in {trajectories_folder}")
    
    num_trajectories_to_plot = min(len(trajectory_files), 100) 
    
    for i, traj_file in enumerate(trajectory_files):
        # Process all files for energy histogram, plot only up to num_trajectories_to_plot
        try:
            # Trajectory data is in SI units (m, m/s)
            df_traj = pd.read_csv(traj_file) 
            if not df_traj.empty:
                if i < num_trajectories_to_plot: 
                    # Convert m to µm for plotting
                    ax.plot(df_traj['x_m'] * 1e6, df_traj['y_m'] * 1e6, linestyle='-', linewidth=0.5, alpha=0.7)
                
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
                    
        except pd.errors.EmptyDataError:
            print(f"Warning: Trajectory file {traj_file} is empty.")
        except Exception as e:
            print(f"Error reading or plotting trajectory file {traj_file}: {e}")

    ax.set_xlabel("X (µm)") # Axis labels remain in µm for visualization
    ax.set_ylabel("Y (µm)") # Axis labels remain in µm for visualization
    ax.set_title(f"Proton Trajectories ({min(len(trajectory_files), num_trajectories_to_plot)} shown, plotted in µm)")
    
    # Set plot limits based on the simulation domain (in µm)
    ax.set_xlim(0, L_total_sim_um)
    ax.set_ylim(0, H_total_sim_um)
    
    ax.legend(loc='upper right')
    ax.set_aspect('equal', adjustable='box') # Important for correct geometric representation
    plt.grid(True, linestyle=':', alpha=0.7)
    
    try:
        plt.savefig(output_plot_file, dpi=300)
        print(f"Plot saved to {output_plot_file}")
    except Exception as e:
        print(f"Error saving plot: {e}")
    
    plt.show() # Show trajectory plot
    plt.close(fig) # Close the trajectory plot figure

    # Plotting the histogram of final energies
    if final_energies_eV:
        print(f"Number of protons considered successful for energy histogram: {successful_protons_count}")
        fig_hist, ax_hist = plt.subplots(figsize=(10, 6))
        ax_hist.hist(final_energies_eV, bins=50, edgecolor='black', alpha=0.75)
        ax_hist.set_xlabel("Final Kinetic Energy (eV)")
        ax_hist.set_ylabel("Number of Protons")
        ax_hist.set_title("Histogram of Final Kinetic Energies of Successful Protons")
        ax_hist.grid(True, linestyle=':', alpha=0.7)
        
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
    main()
