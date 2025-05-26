import numpy as np
import matplotlib.pyplot as plt
import os
import csv
import glob # For finding trajectory files
import sys # Import sys module for command-line arguments

# Constants for energy calculation
M_PROTON = 1.67262192e-27  # kg
Q_ELECTRON = 1.60217663e-19 # Coulombs, for eV conversion

def load_1d_csv(filename):
    if not os.path.exists(filename):
        print(f"Error: File {filename} not found.")
        return None
    return np.loadtxt(filename, delimiter=',')

def load_2d_csv(filename):
    if not os.path.exists(filename):
        print(f"Error: File {filename} not found.")
        return None
    data = np.loadtxt(filename, delimiter=',')
    return data # Remove .T to keep shape (Ny, Nx)

def load_geometry_params_for_plot(filepath):
    """Loads specific geometry parameters from a CSV file for plotting purposes."""
    params = {}
    if not os.path.exists(filepath):
        print(f"Error: Geometry parameters file {filepath} not found.")
        return None
    try:
        with open(filepath, mode='r') as infile:
            reader = csv.reader(infile)
            for rows in reader:
                if len(rows) == 2:
                    key = rows[0]
                    try:
                        params[key] = float(rows[1])
                    except ValueError:
                        print(f"Warning: Could not convert value for key '{key}' to float.")
                        params[key] = rows[1] # Store as string if not float
    except Exception as e:
        print(f"Error reading geometry parameters file {filepath}: {e}")
        return None
    return params

def plot_results(folder_path=None): # Add folder_path argument
    # Determine output folder
    if folder_path:
        output_folder_name = folder_path
        print(f"Using specified folder: {output_folder_name}")
    else:
        output_folder_name = "geometria_Denti_sfasati_profondi_5um" # Default folder
        print(f"No folder specified, using default: {output_folder_name}")
    
    if not os.path.isdir(output_folder_name):
        print(f"Error: The specified folder '{output_folder_name}' does not exist or is not a directory.")
        print("Please ensure the C++ simulation has run and created this folder with CSV files.")
        return None # Return None if folder doesn't exist

    # File paths are now relative to output_folder_name
    potential_file = os.path.join(output_folder_name, "potential.csv")
    ex_file = os.path.join(output_folder_name, "electric_field_x.csv")
    ey_file = os.path.join(output_folder_name, "electric_field_y.csv")
    eps_r_file = os.path.join(output_folder_name, "permittivity.csv")
    x_coords_file = os.path.join(output_folder_name, "x_coordinates.csv")
    y_coords_file = os.path.join(output_folder_name, "y_coordinates.csv")
    geometry_params_file = os.path.join(output_folder_name, "geometry_params.csv")

    # Load data
    V = load_2d_csv(potential_file)
    Ex = load_2d_csv(ex_file)
    Ey = load_2d_csv(ey_file)
    eps_r = load_2d_csv(eps_r_file)
    x_coords = load_1d_csv(x_coords_file)
    y_coords = load_1d_csv(y_coords_file)
    geo_params = load_geometry_params_for_plot(geometry_params_file)

    if V is None or Ex is None or Ey is None or eps_r is None or x_coords is None or y_coords is None:
        print("One or more data files could not be loaded. Aborting plot.")
        return None # Return None on failure
    
    # Attempt to load geometry parameters for profile plot, proceed if available
    y_center_gap_idx = None
    if geo_params:
        y_si_base_h = geo_params.get("y_si_base_height")
        # Try 'initial_y_teeth_height' first, then 'y_teeth_height' as fallback
        y_teeth_h_val = geo_params.get("initial_y_teeth_height", geo_params.get("y_teeth_height"))
        y_vac_gap_thick = geo_params.get("y_vacuum_gap_thick")

        if y_si_base_h is not None and y_teeth_h_val is not None and y_vac_gap_thick is not None:
            y_center_gap_abs = y_si_base_h + y_teeth_h_val + (y_vac_gap_thick / 2.0)
            # Find the closest y-index
            y_center_gap_idx = (np.abs(y_coords - y_center_gap_abs)).argmin()
            print(f"Calculated y-center for profile plot: {y_center_gap_abs:.2f} µm (index: {y_center_gap_idx})")
        else:
            print("Warning: Could not determine vacuum gap center from geometry parameters. Profile plot will be skipped.")
    else:
        print("Warning: Geometry parameters not loaded. Profile plot will be skipped.")


    E_mag = np.sqrt(Ex**2 + Ey**2)

    # Geometry parameters needed for quiver plot mask (No longer needed if using eps_r for mask)
    # y_si_base_height = geo_params["y_si_base_height"]
    # y_teeth_height = geo_params["y_teeth_height"]
    # y_vacuum_gap_thick = geo_params["y_vacuum_gap_thick"]

    # Y-coordinates for vacuum gap boundaries (used in quiver plot mask)
    # y_bot_si_teeth_top = y_si_base_height + y_teeth_height # No longer needed
    # y_top_si_teeth_bottom = y_bot_si_teeth_top + y_vacuum_gap_thick # No longer needed
    
    def draw_detailed_outlines(ax, x_coords_mesh, y_coords_mesh, eps_r_data_mesh, threshold, color_style):
        """Draws geometry outlines based on the permittivity map."""
        # Separate color and linestyle
        # Assuming color_style is like 'w--' or 'k--'
        color = color_style[0]
        linestyle = color_style[1:] if len(color_style) > 1 else '-' # Default to solid if no style part

        ax.contour(x_coords_mesh, y_coords_mesh, eps_r_data_mesh, levels=[threshold], colors=color, linestyles=linestyle, linewidths=0.8)


    # --- Plotting on separate canvases ---

    # Prepare meshgrid for contour plots (potential, E-field, permittivity, and outlines)
    # Note: x_coords, y_coords are 1D. eps_r, V, Ex, Ey are loaded as (Ny, Nx) due to .T
    X_mesh, Y_mesh = np.meshgrid(x_coords, y_coords)
    
    # Define a threshold for distinguishing silicon from vacuum
    # Assuming eps_vac = 1.0 and eps_si = 11.7 (typical values from poisson_solver.cpp)
    eps_vac_assumed = 1.0
    eps_si_assumed = 11.7 
    # It's better if eps_si is read from params or passed, but for now, this is a common value.
    # If 'eps_si' was saved in geometry_params.csv, we could use it:
    # eps_si_from_params = geo_params.get('eps_si', 11.7) # Example
    outline_threshold = (eps_vac_assumed + eps_si_assumed) / 2.0


    # Plot 1: Electric Potential
    plt.figure(figsize=(8, 6)) # New figure for Potential
    ax_V = plt.gca()
    contour_V = ax_V.contourf(X_mesh, Y_mesh, V, levels=50, cmap='viridis') # V is already (Ny, Nx)
    plt.colorbar(contour_V, ax=ax_V, label='Potential (V)')
    ax_V.set_title('Electric Potential (V)')
    ax_V.set_xlabel('x (µm)')
    ax_V.set_ylabel('y (µm)')
    ax_V.set_aspect('equal', adjustable='box')
    draw_detailed_outlines(ax_V, X_mesh, Y_mesh, eps_r, outline_threshold, 'w--') # eps_r is (Ny, Nx)
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder_name, "potential_plot.png"), dpi=300)

    # Plot 2: Electric Field Magnitude
    plt.figure(figsize=(8, 6)) # New figure for E-field Magnitude
    ax_Emag = plt.gca()
    contour_Emag = ax_Emag.contourf(X_mesh, Y_mesh, E_mag, levels=50, cmap='inferno') # E_mag is (Ny, Nx)
    plt.colorbar(contour_Emag, ax=ax_Emag, label='Electric Field Magnitude (V/µm)')
    ax_Emag.set_title('Electric Field Magnitude |E|')
    ax_Emag.set_xlabel('x (µm)')
    ax_Emag.set_ylabel('y (µm)')
    ax_Emag.set_aspect('equal', adjustable='box')
    draw_detailed_outlines(ax_Emag, X_mesh, Y_mesh, eps_r, outline_threshold, 'w--')
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder_name, "efield_magnitude_plot.png"), dpi=300)

    # Plot 3: Permittivity Map
    plt.figure(figsize=(8, 6)) # New figure for Permittivity
    ax_eps = plt.gca()
    contour_eps = ax_eps.contourf(X_mesh, Y_mesh, eps_r, levels=np.linspace(np.min(eps_r), np.max(eps_r), 5), cmap='coolwarm') # eps_r is (Ny, Nx)
    plt.colorbar(contour_eps, ax=ax_eps, label='Relative Permittivity (ε$_r$)')
    ax_eps.set_title('Relative Permittivity Map')
    ax_eps.set_xlabel('x (µm)')
    ax_eps.set_ylabel('y (µm)')
    ax_eps.set_aspect('equal', adjustable='box')
    # Optionally, draw outlines on the permittivity map itself, perhaps with a different color
    draw_detailed_outlines(ax_eps, X_mesh, Y_mesh, eps_r, outline_threshold, 'k--') 
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder_name, "permittivity_map_plot.png"), dpi=300)

    # Plot 4: Electric Field Quiver Plot (Vacuum Only)
    plt.figure(figsize=(10, 7)) # New figure for E-field Quiver
    ax_Evec = plt.gca()
    
    # Downsample for quiver plot to avoid overcrowding
    skip_rate = 5 # Adjust skip rate as needed
    # Ensure skip rate is not too large for the number of points, especially in y
    if y_coords.shape[0] // skip_rate < 2 : # Ensure at least 2 points in y after skipping if possible
        skip_rate_y = max(1, y_coords.shape[0] // 2 if y_coords.shape[0] > 1 else 1)
    else:
        skip_rate_y = skip_rate
    if x_coords.shape[0] // skip_rate < 2 :
        skip_rate_x = max(1, x_coords.shape[0] // 2 if x_coords.shape[0] > 1 else 1)
    else:
        skip_rate_x = skip_rate

    skip = (slice(None, None, skip_rate_y), slice(None, None, skip_rate_x))
    
    X, Y = np.meshgrid(x_coords, y_coords) # Meshgrid for quiver

    # Prepare data for quiver: Ex, Ey, E_mag are already (Ny, Nx)
    Ex_quiver = Ex.copy()
    Ey_quiver = Ey.copy()
    Emag_quiver = E_mag.copy()

    # Create a mask for points outside the vacuum gap using the permittivity map
    # Vacuum is where eps_r is close to eps_vac_assumed (e.g., < outline_threshold)
    # Mask should be True for non-vacuum regions.
    # eps_r is already (Ny, Nx)
    mask_not_vacuum = (eps_r >= outline_threshold) 

    # Apply mask: set values outside vacuum to NaN so they are not plotted
    Ex_quiver[mask_not_vacuum] = np.nan
    Ey_quiver[mask_not_vacuum] = np.nan
    Emag_quiver[mask_not_vacuum] = np.nan
    
    # Plot quiver using the masked and skipped data
    # Only plot if there are non-NaN values to avoid errors with all-NaN slices
    if not np.all(np.isnan(Ex_quiver[skip])):
        ax_Evec.quiver(X_mesh[skip], Y_mesh[skip], 
                       Ex_quiver[skip], Ey_quiver[skip], Emag_quiver[skip], 
                       cmap='viridis', scale=None, scale_units='xy', angles='xy', 
                       headwidth=3, headlength=5, pivot='middle')
    
    ax_Evec.set_title('Electric Field Vectors in Vacuum (Quiver Plot)')
    ax_Evec.set_xlabel('x (µm)')
    ax_Evec.set_ylabel('y (µm)')
    ax_Evec.set_aspect('equal', adjustable='box')
    # Add structure outlines
    draw_detailed_outlines(ax_Evec, X_mesh, Y_mesh, eps_r, outline_threshold, 'k--')
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder_name, "efield_quiver_vacuum_plot.png"), dpi=300) # Renamed save file

    # Plot 5: Profile plot at the center of the vacuum gap
    if y_center_gap_idx is not None:
        V_profile = V[y_center_gap_idx, :]
        Emag_profile = E_mag[y_center_gap_idx, :]

        fig_profile, ax_profile_V = plt.subplots(figsize=(10, 6))

        color_V = 'tab:blue'
        ax_profile_V.set_xlabel('x (µm)')
        ax_profile_V.set_ylabel('Potential (V)', color=color_V)
        ax_profile_V.plot(x_coords, V_profile, color=color_V, linestyle='-', label='Potential (V)')
        ax_profile_V.tick_params(axis='y', labelcolor=color_V)
        ax_profile_V.grid(True, linestyle=':', alpha=0.7)

        ax_profile_Emag = ax_profile_V.twinx()  # instantiate a second axes that shares the same x-axis
        color_Emag = 'tab:red'
        ax_profile_Emag.set_ylabel('Electric Field Magnitude (V/µm)', color=color_Emag)
        ax_profile_Emag.plot(x_coords, Emag_profile, color=color_Emag, linestyle='--', label='|E| (V/µm)')
        ax_profile_Emag.tick_params(axis='y', labelcolor=color_Emag)

        fig_profile.suptitle(f'Profile at y = {y_coords[y_center_gap_idx]:.2f} µm (Center of Vacuum Gap)')
        # To add a combined legend:
        lines_V, labels_V = ax_profile_V.get_legend_handles_labels()
        lines_Emag, labels_Emag = ax_profile_Emag.get_legend_handles_labels()
        ax_profile_Emag.legend(lines_V + lines_Emag, labels_V + labels_Emag, loc='upper right')
        
        fig_profile.tight_layout(rect=[0, 0, 1, 0.96]) # Adjust layout to make space for suptitle
        plt.savefig(os.path.join(output_folder_name, "center_gap_profile_plot.png"), dpi=300)
        print(f"Profile plot saved to {os.path.join(output_folder_name, 'center_gap_profile_plot.png')}")

    # --- Proton Trajectory Analysis ---
    # Initialize h_param_um, it might be set from geo_params if available
    h_param_um = None
    if geo_params and geo_params.get("h") is not None:
        h_param_um = geo_params.get("h") # Assuming 'h' is in µm in geometry_params.csv for this script's context
    elif x_coords is not None and len(x_coords) > 1:
        # Fallback: try to derive from x_coordinates if not in geo_params
        h_param_um = np.abs(x_coords[1] - x_coords[0])
        print(f"Warning: 'h' not found in geometry_params.csv. Using h_param_um derived from x_coordinates: {h_param_um:.2f} µm")
    else:
        print("Warning: 'h' could not be determined from geometry_params.csv or x_coordinates. Proton analysis might be affected or skipped.")

    print("Proton trajectory analysis and histogram plotting has been removed from this script.")

    plt.show() # Show all figures
    return output_folder_name # Return for use in __main__

if __name__ == '__main__':
    cli_folder_path = None
    if len(sys.argv) > 1:
        cli_folder_path = sys.argv[1]
    
    output_dir = plot_results(folder_path=cli_folder_path) # Pass the folder path

    if output_dir: # Check if plot_results returned a valid path
        print(f"Plotting finished. Plots saved to {output_dir} folder.")
    else:
        print("Plotting did not complete successfully.")
