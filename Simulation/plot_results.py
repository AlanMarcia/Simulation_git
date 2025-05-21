import numpy as np
import matplotlib.pyplot as plt
import os
import csv

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

def plot_results():
    # Create output folder
    output_folder = "geometria_Denti_sfasati"
    os.makedirs(output_folder, exist_ok=True)

    # File paths
    potential_file = os.path.join(output_folder, "potential.csv")
    ex_file = os.path.join(output_folder, "electric_field_x.csv")
    ey_file = os.path.join(output_folder, "electric_field_y.csv")
    eps_r_file = os.path.join(output_folder, "permittivity.csv")
    x_coords_file = os.path.join(output_folder, "x_coordinates.csv")
    y_coords_file = os.path.join(output_folder, "y_coordinates.csv")
    # geometry_params_file = os.path.join(output_folder, "geometry_params.csv") # No longer explicitly used for plotting

    # Load data
    V = load_2d_csv(potential_file)
    Ex = load_2d_csv(ex_file)
    Ey = load_2d_csv(ey_file)
    eps_r = load_2d_csv(eps_r_file)
    x_coords = load_1d_csv(x_coords_file)
    y_coords = load_1d_csv(y_coords_file)

    if V is None or Ex is None or Ey is None or eps_r is None or x_coords is None or y_coords is None:
        print("One or more data files could not be loaded. Aborting plot.")
        return

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
    plt.savefig(os.path.join(output_folder, "potential_plot.pdf"))

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
    plt.savefig(os.path.join(output_folder, "efield_magnitude_plot.pdf"))

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
    plt.savefig(os.path.join(output_folder, "permittivity_map_plot.pdf"))

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
    plt.savefig(os.path.join(output_folder, "efield_quiver_vacuum_plot.pdf")) # Renamed save file


    plt.show() # Show all figures

if __name__ == '__main__':
    plot_results()
    print(f"Plotting finished. Plots saved to {output_folder} folder.")
