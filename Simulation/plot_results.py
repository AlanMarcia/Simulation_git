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
    return data.T

def load_geometry_params(filename):
    """Loads geometry parameters from a CSV file."""
    if not os.path.exists(filename):
        print(f"Error: File {filename} not found.")
        return None
    params = {}
    with open(filename, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row) == 2:
                params[row[0]] = float(row[1])
    required_keys = ["h", "x_free_space", "x_structure_len", "y_si_layer_thick", "y_vacuum_gap_thick", "H_total"]
    if not all(key in params for key in required_keys):
        print(f"Error: Missing one or more required keys in {filename}. Required: {required_keys}")
        return None
    return params

def plot_results():
    # File paths
    potential_file = "potential.csv"
    ex_file = "electric_field_x.csv"
    ey_file = "electric_field_y.csv"
    eps_r_file = "permittivity.csv"
    x_coords_file = "x_coordinates.csv"
    y_coords_file = "y_coordinates.csv"
    geometry_params_file = "geometry_params.csv"

    # Load geometry parameters
    geo_params = load_geometry_params(geometry_params_file)
    if geo_params is None:
        print("Could not load geometry parameters. Aborting plot.")
        return

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

    # Geometry for outlines from loaded parameters
    # h_sim = geo_params["h"] # h from simulation, might be useful for other things
    x_free_space = geo_params["x_free_space"]
    x_structure_len = geo_params["x_structure_len"]
    y_si_layer_thick = geo_params["y_si_layer_thick"]
    y_vacuum_gap_thick = geo_params["y_vacuum_gap_thick"]
    H_total_sim = geo_params["H_total"] # Total height from simulation

    x_struct_start = x_free_space
    x_struct_end = x_free_space + x_structure_len
    y_si_bot_end = y_si_layer_thick
    y_vac_end = y_si_layer_thick + y_vacuum_gap_thick # This is the top interface of the vacuum gap
    y_si_top_start = y_si_layer_thick + y_vacuum_gap_thick # This is the bottom interface of the top Si layer

    # --- Plotting on separate canvases ---

    # Plot 1: Electric Potential
    plt.figure(figsize=(8, 6)) # New figure for Potential
    ax_V = plt.gca()
    contour_V = ax_V.contourf(x_coords, y_coords, V.T, levels=50, cmap='viridis') # V.T because contourf expects Z(Y,X)
    plt.colorbar(contour_V, ax=ax_V, label='Potential (V)')
    ax_V.set_title('Electric Potential (V)')
    ax_V.set_xlabel('x (µm)')
    ax_V.set_ylabel('y (µm)')
    ax_V.set_aspect('equal', adjustable='box')
    # Add structure outlines
    ax_V.plot([x_struct_start, x_struct_end], [y_si_bot_end, y_si_bot_end], 'w--', lw=0.8)
    ax_V.plot([x_struct_start, x_struct_end], [y_vac_end, y_vac_end], 'w--', lw=0.8)
    ax_V.plot([x_struct_start, x_struct_start], [0, y_si_bot_end], 'w--', lw=0.8)
    ax_V.plot([x_struct_start, x_struct_start], [y_si_top_start, H_total_sim], 'w--', lw=0.8)
    ax_V.plot([x_struct_end, x_struct_end], [0, y_si_bot_end], 'w--', lw=0.8)
    ax_V.plot([x_struct_end, x_struct_end], [y_si_top_start, H_total_sim], 'w--', lw=0.8)
    plt.tight_layout()
    plt.savefig("potential_plot.png")

    # Plot 2: Electric Field Magnitude
    plt.figure(figsize=(8, 6)) # New figure for E-field Magnitude
    ax_Emag = plt.gca()
    contour_Emag = ax_Emag.contourf(x_coords, y_coords, E_mag.T, levels=50, cmap='inferno') # E_mag.T
    plt.colorbar(contour_Emag, ax=ax_Emag, label='Electric Field Magnitude (V/µm)')
    ax_Emag.set_title('Electric Field Magnitude |E|')
    ax_Emag.set_xlabel('x (µm)')
    ax_Emag.set_ylabel('y (µm)')
    ax_Emag.set_aspect('equal', adjustable='box')
    # Add structure outlines
    ax_Emag.plot([x_struct_start, x_struct_end], [y_si_bot_end, y_si_bot_end], 'w--', lw=0.8)
    ax_Emag.plot([x_struct_start, x_struct_end], [y_vac_end, y_vac_end], 'w--', lw=0.8)
    ax_Emag.plot([x_struct_start, x_struct_start], [0, y_si_bot_end], 'w--', lw=0.8)
    ax_Emag.plot([x_struct_start, x_struct_start], [y_si_top_start, H_total_sim], 'w--', lw=0.8)
    ax_Emag.plot([x_struct_end, x_struct_end], [0, y_si_bot_end], 'w--', lw=0.8)
    ax_Emag.plot([x_struct_end, x_struct_end], [y_si_top_start, H_total_sim], 'w--', lw=0.8)
    plt.tight_layout()
    plt.savefig("efield_magnitude_plot.png")

    # Plot 3: Permittivity Map
    plt.figure(figsize=(8, 6)) # New figure for Permittivity
    ax_eps = plt.gca()
    contour_eps = ax_eps.contourf(x_coords, y_coords, eps_r.T, levels=10, cmap='coolwarm') # eps_r.T
    plt.colorbar(contour_eps, ax=ax_eps, label='Relative Permittivity (ε$_r$)')
    ax_eps.set_title('Relative Permittivity Map')
    ax_eps.set_xlabel('x (µm)')
    ax_eps.set_ylabel('y (µm)')
    ax_eps.set_aspect('equal', adjustable='box')
    # Add structure outlines
    ax_eps.plot([x_struct_start, x_struct_end], [y_si_bot_end, y_si_bot_end], 'w--', lw=0.8)
    ax_eps.plot([x_struct_start, x_struct_end], [y_vac_end, y_vac_end], 'w--', lw=0.8)
    ax_eps.plot([x_struct_start, x_struct_start], [0, y_si_bot_end], 'w--', lw=0.8)
    ax_eps.plot([x_struct_start, x_struct_start], [y_si_top_start, H_total_sim], 'w--', lw=0.8)
    ax_eps.plot([x_struct_end, x_struct_end], [0, y_si_bot_end], 'w--', lw=0.8)
    ax_eps.plot([x_struct_end, x_struct_end], [y_si_top_start, H_total_sim], 'w--', lw=0.8)
    plt.tight_layout()
    plt.savefig("permittivity_map_plot.png")

    plt.show() # Show all figures

if __name__ == '__main__':
    plot_results()
    print("Plotting finished. Plots saved to potential_plot.png, efield_magnitude_plot.png, and permittivity_map_plot.png")
