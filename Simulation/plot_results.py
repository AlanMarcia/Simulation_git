import numpy as np
import matplotlib.pyplot as plt
import os
import csv
import glob # For finding trajectory files
import sys # Import sys module for command-line arguments
import matplotlib.colors as mcolors

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
    y_center_gap_abs = None # Initialize
    if geo_params:
        # Try logic for toothed structures first
        y_si_base_h = geo_params.get("y_si_base_height")
        y_teeth_h_val = geo_params.get("initial_y_teeth_height") or geo_params.get("y_teeth_height") or geo_params.get("y_tooth_height") # Try initial, then plural, then singular
        y_vac_gap_thick_toothed = geo_params.get("y_vacuum_gap_thick")
        pad_bottom = geo_params.get("y_vacuum_padding_bottom", 0.0)

        if y_si_base_h is not None and y_teeth_h_val is not None and y_vac_gap_thick_toothed is not None:
            y_center_gap_abs = pad_bottom + y_si_base_h + y_teeth_h_val + (y_vac_gap_thick_toothed / 2.0)
            print(f"Using toothed geometry logic for profile plot y-center.")
        else:
            # Try piana_rastremata or piana_variabile geometry
            print("Toothed geometry parameters not fully found. Trying 'piana_rastremata' or 'piana_variabile' logic for profile plot y-center.")
            y_si_thick_start = geo_params.get("y_si_layer_thick_start")
            y_si_thick_end = geo_params.get("y_si_layer_thick_end")
            y_si_thick_right = geo_params.get("y_si_layer_thick_right")
            y_si_thick_left = geo_params.get("y_si_layer_thick_left")
            y_vacuum_gap_thick_var = geo_params.get("y_vacuum_gap_thick")
            
            # For piana_rastremata or piana_variabile: use thickness at the start (x=0)
            if y_si_thick_start is not None and y_si_thick_end is not None and y_vacuum_gap_thick_var is not None:
                # For piana_rastremata: thickness varies linearly, use thickness at start (x=0)
                y_center_gap_abs = pad_bottom + y_si_thick_start + (y_vacuum_gap_thick_var / 2.0)
                print(f"Using 'piana_rastremata' logic for profile plot y-center (thickness at x=0: {y_si_thick_start:.2f} μm).")
            elif y_si_thick_right is not None and y_si_thick_left is not None and y_vacuum_gap_thick_var is not None:
                # For piana_variabile: thickness varies linearly, use thickness at start (x=0, left side)
                y_center_gap_abs = pad_bottom + y_si_thick_left + (y_vacuum_gap_thick_var / 2.0)
                print(f"Using 'piana_variabile' logic for profile plot y-center (thickness at x=0: {y_si_thick_left:.2f} μm).")
            else:
                # Fallback logic for simple "geometria_piana"
                print("Rastremata/variabile geometry parameters not found. Trying simple 'geometria_piana' logic for profile plot y-center.")
                y_si_layer_thick_piana = y_si_thick_start
                y_vacuum_gap_thick_piana = geo_params.get("y_vacuum_gap_thick")

                if y_si_layer_thick_piana is not None and y_vacuum_gap_thick_piana is not None:
                    y_center_gap_abs = pad_bottom + y_si_layer_thick_piana + (y_vacuum_gap_thick_piana / 2.0)
                    print(f"Using 'geometria_piana' logic for profile plot y-center.")
                else:
                    print("Warning: Could not determine vacuum gap center from available geometry parameters. Profile plot will be skipped.")
        
        if y_center_gap_abs is not None:
            # Find the closest y-index
            if y_coords is not None and len(y_coords) > 0:
                y_center_gap_idx = (np.abs(y_coords - y_center_gap_abs)).argmin()
                print(f"Calculated y-center for profile plot: {y_center_gap_abs:.2f} ┬╡m (index: {y_center_gap_idx})")
            else:
                print("Warning: y_coords not available for determining profile plot index.")
                y_center_gap_idx = None # Ensure it's None if y_coords are missing
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
        if threshold is None: # If no valid threshold, do not attempt to draw outlines
            # print("Skipping outlines as threshold is None.")
            return

        # Separate color and linestyle
        # Assuming color_style is like 'w--' or 'k--'
        color = color_style[0]
        linestyle = color_style[1:] if len(color_style) > 1 else '-' # Default to solid if no style part

        ax.contour(x_coords_mesh, y_coords_mesh, eps_r_data_mesh, levels=[threshold], colors=color, linestyles=linestyle, linewidths=0.8)


    # --- Plotting on separate canvases ---

    # Prepare meshgrid for contour plots (potential, E-field, permittivity, and outlines)
    # Note: x_coords, y_coords are 1D. eps_r, V, Ex, Ey are loaded as (Ny, Nx) due to .T
    X_mesh, Y_mesh = np.meshgrid(x_coords, y_coords)
    
    # Define thresholds for distinguishing aluminum, silicon, and vacuum
    # Dynamically calculate based on loaded eps_r data
    outline_threshold_silicon = None
    outline_threshold_aluminum = None
    
    if eps_r is not None:
        unique_eps_values = np.unique(eps_r)
        unique_eps_values.sort() # Ensure sorted for min/max logic
        
        if len(unique_eps_values) >= 3:
            # We have vacuum (lowest), silicon (middle), aluminum (highest)
            val_vacuum = unique_eps_values[0]
            val_silicon = unique_eps_values[1]
            val_aluminum = unique_eps_values[-1]
            
            # Threshold between vacuum and silicon
            outline_threshold_silicon = (val_vacuum + val_silicon) / 2.0
            # Threshold between silicon and aluminum
            outline_threshold_aluminum = (val_silicon + val_aluminum) / 2.0
            
            print(f"Detected 3 materials:")
            print(f"  Vacuum: ε_r = {val_vacuum:.2f}")
            print(f"  Silicon: ε_r = {val_silicon:.2f}")
            print(f"  Aluminum: ε_r = {val_aluminum:.2f}")
            print(f"  Silicon outline threshold: {outline_threshold_silicon:.2f}")
            print(f"  Aluminum outline threshold: {outline_threshold_aluminum:.2f}")
            
        elif len(unique_eps_values) >= 2:
            # Assume the two most extreme values (after sorting) represent vacuum and material
            val_low = unique_eps_values[0]
            val_high = unique_eps_values[-1]
            # Check if these values are reasonably distinct to represent two phases
            if val_high > val_low + 0.5: # Heuristic: difference must be at least 0.5
                outline_threshold_silicon = (val_low + val_high) / 2.0
                print(f"Dynamically calculated outline threshold: {outline_threshold_silicon:.2f} (from eps_r min/max: {val_low:.2f}, {val_high:.2f})")
            else:
                print(f"Warning: Unique permittivity values ({unique_eps_values}) are too close. Could not reliably determine outline threshold.")
        elif len(unique_eps_values) == 1:
            print(f"Warning: Permittivity data contains only one unique value ({unique_eps_values[0]:.2f}). No outlines will be drawn.")
        else: # Should not happen if eps_r is loaded and not empty
            print("Warning: Could not process permittivity data for outlines.")
    else:
        print("Warning: Permittivity data (eps_r) is None. Cannot draw outlines.")

    # eps_vac_assumed = 1.0 # No longer primary method for threshold
    # eps_si_assumed = 11.7 
    # outline_threshold = (eps_vac_assumed + eps_si_assumed) / 2.0 # Old method


    # Plot 1: Electric Potential
    plt.figure(figsize=(8, 6)) # New figure for Potential
    ax_V = plt.gca()
    contour_V = ax_V.contourf(X_mesh, Y_mesh, V, levels=50, cmap='viridis') # V is already (Ny, Nx)
    plt.colorbar(contour_V, ax=ax_V, label='Potential (V)')
    ax_V.set_title('Electric Potential (V)')
    ax_V.set_xlabel('x (μm)')
    ax_V.set_ylabel('y (μm)')
    ax_V.set_aspect('equal', adjustable='box')
    draw_detailed_outlines(ax_V, X_mesh, Y_mesh, eps_r, outline_threshold_silicon, 'w--') # Silicon outline
    draw_detailed_outlines(ax_V, X_mesh, Y_mesh, eps_r, outline_threshold_aluminum, 'r-') # Aluminum outline
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder_name, "potential_plot.png"), dpi=300)

    # Plot 2: Electric Field Magnitude (Vacuum Only)
    plt.figure(figsize=(8,3)) # New figure for E-field Magnitude
    ax_Emag = plt.gca()
    
    # Create a masked version of E_mag for vacuum only
    E_mag_vacuum = E_mag.copy()
    if outline_threshold_silicon is not None:
        # Mask regions where eps_r >= threshold (i.e., silicon and aluminum)
        mask_not_vacuum = (eps_r >= outline_threshold_silicon)
        E_mag_vacuum[mask_not_vacuum] = np.nan  # Set non-vacuum regions to NaN
    else:
        print("Warning: Could not determine vacuum regions for E_mag plot. Plotting everywhere.")
    
    contour_Emag = ax_Emag.contourf(X_mesh, Y_mesh, E_mag_vacuum, levels=30, cmap='gnuplot') # E_mag_vacuum is (Ny, Nx), only vacuum
    plt.colorbar(contour_Emag, ax=ax_Emag, label='Electric Field Magnitude (V/μm)', orientation='horizontal')
    ax_Emag.set_title('Electric Field Magnitude |E| (Vacuum Only)')
    ax_Emag.set_xlabel('x (μm)')
    ax_Emag.set_ylabel('y (μm)')
    ax_Emag.set_aspect('equal', adjustable='box')
    draw_detailed_outlines(ax_Emag, X_mesh, Y_mesh, eps_r, outline_threshold_silicon, 'k--')
    draw_detailed_outlines(ax_Emag, X_mesh, Y_mesh, eps_r, outline_threshold_aluminum, 'r-')
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder_name, "efield_magnitude_plot.png"), dpi=300)

    # Plot 3: Permittivity Map
    plt.figure(figsize=(8, 6)) # New figure for Permittivity
    ax_eps = plt.gca()
    
    # Check if permittivity has variation
    eps_min = np.min(eps_r)
    eps_max = np.max(eps_r)
    
    if eps_max > eps_min + 1e-6:  # Has variation
        contour_eps = ax_eps.contourf(X_mesh, Y_mesh, eps_r, levels=np.linspace(eps_min, eps_max, 5), cmap='RdBu')
        plt.colorbar(contour_eps, ax=ax_eps, label='Relative Permittivity (εᵣ)')
    else:  # Constant value - use imshow instead
        im_eps = ax_eps.imshow(eps_r, extent=[X_mesh.min(), X_mesh.max(), Y_mesh.min(), Y_mesh.max()], 
                               origin='lower', cmap='RdBu', aspect='auto')
        plt.colorbar(im_eps, ax=ax_eps, label='Relative Permittivity (εᵣ)')
        ax_eps.text(0.5, 0.95, f'Uniform εᵣ = {eps_min:.2f}', 
                   transform=ax_eps.transAxes, ha='center', va='top',
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    ax_eps.set_title('Relative Permittivity Map')
    ax_eps.set_xlabel('x (μm)')
    ax_eps.set_ylabel('y (μm)')
    ax_eps.set_aspect('equal', adjustable='box')
    # Optionally, draw outlines on the permittivity map itself, perhaps with a different color
    draw_detailed_outlines(ax_eps, X_mesh, Y_mesh, eps_r, outline_threshold_silicon, 'k--')
    draw_detailed_outlines(ax_eps, X_mesh, Y_mesh, eps_r, outline_threshold_aluminum, 'r-')
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder_name, "permittivity_map_plot.png"), dpi=300)

    # Plot 4: Electric Field Quiver Plot (Vacuum Only)
    plt.figure(figsize=(10, 7)) # New figure for E-field Quiver
    ax_Evec = plt.gca()
    
    # Downsample for quiver plot to avoid overcrowding
    skip_rate = 3 # Adjust skip rate as needed
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
    # Vacuum is where eps_r is close to eps_vac_assumed (e.g., < outline_threshold_silicon)
    # Mask should be True for non-vacuum regions.
    # eps_r is already (Ny, Nx)
    if outline_threshold_silicon is not None:
        mask_not_vacuum = (eps_r >= outline_threshold_silicon) 
    else: # If no threshold, assume all is vacuum for quiver (or handle differently)
        mask_not_vacuum = np.zeros_like(eps_r, dtype=bool) # No mask, plot everywhere
        print("Warning: Quiver plot mask could not be determined due to missing outline_threshold. Plotting vectors everywhere.")


    # Apply mask: set values outside vacuum to NaN so they are not plotted
    Ex_quiver[mask_not_vacuum] = np.nan
    Ey_quiver[mask_not_vacuum] = np.nan
    Emag_quiver[mask_not_vacuum] = np.nan
    
    # Plot quiver using the masked and skipped data
    # Only plot if there are non-NaN values to avoid errors with all-NaN slices
    if not np.all(np.isnan(Ex_quiver[skip])):
        ax_Evec.quiver(X_mesh[skip], Y_mesh[skip], 
                       Ex_quiver[skip], Ey_quiver[skip], Emag_quiver[skip], 
                       cmap='viridis', scale=20, scale_units='xy', angles='xy', 
                       headwidth=3, headlength=5, pivot='middle')
    
    ax_Evec.set_title('Electric Field Vectors in Vacuum (Quiver Plot)')
    ax_Evec.set_xlabel('x (μm)')
    ax_Evec.set_ylabel('y (μm)')
    ax_Evec.set_aspect('equal', adjustable='box')
    ax_Evec.set_xlim(15,60)
    ax_Evec.set_ylim(70,150)
    # Add structure outlines
    draw_detailed_outlines(ax_Evec, X_mesh, Y_mesh, eps_r, outline_threshold_silicon, 'k--')
    draw_detailed_outlines(ax_Evec, X_mesh, Y_mesh, eps_r, outline_threshold_aluminum, 'r-')
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder_name, "efield_quiver_vacuum_plot.png"), dpi=300) # Renamed save file

    # Plot 5: Profile plot at the center of the vacuum gap
    if y_center_gap_idx is not None:
        V_profile = V[y_center_gap_idx, :]
        Ex_profile = Ex[y_center_gap_idx, :]
        Ey_profile = Ey[y_center_gap_idx, :]
        Emag_profile = E_mag[y_center_gap_idx, :]

        # --- Plot 5a: Combined Potential and E-field magnitude profile ---
        fig_profile, ax_profile_V = plt.subplots(figsize=(12, 6))

        color_V = 'tab:blue'
        ax_profile_V.set_xlabel('x (μm)', fontsize=12)
        ax_profile_V.set_ylabel('Potential (V)', color=color_V, fontsize=12)
        ax_profile_V.plot(x_coords, V_profile, color=color_V, linestyle='-', linewidth=2, label='Potential (V)')
        ax_profile_V.tick_params(axis='y', labelcolor=color_V)
        ax_profile_V.grid(True, linestyle=':', alpha=0.7)

        ax_profile_Emag = ax_profile_V.twinx()  # instantiate a second axes that shares the same x-axis
        color_Emag = 'tab:red'
        ax_profile_Emag.set_ylabel('Electric Field Magnitude (V/μm)', color=color_Emag, fontsize=12)
        ax_profile_Emag.plot(x_coords, Emag_profile, color=color_Emag, linestyle='--', linewidth=2, label='|E| (V/μm)')
        ax_profile_Emag.tick_params(axis='y', labelcolor=color_Emag)

        fig_profile.suptitle(f'Potential and E-field Profile at y = {y_coords[y_center_gap_idx]:.2f} μm (Center of Acceleration Channel)', fontsize=14, fontweight='bold')
        # To add a combined legend:
        lines_V, labels_V = ax_profile_V.get_legend_handles_labels()
        lines_Emag, labels_Emag = ax_profile_Emag.get_legend_handles_labels()
        ax_profile_Emag.legend(lines_V + lines_Emag, labels_V + labels_Emag, loc='upper right', fontsize=10)
        
        fig_profile.tight_layout(rect=[0, 0, 1, 0.96]) # Adjust layout to make space for suptitle
        plt.savefig(os.path.join(output_folder_name, "center_gap_profile_plot.png"), dpi=300)
        print(f"Combined profile plot saved to {os.path.join(output_folder_name, 'center_gap_profile_plot.png')}")
        
        # --- Plot 5b: Dedicated E-field components profile ---
        fig_efield, (ax_ex, ax_ey, ax_emag) = plt.subplots(3, 1, figsize=(12, 10), sharex=True)
        
        # Ex component
        ax_ex.plot(x_coords, Ex_profile, color='darkblue', linestyle='-', linewidth=2)
        ax_ex.set_ylabel('Ex (V/μm)', fontsize=11, fontweight='bold')
        ax_ex.set_title('Electric Field Components along Acceleration Channel Center', fontsize=13, fontweight='bold')
        ax_ex.grid(True, linestyle=':', alpha=0.7)
        ax_ex.axhline(y=0, color='black', linestyle='--', linewidth=0.8, alpha=0.5)
        
        # Ey component
        ax_ey.plot(x_coords, Ey_profile, color='darkgreen', linestyle='-', linewidth=2)
        ax_ey.set_ylabel('Ey (V/μm)', fontsize=11, fontweight='bold')
        ax_ey.grid(True, linestyle=':', alpha=0.7)
        ax_ey.axhline(y=0, color='black', linestyle='--', linewidth=0.8, alpha=0.5)
        
        # E magnitude
        ax_emag.plot(x_coords, Emag_profile, color='darkred', linestyle='-', linewidth=2)
        ax_emag.set_ylabel('|E| (V/μm)', fontsize=11, fontweight='bold')
        ax_emag.set_xlabel('x (μm)', fontsize=12)
        ax_emag.grid(True, linestyle=':', alpha=0.7)
        
        # Add info text box with statistics
        ex_max = np.max(np.abs(Ex_profile))
        ey_max = np.max(np.abs(Ey_profile))
        emag_max = np.max(Emag_profile)
        emag_mean = np.mean(Emag_profile)
        
        stats_text = f'Max |Ex|: {ex_max:.3f} V/μm\nMax |Ey|: {ey_max:.3f} V/μm\nMax |E|: {emag_max:.3f} V/μm\nMean |E|: {emag_mean:.3f} V/μm'
        ax_emag.text(0.98, 0.97, stats_text, transform=ax_emag.transAxes, fontsize=9,
                    verticalalignment='top', horizontalalignment='right',
                    bbox=dict(boxstyle='round,pad=0.5', fc='lightyellow', alpha=0.8))
        
        fig_efield.suptitle(f'at y = {y_coords[y_center_gap_idx]:.2f} μm', fontsize=11, style='italic')
        fig_efield.tight_layout(rect=[0, 0, 1, 0.98])
        plt.savefig(os.path.join(output_folder_name, "efield_components_profile_plot.png"), dpi=300)
        print(f"E-field components profile plot saved to {os.path.join(output_folder_name, 'efield_components_profile_plot.png')}")
        
        # --- Plot 5c: E-field magnitude with highlighted regions ---
        fig_emag_detail, ax_emag_detail = plt.subplots(figsize=(14, 6))
        
        ax_emag_detail.plot(x_coords, Emag_profile, color='crimson', linestyle='-', linewidth=2.5, label='|E| magnitude')
        ax_emag_detail.fill_between(x_coords, 0, Emag_profile, alpha=0.3, color='coral')
        
        # Add horizontal line for mean
        ax_emag_detail.axhline(y=emag_mean, color='navy', linestyle='--', linewidth=1.5, label=f'Mean |E| = {emag_mean:.3f} V/μm', alpha=0.7)
        
        ax_emag_detail.set_xlabel('x (μm)', fontsize=13)
        ax_emag_detail.set_ylabel('Electric Field Magnitude |E| (V/μm)', fontsize=13)
        ax_emag_detail.set_title(f'Electric Field Magnitude Profile at Channel Center (y = {y_coords[y_center_gap_idx]:.2f} μm)', 
                                fontsize=14, fontweight='bold')
        ax_emag_detail.grid(True, linestyle=':', alpha=0.6)
        ax_emag_detail.legend(loc='best', fontsize=11)
        
        # Add annotations for peak field
        idx_max = np.argmax(Emag_profile)
        x_max = x_coords[idx_max]
        ax_emag_detail.annotate(f'Peak: {emag_max:.3f} V/μm\nat x = {x_max:.1f} μm',
                               xy=(x_max, emag_max), xytext=(x_max + 20, emag_max * 0.8),
                               arrowprops=dict(arrowstyle='->', color='black', lw=1.5),
                               fontsize=10, bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.7))
        
        fig_emag_detail.tight_layout()
        plt.savefig(os.path.join(output_folder_name, "efield_magnitude_detail_plot.png"), dpi=300)
        print(f"Detailed E-field magnitude plot saved to {os.path.join(output_folder_name, 'efield_magnitude_detail_plot.png')}")

    # --- Proton Trajectory Analysis ---
    # Initialize h_param_um, it might be set from geo_params if available
    h_param_um = None
    if geo_params and geo_params.get("h") is not None:
        h_param_um = geo_params.get("h") # Assuming 'h' is in μm in geometry_params.csv for this script's context
    elif x_coords is not None and len(x_coords) > 1:
        # Fallback: try to derive from x_coordinates if not in geo_params
        h_param_um = np.abs(x_coords[1] - x_coords[0])
        print(f"Warning: 'h' not found in geometry_params.csv. Using h_param_um derived from x_coordinates: {h_param_um:.2f} μm")
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
