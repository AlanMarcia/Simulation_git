#!/usr/bin/env python3
"""
Plot FDTD results from C++ simulation
Reads CSV files and generates visualization plots
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.gridspec import GridSpec
import os

# Directory with results
RESULTS_DIR = "fdtd_results"

# Geometry parameters (must match C++ code)
WAVELENGTH_UM = 2.0
L_TOTAL = 360.0
X_FREE_SPACE = 30.0
X_STRUCTURE_LEN = 300.0
Y_VACUUM_PADDING_BOTTOM = 20.0
Y_VACUUM_PADDING_TOP = 20.0
Y_SI_BASE_HEIGHT = 100.0
Y_TOOTH_HEIGHT = 30.0
Y_VACUUM_GAP_THICK = 30.0

H_TOTAL = (Y_VACUUM_PADDING_BOTTOM + 
           Y_SI_BASE_HEIGHT * 2.0 + 
           Y_TOOTH_HEIGHT * 2.0 + 
           Y_VACUUM_GAP_THICK + 
           Y_VACUUM_PADDING_TOP)

POINTS_PER_WAVELENGTH = 20
DX_UM = WAVELENGTH_UM / POINTS_PER_WAVELENGTH
DY_UM = DX_UM

NX = int(L_TOTAL / DX_UM)
NY = int(H_TOTAL / DY_UM)

print("="*70)
print("FDTD RESULTS VISUALIZATION")
print("="*70)
print(f"\nExpected grid size: {NX} x {NY}")
print(f"Resolution: {DX_UM} µm")

# Create coordinate arrays
x_coords = np.linspace(0, L_TOTAL, NX)
y_coords = np.linspace(0, H_TOTAL, NY)

# ============================================================================
# LOAD DATA
# ============================================================================
print("\nLoading data files...")

# Load permittivity map
eps_file = os.path.join(RESULTS_DIR, "permittivity_map.csv")
if os.path.exists(eps_file):
    eps_r = np.loadtxt(eps_file, delimiter=',').T  # Transpose to get (NX, NY)
    print(f"✓ Loaded permittivity map: {eps_r.shape}")
else:
    print(f"✗ Warning: {eps_file} not found")
    eps_r = np.ones((NX, NY))

# Load snapshots
snapshots = []
snapshot_times = []
for i in range(4):
    snapshot_file = os.path.join(RESULTS_DIR, f"Ez_snapshot_{i}.csv")
    if os.path.exists(snapshot_file):
        Ez = np.loadtxt(snapshot_file, delimiter=',').T  # Transpose to get (NX, NY)
        snapshots.append(Ez)
        print(f"✓ Loaded snapshot {i}: {Ez.shape}, max |Ez| = {np.max(np.abs(Ez)):.3e} V/m")
    else:
        print(f"✗ Warning: {snapshot_file} not found")

# Load final field
final_file = os.path.join(RESULTS_DIR, "Ez_final.csv")
if os.path.exists(final_file):
    Ez_final = np.loadtxt(final_file, delimiter=',').T  # Transpose to get (NX, NY)
    print(f"✓ Loaded final field: {Ez_final.shape}, max |Ez| = {np.max(np.abs(Ez_final)):.3e} V/m")
else:
    print(f"✗ Warning: {final_file} not found")
    Ez_final = None

# Load field along gap
gap_file = os.path.join(RESULTS_DIR, "field_along_gap.csv")
if os.path.exists(gap_file):
    gap_data = np.loadtxt(gap_file, delimiter=',', skiprows=1)
    x_gap = gap_data[:, 0]
    Ez_gap = gap_data[:, 1]
    print(f"✓ Loaded gap field: {len(x_gap)} points")
else:
    print(f"✗ Warning: {gap_file} not found")
    x_gap = None
    Ez_gap = None

# Load diagnostics
diag_file = os.path.join(RESULTS_DIR, "diagnostics.csv")
if os.path.exists(diag_file):
    diag_data = np.loadtxt(diag_file, delimiter=',', skiprows=1)
    time_ps = diag_data[:, 0]
    max_Ez_vs_time = diag_data[:, 1]
    print(f"✓ Loaded diagnostics: {len(time_ps)} time points")
else:
    print(f"✗ Warning: {diag_file} not found")
    time_ps = None
    max_Ez_vs_time = None

# ============================================================================
# PLOT 1: GEOMETRY
# ============================================================================
print("\nGenerating plots...")

fig1 = plt.figure(figsize=(14, 6))
plt.imshow(eps_r.T, extent=[0, L_TOTAL, 0, H_TOTAL], 
           origin='lower', cmap='RdYlBu_r', aspect='auto')
plt.colorbar(label='Relative Permittivity $\\epsilon_r$')
plt.xlabel('x [µm]')
plt.ylabel('y [µm]')
plt.title('FDTD Geometry - Denti Uguali Structure')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(RESULTS_DIR, 'geometry_permittivity.png'), dpi=150)
print("✓ Saved: geometry_permittivity.png")

# ============================================================================
# PLOT 2: FIELD SNAPSHOTS
# ============================================================================
if snapshots:
    fig2 = plt.figure(figsize=(16, 10))
    snapshot_labels = ['25%', '50%', '75%', '100%']
    
    # Find global max for consistent colorscale
    vmax = max([np.max(np.abs(snap)) for snap in snapshots])
    
    for idx, Ez in enumerate(snapshots):
        plt.subplot(2, 2, idx+1)
        im = plt.imshow(Ez.T, extent=[0, L_TOTAL, 0, H_TOTAL],
                       origin='lower', cmap='RdBu', aspect='auto',
                       vmin=-vmax, vmax=vmax)
        plt.colorbar(im, label='$E_z$ [V/m]')
        plt.xlabel('x [µm]')
        plt.ylabel('y [µm]')
        plt.title(f'Electric Field $E_z$ at {snapshot_labels[idx]} of simulation\nmax |$E_z$| = {np.max(np.abs(Ez)):.3e} V/m')
        plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, 'field_snapshots.png'), dpi=150)
    print("✓ Saved: field_snapshots.png")

# ============================================================================
# PLOT 3: FINAL FIELD
# ============================================================================
if Ez_final is not None:
    fig3 = plt.figure(figsize=(16, 7))
    
    vmax = np.max(np.abs(Ez_final))
    im = plt.imshow(Ez_final.T, extent=[0, L_TOTAL, 0, H_TOTAL],
                   origin='lower', cmap='RdBu', aspect='auto',
                   vmin=-vmax, vmax=vmax)
    
    # Overlay geometry contour
    X, Y = np.meshgrid(x_coords, y_coords)
    plt.contour(X, Y, eps_r.T, levels=[1.5], 
               colors='black', linewidths=1, alpha=0.5)
    
    plt.colorbar(im, label='$E_z$ [V/m]')
    plt.xlabel('x [µm]')
    plt.ylabel('y [µm]')
    plt.title(f'Final Electric Field $E_z$ Distribution\nmax |$E_z$| = {vmax:.3e} V/m')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, 'field_final.png'), dpi=150)
    print("✓ Saved: field_final.png")

# ============================================================================
# PLOT 4: FIELD ALONG GAP
# ============================================================================
if x_gap is not None and Ez_gap is not None:
    fig4 = plt.figure(figsize=(14, 6))
    
    plt.subplot(2, 1, 1)
    plt.plot(x_gap, Ez_gap / 1e6, 'b-', linewidth=1)
    plt.xlabel('x [µm]')
    plt.ylabel('$E_z$ [MV/m]')
    plt.title('Electric Field Along Vacuum Gap Center')
    plt.grid(True, alpha=0.3)
    plt.axvline(X_FREE_SPACE, color='r', linestyle='--', alpha=0.5, label='Structure start')
    plt.axvline(X_FREE_SPACE + X_STRUCTURE_LEN, color='r', linestyle='--', alpha=0.5, label='Structure end')
    plt.legend()
    
    plt.subplot(2, 1, 2)
    plt.plot(x_gap, np.abs(Ez_gap) / 1e6, 'r-', linewidth=1)
    plt.xlabel('x [µm]')
    plt.ylabel('|$E_z$| [MV/m]')
    plt.title('Electric Field Magnitude Along Vacuum Gap Center')
    plt.grid(True, alpha=0.3)
    plt.axvline(X_FREE_SPACE, color='r', linestyle='--', alpha=0.5)
    plt.axvline(X_FREE_SPACE + X_STRUCTURE_LEN, color='r', linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, 'field_along_gap.png'), dpi=150)
    print("✓ Saved: field_along_gap.png")

# ============================================================================
# PLOT 5: DIAGNOSTICS - MAX FIELD VS TIME
# ============================================================================
if time_ps is not None and max_Ez_vs_time is not None:
    fig5 = plt.figure(figsize=(14, 6))
    
    plt.subplot(1, 2, 1)
    plt.plot(time_ps, max_Ez_vs_time / 1e6, 'b-', linewidth=1.5)
    plt.xlabel('Time [ps]')
    plt.ylabel('max |$E_z$| [MV/m]')
    plt.title('Maximum Electric Field vs Time')
    plt.grid(True, alpha=0.3)
    
    plt.subplot(1, 2, 2)
    plt.semilogy(time_ps, max_Ez_vs_time, 'r-', linewidth=1.5)
    plt.xlabel('Time [ps]')
    plt.ylabel('max |$E_z$| [V/m]')
    plt.title('Maximum Electric Field vs Time (log scale)')
    plt.grid(True, alpha=0.3, which='both')
    
    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, 'field_max_vs_time.png'), dpi=150)
    print("✓ Saved: field_max_vs_time.png")

# ============================================================================
# PLOT 6: DETAILED VIEW OF GAP REGION
# ============================================================================
if Ez_final is not None:
    fig6 = plt.figure(figsize=(16, 8))
    
    # Define gap region
    y_gap_center = H_TOTAL / 2.0
    y_gap_range = 80.0  # Show ±40 µm around gap center
    x_structure_range = [X_FREE_SPACE - 10, X_FREE_SPACE + X_STRUCTURE_LEN + 10]
    
    # Find indices - Note: after transpose, Ez_final has shape (NX, NY)
    x_indices = (x_coords >= x_structure_range[0]) & (x_coords <= x_structure_range[1])
    y_indices = (y_coords >= y_gap_center - y_gap_range/2) & (y_coords <= y_gap_center + y_gap_range/2)
    
    # Extract zoom region
    Ez_zoom = Ez_final[x_indices, :][:, y_indices]
    x_zoom = x_coords[x_indices]
    y_zoom = y_coords[y_indices]
    
    vmax_zoom = np.max(np.abs(Ez_zoom))
    
    plt.subplot(1, 2, 1)
    im = plt.imshow(Ez_zoom.T, extent=[x_zoom[0], x_zoom[-1], y_zoom[0], y_zoom[-1]],
                   origin='lower', cmap='RdBu', aspect='auto',
                   vmin=-vmax_zoom, vmax=vmax_zoom)
    
    # Overlay geometry
    eps_zoom = eps_r[x_indices, :][:, y_indices]
    X_zoom, Y_zoom = np.meshgrid(x_zoom, y_zoom)
    plt.contour(X_zoom, Y_zoom, eps_zoom.T, levels=[1.5], 
               colors='black', linewidths=1.5)
    
    plt.colorbar(im, label='$E_z$ [V/m]')
    plt.xlabel('x [µm]')
    plt.ylabel('y [µm]')
    plt.title(f'Zoomed View: Gap Region\nmax |$E_z$| = {vmax_zoom:.3e} V/m')
    plt.grid(True, alpha=0.3)
    
    plt.subplot(1, 2, 2)
    # Plot field magnitude with log scale
    Ez_zoom_abs = np.abs(Ez_zoom)
    vmin_log = max(1e2, np.min(Ez_zoom_abs[Ez_zoom_abs > 0]) if np.any(Ez_zoom_abs > 0) else 1e2)
    vmax_log = max(vmin_log * 10, vmax_zoom)  # Ensure vmax > vmin
    
    im2 = plt.imshow(Ez_zoom_abs.T, extent=[x_zoom[0], x_zoom[-1], y_zoom[0], y_zoom[-1]],
                    origin='lower', cmap='hot', aspect='auto',
                    norm=colors.LogNorm(vmin=vmin_log, vmax=vmax_log))
    
    plt.contour(X_zoom, Y_zoom, eps_zoom.T, levels=[1.5], 
               colors='cyan', linewidths=1.5, alpha=0.7)
    
    plt.colorbar(im2, label='|$E_z$| [V/m]')
    plt.xlabel('x [µm]')
    plt.ylabel('y [µm]')
    plt.title('Field Magnitude (log scale)')
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, 'field_gap_zoom.png'), dpi=150)
    print("✓ Saved: field_gap_zoom.png")

# ============================================================================
# SUMMARY
# ============================================================================
print("\n" + "="*70)
print("VISUALIZATION COMPLETED")
print("="*70)
print(f"\nAll plots saved in: {RESULTS_DIR}/")
print("\nGenerated files:")
print("  1. geometry_permittivity.png  - Structure geometry")
print("  2. field_snapshots.png        - Field evolution (4 snapshots)")
print("  3. field_final.png            - Final field distribution")
print("  4. field_along_gap.png        - Field profile along gap")
print("  5. field_max_vs_time.png      - Maximum field vs time")
print("  6. field_gap_zoom.png         - Detailed gap region view")

if Ez_final is not None:
    print(f"\nFinal field statistics:")
    print(f"  Max |Ez| = {np.max(np.abs(Ez_final)):.3e} V/m")
    print(f"  Mean |Ez| = {np.mean(np.abs(Ez_final)):.3e} V/m")
    
    # Calculate field in gap
    y_gap_idx = int((H_TOTAL / 2.0) / DY_UM)
    x_struct_start = int(X_FREE_SPACE / DX_UM)
    x_struct_end = int((X_FREE_SPACE + X_STRUCTURE_LEN) / DX_UM)
    
    Ez_in_gap = Ez_final[x_struct_start:x_struct_end, y_gap_idx]
    print(f"  Max |Ez| in gap = {np.max(np.abs(Ez_in_gap)):.3e} V/m")
    print(f"  Mean |Ez| in gap = {np.mean(np.abs(Ez_in_gap)):.3e} V/m")

print("\n" + "="*70)
