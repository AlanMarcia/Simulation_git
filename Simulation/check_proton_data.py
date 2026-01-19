#!/usr/bin/env python3
"""
Script to check proton initial positions and final energies
"""
import pandas as pd
import numpy as np
import sys

# Constants
M_PROTON = 1.67262192e-27  # kg
Q_PROTON = 1.60217663e-19  # C
KEV_TO_JOULES = 1.60217663e-16  # 1 keV in Joules

def analyze_proton_data(folder):
    # Load trajectory data
    csv_file = f"{folder}/all_proton_trajectories.csv"
    print(f"\nAnalyzing {csv_file}...")
    
    df = pd.read_csv(csv_file)
    
    # Get initial positions (time_s = 0)
    initial_data = df[df['time_s'] == 0.0].copy()
    
    print("\n=== INITIAL POSITIONS (at t=0) ===")
    print(f"Number of protons: {len(initial_data)}")
    print(f"X position: {initial_data['x_m'].iloc[0] * 1e6:.2f} um (all same)")
    print(f"Y range: [{initial_data['y_m'].min() * 1e6:.4f}, {initial_data['y_m'].max() * 1e6:.4f}] um")
    print(f"Initial velocity vx: {initial_data['vx_m_per_s'].iloc[0]:.3e} m/s")
    
    # Calculate initial kinetic energy
    v_initial = np.sqrt(initial_data['vx_m_per_s']**2 + initial_data['vy_m_per_s']**2)
    E_initial_J = 0.5 * M_PROTON * v_initial**2
    E_initial_keV = E_initial_J / KEV_TO_JOULES
    print(f"Initial kinetic energy: {E_initial_keV.iloc[0]:.2f} keV")
    
    # Get final positions for each proton (last entry for each proton_id)
    final_data = df.groupby('proton_id').last().reset_index()
    
    print("\n=== FINAL STATE ===")
    print(f"Protons analyzed: {len(final_data)}")
    
    # Calculate final kinetic energy
    v_final = np.sqrt(final_data['vx_m_per_s']**2 + final_data['vy_m_per_s']**2)
    E_final_J = 0.5 * M_PROTON * v_final**2
    E_final_keV = E_final_J / KEV_TO_JOULES
    
    # Filter successful protons (reached the end)
    successful = final_data[final_data['x_m'] >= df['x_m'].max() * 0.95]
    
    if len(successful) > 0:
        v_success = np.sqrt(successful['vx_m_per_s']**2 + successful['vy_m_per_s']**2)
        E_success_J = 0.5 * M_PROTON * v_success**2
        E_success_keV = E_success_J / KEV_TO_JOULES
        
        print(f"\nSuccessful protons (reached end): {len(successful)}")
        print(f"Final energy range: [{E_success_keV.min():.2f}, {E_success_keV.max():.2f}] keV")
        print(f"Mean final energy: {E_success_keV.mean():.2f} keV")
        print(f"Median final energy: {E_success_keV.median():.2f} keV")
        
        # Check if any exceed theoretical maximum
        E_gained = E_success_keV - E_initial_keV.iloc[0]
        print(f"\nEnergy gained: [{E_gained.min():.2f}, {E_gained.max():.2f}] keV")
        print(f"Mean energy gained: {E_gained.mean():.2f} keV")
        
        if E_success_keV.max() > 150:
            print(f"\n⚠️  WARNING: Some protons exceed 150 keV!")
            print(f"   With 40 keV initial + 100 kV potential, max should be ~140 keV")
            print(f"   Maximum reached: {E_success_keV.max():.2f} keV")
        else:
            print(f"\n✓ Energy conservation looks good (all under 150 keV)")
    else:
        print("No successful protons found!")
    
    # Final position statistics
    print(f"\nFinal X range: [{final_data['x_m'].min() * 1e6:.2f}, {final_data['x_m'].max() * 1e6:.2f}] um")
    print(f"Final Y range: [{final_data['y_m'].min() * 1e6:.2f}, {final_data['y_m'].max() * 1e6:.2f}] um")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        folder = sys.argv[1]
    else:
        folder = "geometria_Denti_sfasati_profondi"
    
    analyze_proton_data(folder)
