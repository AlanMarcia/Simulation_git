#!/usr/bin/env python3
"""
Script per analizzare i risultati delle simulazioni batch.
Genera statistiche e confronti tra diversi campi elettrici.

Uso: python analyze_batch_results.py [directory_base]
"""

import sys
import os
import csv
from pathlib import Path
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt

def load_trajectories(csv_path):
    """Carica le traiettorie da un file CSV"""
    trajectories = defaultdict(list)
    try:
        with open(csv_path, newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                pid = int(row["proton_id"])
                data = {
                    'time': float(row["time_s"]),
                    'x': float(row["pos_x_m"]),
                    'y': float(row["pos_y_m"]),
                    'vx': float(row["vel_x_m_s"]),
                    'vy': float(row["vel_y_m_s"]),
                }
                trajectories[pid].append(data)
    except Exception as e:
        print(f"Error reading {csv_path}: {e}")
    
    return trajectories

def calculate_statistics(trajectories):
    """Calcola statistiche dalle traiettorie"""
    if not trajectories:
        return {}
    
    stats = {
        'n_protons': len(trajectories),
        'final_x': [],
        'final_y': [],
        'final_vx': [],
        'final_vy': [],
        'deflection_y': [],
        'total_displacement': [],
        'max_deflection_y': [],
    }
    
    for traj in trajectories.values():
        if traj:
            final_point = traj[-1]
            initial_point = traj[0]
            
            stats['final_x'].append(final_point['x'])
            stats['final_y'].append(final_point['y'])
            stats['final_vx'].append(final_point['vx'])
            stats['final_vy'].append(final_point['vy'])
            stats['deflection_y'].append(final_point['y'] - initial_point['y'])
            
            # Calcola lo spostamento totale
            dx = final_point['x'] - initial_point['x']
            dy = final_point['y'] - initial_point['y']
            displacement = np.sqrt(dx**2 + dy**2)
            stats['total_displacement'].append(displacement)
            
            # Calcola la massima deflessione in y
            y_values = [p['y'] for p in traj]
            max_defl = max(abs(max(y_values) - initial_point['y']), 
                          abs(min(y_values) - initial_point['y']))
            stats['max_deflection_y'].append(max_defl)
    
    # Converte liste in statistiche
    for key in stats:
        if isinstance(stats[key], list) and len(stats[key]) > 0:
            arr = np.array(stats[key])
            stats[key + '_mean'] = np.mean(arr)
            stats[key + '_std'] = np.std(arr)
            stats[key + '_min'] = np.min(arr)
            stats[key + '_max'] = np.max(arr)
    
    return stats

def print_statistics(field_name, stats):
    """Stampa le statistiche in formato leggibile"""
    print(f"\n{'='*60}")
    print(f"  {field_name}")
    print(f"{'='*60}")
    
    if not stats:
        print("No data available")
        return
    
    print(f"Number of protons: {stats.get('n_protons', 'N/A')}")
    print()
    
    print("Final Position (X):")
    print(f"  Mean: {stats.get('final_x_mean', 'N/A'):.6e} m")
    print(f"  Std:  {stats.get('final_x_std', 'N/A'):.6e} m")
    print()
    
    print("Deflection in Y:")
    print(f"  Mean: {stats.get('deflection_y_mean', 'N/A')*1e6:.3f} μm")
    print(f"  Std:  {stats.get('deflection_y_std', 'N/A')*1e6:.3f} μm")
    print(f"  Min:  {stats.get('deflection_y_min', 'N/A')*1e6:.3f} μm")
    print(f"  Max:  {stats.get('deflection_y_max', 'N/A')*1e6:.3f} μm")
    print()
    
    print("Maximum Deflection in Y:")
    print(f"  Mean: {stats.get('max_deflection_y_mean', 'N/A')*1e6:.3f} μm")
    print(f"  Max:  {stats.get('max_deflection_y_max', 'N/A')*1e6:.3f} μm")
    print()
    
    print("Final Velocity (X):")
    print(f"  Mean: {stats.get('final_vx_mean', 'N/A'):.6e} m/s")
    print(f"  Std:  {stats.get('final_vx_std', 'N/A'):.6e} m/s")
    print()
    
    print("Total Displacement:")
    print(f"  Mean: {stats.get('total_displacement_mean', 'N/A')*1e6:.3f} μm")

def plot_comparison(results, output_dir):
    """Crea grafici di confronto tra i diversi campi"""
    
    field_names = []
    deflections = []
    energies = []
    
    for field_name, stats in sorted(results.items()):
        if stats:
            field_names.append(field_name)
            deflections.append(stats.get('deflection_y_max', 0) * 1e6)
            
            # Stima l'energia del campo dal nome se possibile
            # Es: field_neg100V -> 100V
            try:
                voltage = int(field_name.replace('field_neg', '').replace('V', ''))
                energies.append(voltage)
            except:
                energies.append(0)
    
    if not field_names:
        print("No data to plot")
        return
    
    # Ordina per voltaggio
    if any(energies):
        sorted_data = sorted(zip(energies, field_names, deflections))
        energies_sorted, field_names_sorted, deflections_sorted = zip(*sorted_data)
    else:
        field_names_sorted = field_names
        deflections_sorted = deflections
    
    # Crea figura
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 5))
    
    # Plot 1: Deflessione vs Campo
    ax1.bar(range(len(field_names_sorted)), deflections_sorted, color='steelblue', alpha=0.7)
    ax1.set_xlabel("Field Configuration", fontsize=11)
    ax1.set_ylabel("Max Deflection in Y [μm]", fontsize=11)
    ax1.set_title("Maximum Y Deflection Comparison", fontsize=12, fontweight='bold')
    ax1.set_xticks(range(len(field_names_sorted)))
    ax1.set_xticklabels(field_names_sorted, rotation=45, ha='right')
    ax1.grid(True, alpha=0.3, axis='y')
    
    # Plot 2: Deflessione vs Voltaggio (se disponibile)
    if any(energies):
        ax2.plot(energies_sorted, deflections_sorted, 'o-', linewidth=2, markersize=8, color='darkred')
        ax2.set_xlabel("Applied Voltage [V]", fontsize=11)
        ax2.set_ylabel("Max Deflection in Y [μm]", fontsize=11)
        ax2.set_title("Deflection vs Voltage", fontsize=12, fontweight='bold')
        ax2.grid(True, alpha=0.3)
    else:
        ax2.text(0.5, 0.5, 'Voltage data not available', 
                ha='center', va='center', transform=ax2.transAxes)
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/batch_comparison.png", dpi=150, bbox_inches='tight')
    print(f"\nComparison plot saved: {output_dir}/batch_comparison.png")
    plt.close()

def main():
    base_output_dir = sys.argv[1] if len(sys.argv) > 1 else "proton_trajectories_output"
    
    if not os.path.exists(base_output_dir):
        print(f"Error: Directory {base_output_dir} does not exist")
        sys.exit(1)
    
    # Trova tutte le sottocartelle
    subdirs = sorted([d for d in os.listdir(base_output_dir) 
                      if os.path.isdir(os.path.join(base_output_dir, d))])
    
    if not subdirs:
        print(f"No subdirectories found in {base_output_dir}")
        sys.exit(1)
    
    print(f"\nAnalyzing {len(subdirs)} simulations...\n")
    
    results = {}
    
    for field_name in subdirs:
        csv_path = os.path.join(base_output_dir, field_name, "proton_trajectories.csv")
        trajectories = load_trajectories(csv_path)
        stats = calculate_statistics(trajectories)
        results[field_name] = stats
        print_statistics(field_name, stats)
    
    # Genera grafici di confronto
    plot_comparison(results, base_output_dir)
    
    print("\n" + "="*60)
    print("Analysis complete!")
    print("="*60 + "\n")

if __name__ == "__main__":
    main()
