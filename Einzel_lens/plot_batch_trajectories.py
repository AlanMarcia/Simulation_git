#!/usr/bin/env python3
"""
Script per visualizzare le traiettorie di protoni da tutte le simulazioni.
Uso: python plot_batch_trajectories.py [directory_base]
"""

import sys
import os
import csv
from collections import defaultdict
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def load_trajectories(csv_path):
    """Carica le traiettorie da un file CSV"""
    trajectories = defaultdict(list)
    try:
        with open(csv_path, newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            required = {"proton_id", "pos_x_m", "pos_y_m"}
            if not required.issubset(reader.fieldnames or {}):
                return trajectories
            for row in reader:
                pid = int(row["proton_id"])
                x = float(row["pos_x_m"])
                y = float(row["pos_y_m"])
                trajectories[pid].append((x, y))
    except Exception as e:
        print(f"Error reading {csv_path}: {e}")
    
    return trajectories

def plot_field_trajectories(output_dir, field_name):
    """Plotta le traiettorie per un singolo campo"""
    csv_path = os.path.join(output_dir, field_name, "proton_trajectories.csv")
    
    if not os.path.exists(csv_path):
        print(f"Warning: {csv_path} not found")
        return None
    
    trajectories = load_trajectories(csv_path)
    
    if not trajectories:
        print(f"No trajectories found in {csv_path}")
        return None
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    for pid in sorted(trajectories.keys()):
        pts = trajectories[pid]
        xs = [p[0] * 1e6 for p in pts]
        ys = [p[1] * 1e6 for p in pts]
        ax.plot(xs, ys, linewidth=1.0, alpha=0.7)
    
    ax.set_xlabel("x [μm]", fontsize=10)
    ax.set_ylabel("y [μm]", fontsize=10)
    ax.set_title(f"Trajectories - {field_name}", fontsize=12, fontweight='bold')
    ax.grid(True, linestyle="--", alpha=0.3)
    ax.set_aspect('equal', adjustable='box')
    
    return fig

def main():
    base_output_dir = sys.argv[1] if len(sys.argv) > 1 else "proton_trajectories_output"
    
    if not os.path.exists(base_output_dir):
        print(f"Error: Directory {base_output_dir} does not exist")
        sys.exit(1)
    
    # Trova tutte le sottocartelle (una per ogni campo)
    subdirs = sorted([d for d in os.listdir(base_output_dir) 
                      if os.path.isdir(os.path.join(base_output_dir, d))])
    
    if not subdirs:
        print(f"No subdirectories found in {base_output_dir}")
        sys.exit(1)
    
    print(f"Found {len(subdirs)} field simulations")
    
    # Crea un grafico separato per ogni campo
    for field_name in subdirs:
        csv_path = os.path.join(base_output_dir, field_name, "proton_trajectories.csv")
        trajectories = load_trajectories(csv_path)
        
        if not trajectories:
            print(f"Warning: No trajectories found for {field_name}")
            continue
        
        # Crea figura per questo campo
        fig, ax = plt.subplots(figsize=(10, 8))
        
        for pid in sorted(trajectories.keys()):
            pts = trajectories[pid]
            xs = [p[0] * 1e6 for p in pts]
            ys = [p[1] * 1e6 for p in pts]
            ax.plot(xs, ys, linewidth=1.0, alpha=0.7)
        
        ax.set_xlabel("x [μm]", fontsize=12)
        ax.set_ylabel("y [μm]", fontsize=12)
        ax.set_title(f"Proton Trajectories - {field_name}", fontsize=14, fontweight='bold')
        ax.grid(True, linestyle="--", alpha=0.3)
        ax.set_aspect('equal', adjustable='box')
        
        plt.tight_layout()
        
        # Salva il grafico con nome del campo
        output_filename = f"{field_name}_trajectories.png"
       # plt.savefig(output_filename, dpi=300, bbox_inches='tight')
        print(f"Saved: {output_filename}")
        plt.show()
        
        plt.close()
    
    print(f"\nProcessing complete! {len(subdirs)} plot(s) generated.")
    
    # Stampa un riepilogo
    print("\nSummary:")
    for field_name in subdirs:
        csv_path = os.path.join(base_output_dir, field_name, "proton_trajectories.csv")
        trajectories = load_trajectories(csv_path)
        print(f"  {field_name}: {len(trajectories)} protons")

if __name__ == "__main__":
    main()
