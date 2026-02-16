#!/usr/bin/env python3
"""
Script per calcolare la lunghezza focale della lente Einzel per ogni voltaggio.
La focale è definita come la posizione x dove le traiettorie dei protoni
attraversano l'asse y (y ≈ 0) dopo essere passate attraverso la lente.

Uso: python calculate_focal_length.py [directory_base]
"""

import sys
import os
import csv
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
                x = float(row["pos_x_m"])
                y = float(row["pos_y_m"])
                trajectories[pid].append((x, y))
    except Exception as e:
        print(f"Error reading {csv_path}: {e}")
    
    return trajectories

def find_focal_point(trajectory, lens_exit_x=0.0):
    """
    Trova il punto focale di una traiettoria (dove attraversa y=0 dopo la lente).
    
    Args:
        trajectory: lista di tuple (x, y)
        lens_exit_x: posizione x approssimativa di uscita dalla lente
    
    Returns:
        x_focal: posizione x del punto focale, o None se non trovato
    """
    # Filtra solo i punti dopo l'uscita dalla lente
    post_lens = [(x, y) for x, y in trajectory if x > lens_exit_x]
    
    if len(post_lens) < 2:
        return None
    
    # Trova dove |y| è minimo (più vicino a zero)
    min_idx = np.argmin([abs(y) for x, y in post_lens])
    
    # Se siamo al bordo, non possiamo interpolare
    if min_idx == 0 or min_idx == len(post_lens) - 1:
        return post_lens[min_idx][0]
    
    # Interpolazione lineare per trovare il punto esatto dove y=0
    x1, y1 = post_lens[min_idx - 1]
    x2, y2 = post_lens[min_idx]
    x3, y3 = post_lens[min_idx + 1]
    
    # Trova il minimo locale di |y|
    # Se y cambia segno, interpola linearmente
    if y1 * y2 < 0:  # Attraversamento tra punto 1 e 2
        x_focal = x1 + (x2 - x1) * abs(y1) / (abs(y1) + abs(y2))
        return x_focal
    elif y2 * y3 < 0:  # Attraversamento tra punto 2 e 3
        x_focal = x2 + (x3 - x2) * abs(y2) / (abs(y2) + abs(y3))
        return x_focal
    else:
        # Nessun attraversamento, usa il punto più vicino
        return x2

def calculate_focal_length_for_field(trajectories, lens_exit_x=0.0):
    """
    Calcola la lunghezza focale media per un campo.
    
    Returns:
        focal_lengths: lista delle lunghezze focali per ogni protone
        mean_focal: lunghezza focale media
        std_focal: deviazione standard
    """
    focal_points = []
    
    for pid, trajectory in trajectories.items():
        if len(trajectory) < 10:  # Traiettoria troppo corta
            continue
        
        x_focal = find_focal_point(trajectory, lens_exit_x)
        if x_focal is not None:
            focal_points.append(x_focal)
    
    if not focal_points:
        return [], None, None
    
    focal_array = np.array(focal_points)
    mean_focal = np.mean(focal_array)
    std_focal = np.std(focal_array)
    
    return focal_points, mean_focal, std_focal

def extract_voltage_from_name(field_name):
    """Estrae il voltaggio dal nome del campo (es. field_neg1000V -> 1000)"""
    try:
        # Rimuovi 'field_neg' e 'V'
        voltage_str = field_name.replace('field_neg', '').replace('V', '')
        return int(voltage_str)
    except:
        return None

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
    
    print(f"\n{'='*70}")
    print(f"  Focal Length Analysis - Einzel Lens")
    print(f"{'='*70}\n")
    
    # Parametro: posizione x approssimativa di uscita dalla lente
    # (puoi modificarlo in base alla tua geometria)
    lens_exit_x = 100e-6  # 100 micrometri
    
    results = {}
    voltages = []
    focal_lengths = []
    focal_stds = []
    
    for field_name in subdirs:
        csv_path = os.path.join(base_output_dir, field_name, "proton_trajectories.csv")
        trajectories = load_trajectories(csv_path)
        
        if not trajectories:
            print(f"Warning: No trajectories for {field_name}")
            continue
        
        focal_points, mean_focal, std_focal = calculate_focal_length_for_field(
            trajectories, lens_exit_x
        )
        
        voltage = extract_voltage_from_name(field_name)
        
        if mean_focal is not None:
            results[field_name] = {
                'voltage': voltage,
                'focal_points': focal_points,
                'mean_focal': mean_focal,
                'std_focal': std_focal,
                'n_protons': len(focal_points)
            }
            
            if voltage is not None:
                voltages.append(voltage)
                focal_lengths.append(mean_focal * 1e6)  # Converti in micrometri
                focal_stds.append(std_focal * 1e6 if std_focal else 0)
            
            print(f"{field_name}:")
            print(f"  Voltage: {voltage} V" if voltage else f"  Field: {field_name}")
            print(f"  Focal length: {mean_focal*1e6:.2f} ± {std_focal*1e6:.2f} μm")
            print(f"  Number of protons: {len(focal_points)}")
            print()
    
    # Ordina per voltaggio
    if voltages:
        sorted_data = sorted(zip(voltages, focal_lengths, focal_stds))
        voltages_sorted, focal_lengths_sorted, focal_stds_sorted = zip(*sorted_data)
        
        # Crea grafico
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))
        
        # Grafico 1: Lunghezza focale vs voltaggio
        ax1.errorbar(voltages_sorted, focal_lengths_sorted, yerr=focal_stds_sorted,
                     fmt='o-', linewidth=2, markersize=8, capsize=5, capthick=2,
                     color='darkblue', ecolor='lightblue', label='Focal length')
        ax1.set_xlabel("Applied Voltage [V]", fontsize=12, fontweight='bold')
        ax1.set_ylabel("Focal Length [μm]", fontsize=12, fontweight='bold')
        ax1.set_title("Focal Length vs Applied Voltage", fontsize=14, fontweight='bold')
        ax1.grid(True, alpha=0.3, linestyle='--')
        ax1.legend()
        
        # Grafico 2: Potere diottrico (1/f) vs voltaggio
        focal_lengths_m = np.array(focal_lengths_sorted) * 1e-6
        diopters = 1.0 / focal_lengths_m  # In 1/m
        
        ax2.plot(voltages_sorted, diopters, 'o-', linewidth=2, markersize=8,
                color='darkred', label='Optical power')
        ax2.set_xlabel("Applied Voltage [V]", fontsize=12, fontweight='bold')
        ax2.set_ylabel("Optical Power [1/m]", fontsize=12, fontweight='bold')
        ax2.set_title("Optical Power vs Applied Voltage", fontsize=14, fontweight='bold')
        ax2.grid(True, alpha=0.3, linestyle='--')
        ax2.legend()
        
        plt.tight_layout()
        plt.savefig("focal_length_analysis.png", dpi=300, bbox_inches='tight')
        print(f"\n{'='*70}")
        print(f"✓ Focal length plot saved: focal_length_analysis.png")
        print(f"{'='*70}\n")
        
        # Salva anche i dati in CSV
        with open("focal_length_data.csv", "w", newline='') as f:
            writer = csv.writer(f)
            writer.writerow(["Voltage_V", "Focal_Length_um", "Std_Dev_um", "Optical_Power_1_per_m"])
            for v, f_len, f_std in zip(voltages_sorted, focal_lengths_sorted, focal_stds_sorted):
                f_m = f_len * 1e-6
                diopter = 1.0 / f_m
                writer.writerow([v, f_len, f_std, diopter])
        
        print(f"✓ Focal length data saved: focal_length_data.csv\n")
        
        # Mostra grafico
        plt.show()
    else:
        print("Warning: No voltage information found in field names")

if __name__ == "__main__":
    main()
