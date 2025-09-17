#!/usr/bin/env python3
"""
Script per analizzare i risultati della simulazione PIC DLA
Legge i file CSV generati dal simulatore PIC e crea grafici e analisi
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.animation import FuncAnimation
import seaborn as sns
import glob
import os
import argparse
from pathlib import Path

# Configurazione stile grafici
plt.style.use('seaborn-v0_8')
sns.set_palette("husl")

class PICAnalyzer:
    def __init__(self, data_dir="."):
        """
        Inizializza l'analizzatore PIC
        
        Args:
            data_dir (str): Directory contenente i file CSV
        """
        self.data_dir = Path(data_dir)
        self.field_files = []
        self.particle_files = []
        self.time_steps = []
        
        # Scansiona i file disponibili
        self.scan_files()
        
    def scan_files(self):
        """Scansiona la directory per trovare tutti i file PIC"""
        # File dei campi
        field_pattern = str(self.data_dir / "pic_dla_*_fields.csv")
        self.field_files = sorted(glob.glob(field_pattern))
        
        # File delle particelle
        particle_pattern = str(self.data_dir / "pic_dla_*_particles.csv")
        self.particle_files = sorted(glob.glob(particle_pattern))
        
        # Estrai i time steps
        for file in self.particle_files:
            if "step_" in file:
                step = file.split("step_")[1].split("_")[0]
                if step.isdigit():
                    self.time_steps.append(int(step))
        
        self.time_steps = sorted(self.time_steps)
        
        print(f"Trovati {len(self.field_files)} file di campi")
        print(f"Trovati {len(self.particle_files)} file di particelle")
        print(f"Time steps disponibili: {self.time_steps}")
        
    def load_fields(self, time_step=None, final=False):
        """
        Carica i dati dei campi elettromagnetici
        
        Args:
            time_step (int): Specifico time step da caricare
            final (bool): Se True, carica il file finale
            
        Returns:
            pandas.DataFrame: Dati dei campi
        """
        if final:
            filename = str(self.data_dir / "pic_dla_final_fields.csv")
        elif time_step is not None:
            filename = str(self.data_dir / f"pic_dla_step_{time_step}_fields.csv")
        else:
            # Prendi il primo file disponibile
            filename = self.field_files[0] if self.field_files else None
            
        if filename and os.path.exists(filename):
            df = pd.read_csv(filename)
            print(f"Caricato file campi: {filename}")
            print(f"Dimensioni griglia: {len(df['x'].unique())} x {len(df['y'].unique())}")
            return df
        else:
            print(f"File non trovato: {filename}")
            return None
            
    def load_particles(self, time_step=None, final=False):
        """
        Carica i dati delle particelle
        
        Args:
            time_step (int): Specifico time step da caricare
            final (bool): Se True, carica il file finale
            
        Returns:
            pandas.DataFrame: Dati delle particelle
        """
        if final:
            filename = str(self.data_dir / "pic_dla_final_particles.csv")
        elif time_step is not None:
            filename = str(self.data_dir / f"pic_dla_step_{time_step}_particles.csv")
        else:
            # Prendi il primo file disponibile
            filename = self.particle_files[0] if self.particle_files else None
            
        if filename and os.path.exists(filename):
            df = pd.read_csv(filename)
            print(f"Caricato file particelle: {filename}")
            print(f"Particelle totali: {len(df)}")
            print(f"Particelle attive: {df['active'].sum()}")
            return df
        else:
            print(f"File non trovato: {filename}")
            return None
    
    def plot_geometry(self, fields_df):
        """
        Visualizza la geometria del dispositivo DLA
        
        Args:
            fields_df (DataFrame): Dati dei campi contenenti eps_r
        """
        if fields_df is None:
            print("Nessun dato dei campi disponibile")
            return
            
        fig, ax = plt.subplots(1, 1, figsize=(15, 6))
        
        # Reshape per la visualizzazione 2D
        x_unique = sorted(fields_df['x'].unique())
        y_unique = sorted(fields_df['y'].unique())
        X, Y = np.meshgrid(x_unique, y_unique)
        
        eps_matrix = np.zeros((len(y_unique), len(x_unique)))
        for i, y in enumerate(y_unique):
            for j, x in enumerate(x_unique):
                eps_val = fields_df[(fields_df['x'] == x) & (fields_df['y'] == y)]['eps_r'].iloc[0]
                eps_matrix[i, j] = eps_val
        
        # Plot della geometria
        im = ax.imshow(eps_matrix, extent=[min(x_unique), max(x_unique), 
                                         min(y_unique), max(y_unique)], 
                      cmap='coolwarm', aspect='auto', origin='lower')
        
        ax.set_xlabel('Posizione x [μm]')
        ax.set_ylabel('Posizione y [μm]')
        ax.set_title('Geometria DLA - Distribuzione Permittività')
        
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Permittività relativa ε_r')
        
        plt.tight_layout()
        plt.savefig(self.data_dir / 'geometria_dla.png', dpi=300, bbox_inches='tight')
        plt.show()
        
    def plot_electric_field(self, fields_df):
        """
        Visualizza il campo elettrico
        
        Args:
            fields_df (DataFrame): Dati dei campi
        """
        if fields_df is None:
            print("Nessun dato dei campi disponibile")
            return
            
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        
        x_unique = sorted(fields_df['x'].unique())
        y_unique = sorted(fields_df['y'].unique())
        X, Y = np.meshgrid(x_unique, y_unique)
        
        # Campo Ex
        Ex_matrix = np.zeros((len(y_unique), len(x_unique)))
        Ey_matrix = np.zeros((len(y_unique), len(x_unique)))
        E_mag_matrix = np.zeros((len(y_unique), len(x_unique)))
        
        for i, y in enumerate(y_unique):
            for j, x in enumerate(x_unique):
                row = fields_df[(fields_df['x'] == x) & (fields_df['y'] == y)]
                if not row.empty:
                    Ex_matrix[i, j] = row['Ex'].iloc[0]
                    Ey_matrix[i, j] = row['Ey'].iloc[0]
                    E_mag_matrix[i, j] = np.sqrt(Ex_matrix[i, j]**2 + Ey_matrix[i, j]**2)
        
        # Plot Ex
        im1 = axes[0,0].imshow(Ex_matrix, extent=[min(x_unique), max(x_unique), 
                                                 min(y_unique), max(y_unique)], 
                              cmap='RdBu', aspect='auto', origin='lower')
        axes[0,0].set_title('Campo Elettrico Ex [V/m]')
        axes[0,0].set_xlabel('x [μm]')
        axes[0,0].set_ylabel('y [μm]')
        plt.colorbar(im1, ax=axes[0,0])
        
        # Plot Ey
        im2 = axes[0,1].imshow(Ey_matrix, extent=[min(x_unique), max(x_unique), 
                                                 min(y_unique), max(y_unique)], 
                              cmap='RdBu', aspect='auto', origin='lower')
        axes[0,1].set_title('Campo Elettrico Ey [V/m]')
        axes[0,1].set_xlabel('x [μm]')
        axes[0,1].set_ylabel('y [μm]')
        plt.colorbar(im2, ax=axes[0,1])
        
        # Plot magnitudine
        im3 = axes[1,0].imshow(E_mag_matrix, extent=[min(x_unique), max(x_unique), 
                                                    min(y_unique), max(y_unique)], 
                              cmap='viridis', aspect='auto', origin='lower')
        axes[1,0].set_title('Magnitudine Campo Elettrico |E| [V/m]')
        axes[1,0].set_xlabel('x [μm]')
        axes[1,0].set_ylabel('y [μm]')
        plt.colorbar(im3, ax=axes[1,0])
        
        # Plot vettoriale (sottocampionato)
        skip = max(1, len(x_unique) // 30)  # Mostra 1 freccia ogni 30 punti
        X_sub = X[::skip, ::skip]
        Y_sub = Y[::skip, ::skip]
        Ex_sub = Ex_matrix[::skip, ::skip]
        Ey_sub = Ey_matrix[::skip, ::skip]
        
        axes[1,1].quiver(X_sub, Y_sub, Ex_sub, Ey_sub, 
                        scale_units='xy', angles='xy', scale=1e8, alpha=0.7)
        axes[1,1].set_title('Campo Elettrico Vettoriale')
        axes[1,1].set_xlabel('x [μm]')
        axes[1,1].set_ylabel('y [μm]')
        axes[1,1].set_aspect('equal')
        
        plt.tight_layout()
        plt.savefig(self.data_dir / 'campo_elettrico.png', dpi=300, bbox_inches='tight')
        plt.show()
        
    def plot_particle_trajectories(self, final=True):
        """
        Visualizza le traiettorie delle particelle
        
        Args:
            final (bool): Se True, usa il file finale
        """
        # Carica tutti i file delle particelle per le traiettorie
        trajectory_data = {}
        
        if final and os.path.exists(self.data_dir / "pic_dla_final_particles.csv"):
            # Solo snapshot finale
            df = self.load_particles(final=True)
            if df is not None:
                self.plot_single_particle_snapshot(df, "finale")
        else:
            # Traiettorie complete se disponibili
            for step in self.time_steps[:10]:  # Primi 10 step per non sovraccaricare
                df = self.load_particles(time_step=step)
                if df is not None:
                    trajectory_data[step] = df
            
            if trajectory_data:
                self.plot_trajectory_evolution(trajectory_data)
    
    def plot_single_particle_snapshot(self, particles_df, label=""):
        """
        Plot snapshot delle particelle in un singolo momento
        """
        fig, axes = plt.subplots(1, 3, figsize=(18, 6))
        
        # Filtra particelle attive
        active_particles = particles_df[particles_df['active'] == 1]
        inactive_particles = particles_df[particles_df['active'] == 0]
        
        # Plot posizioni
        axes[0].scatter(active_particles['x_um'], active_particles['y_um'], 
                       c='blue', s=1, alpha=0.7, label=f'Attive ({len(active_particles)})')
        axes[0].scatter(inactive_particles['x_um'], inactive_particles['y_um'], 
                       c='red', s=1, alpha=0.3, label=f'Inattive ({len(inactive_particles)})')
        axes[0].set_xlabel('Posizione x [μm]')
        axes[0].set_ylabel('Posizione y [μm]')
        axes[0].set_title(f'Posizioni Particelle {label}')
        axes[0].legend()
        axes[0].grid(True, alpha=0.3)
        
        # Plot velocità
        axes[1].scatter(active_particles['vx_mm_per_s'], active_particles['vy_mm_per_s'], 
                       c=active_particles['energy_keV'], s=2, alpha=0.7, cmap='plasma')
        axes[1].set_xlabel('Velocità vx [Mm/s]')
        axes[1].set_ylabel('Velocità vy [Mm/s]')
        axes[1].set_title(f'Spazio delle Velocità {label}')
        cbar = plt.colorbar(axes[1].collections[0], ax=axes[1])
        cbar.set_label('Energia [keV]')
        
        # Istogramma energie
        if len(active_particles) > 0:
            axes[2].hist(active_particles['energy_keV'], bins=50, alpha=0.7, 
                        color='blue', edgecolor='black')
            axes[2].axvline(active_particles['energy_keV'].mean(), color='red', 
                           linestyle='--', label=f'Media: {active_particles["energy_keV"].mean():.1f} keV')
            axes[2].set_xlabel('Energia [keV]')
            axes[2].set_ylabel('Numero di Particelle')
            axes[2].set_title(f'Distribuzione Energie {label}')
            axes[2].legend()
            axes[2].grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(self.data_dir / f'particelle_{label.replace(" ", "_")}.png', 
                   dpi=300, bbox_inches='tight')
        plt.show()
        
    def plot_trajectory_evolution(self, trajectory_data):
        """
        Plot dell'evoluzione delle traiettorie nel tempo
        """
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        
        times = sorted(trajectory_data.keys())
        colors = plt.cm.viridis(np.linspace(0, 1, len(times)))
        
        for i, (time_step, df) in enumerate(trajectory_data.items()):
            active = df[df['active'] == 1]
            if len(active) == 0:
                continue
                
            color = colors[i]
            alpha = 0.3 + 0.7 * i / len(times)  # Fade in nel tempo
            
            # Posizioni
            axes[0,0].scatter(active['x_um'], active['y_um'], 
                             c=[color], s=1, alpha=alpha, label=f'Step {time_step}')
            
            # Energie nel tempo
            axes[0,1].scatter([time_step] * len(active), active['energy_keV'], 
                             c=[color], s=1, alpha=alpha)
            
            # Velocità longitudinale
            axes[1,0].scatter([time_step] * len(active), active['vx_mm_per_s'], 
                             c=[color], s=1, alpha=alpha)
            
            # Velocità trasversale
            axes[1,1].scatter([time_step] * len(active), active['vy_mm_per_s'], 
                             c=[color], s=1, alpha=alpha)
        
        axes[0,0].set_xlabel('x [μm]')
        axes[0,0].set_ylabel('y [μm]')
        axes[0,0].set_title('Evoluzione Posizioni')
        axes[0,0].grid(True, alpha=0.3)
        
        axes[0,1].set_xlabel('Time Step')
        axes[0,1].set_ylabel('Energia [keV]')
        axes[0,1].set_title('Evoluzione Energetica')
        axes[0,1].grid(True, alpha=0.3)
        
        axes[1,0].set_xlabel('Time Step')
        axes[1,0].set_ylabel('vx [Mm/s]')
        axes[1,0].set_title('Velocità Longitudinale')
        axes[1,0].grid(True, alpha=0.3)
        
        axes[1,1].set_xlabel('Time Step')
        axes[1,1].set_ylabel('vy [Mm/s]')
        axes[1,1].set_title('Velocità Trasversale')
        axes[1,1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(self.data_dir / 'evoluzione_particelle.png', dpi=300, bbox_inches='tight')
        plt.show()
        
    def analyze_energy_gain(self):
        """
        Analizza il guadagno energetico delle particelle
        """
        # Confronta stato iniziale e finale
        initial_df = self.load_particles(time_step=0)
        final_df = self.load_particles(final=True)
        
        if initial_df is None or final_df is None:
            print("Impossibile caricare dati iniziali e finali")
            return
            
        # Statistiche energetiche
        initial_active = initial_df[initial_df['active'] == 1]
        final_active = final_df[final_df['active'] == 1]
        
        fig, axes = plt.subplots(2, 2, figsize=(16, 10))
        
        # Confronto distribuzioni energia
        axes[0,0].hist(initial_active['energy_keV'], bins=30, alpha=0.7, 
                      label=f'Iniziale (μ={initial_active["energy_keV"].mean():.1f} keV)', 
                      color='blue')
        if len(final_active) > 0:
            axes[0,0].hist(final_active['energy_keV'], bins=30, alpha=0.7, 
                          label=f'Finale (μ={final_active["energy_keV"].mean():.1f} keV)', 
                          color='red')
        axes[0,0].set_xlabel('Energia [keV]')
        axes[0,0].set_ylabel('Numero di Particelle')
        axes[0,0].set_title('Confronto Distribuzioni Energetiche')
        axes[0,0].legend()
        axes[0,0].grid(True, alpha=0.3)
        
        # Efficienza di trasporto
        survival_rate = len(final_active) / len(initial_active) * 100
        axes[0,1].bar(['Iniziali', 'Finali Attive', 'Perse'], 
                     [len(initial_active), len(final_active), 
                      len(initial_active) - len(final_active)],
                     color=['blue', 'green', 'red'], alpha=0.7)
        axes[0,1].set_ylabel('Numero di Particelle')
        axes[0,1].set_title(f'Efficienza Trasporto: {survival_rate:.1f}%')
        axes[0,1].grid(True, alpha=0.3)
        
        # Guadagno energetico vs posizione finale
        if len(final_active) > 0:
            energy_gain = final_active['energy_keV'] - initial_active['energy_keV'].mean()
            axes[1,0].scatter(final_active['x_um'], energy_gain, 
                             c=final_active['energy_keV'], s=10, alpha=0.7, cmap='plasma')
            axes[1,0].set_xlabel('Posizione finale x [μm]')
            axes[1,0].set_ylabel('Guadagno energetico [keV]')
            axes[1,0].set_title('Guadagno Energetico vs Posizione')
            axes[1,0].grid(True, alpha=0.3)
            cbar = plt.colorbar(axes[1,0].collections[0], ax=axes[1,0])
            cbar.set_label('Energia finale [keV]')
            
            # Statistiche riassuntive
            stats_text = f"""
            Statistiche Finali:
            
            Particelle sopravvissute: {len(final_active)}/{len(initial_active)} ({survival_rate:.1f}%)
            
            Energia media iniziale: {initial_active['energy_keV'].mean():.1f} ± {initial_active['energy_keV'].std():.1f} keV
            Energia media finale: {final_active['energy_keV'].mean():.1f} ± {final_active['energy_keV'].std():.1f} keV
            
            Guadagno medio: {final_active['energy_keV'].mean() - initial_active['energy_keV'].mean():.1f} keV
            Massimo guadagno: {energy_gain.max():.1f} keV
            
            Posizione media finale: {final_active['x_um'].mean():.1f} μm
            """
            
            axes[1,1].text(0.05, 0.95, stats_text, transform=axes[1,1].transAxes, 
                           verticalalignment='top', fontsize=10, 
                           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
            axes[1,1].set_title('Statistiche Riassuntive')
            axes[1,1].axis('off')
        
        plt.tight_layout()
        plt.savefig(self.data_dir / 'analisi_energia.png', dpi=300, bbox_inches='tight')
        plt.show()
        
        # Stampa risultati
        print("\n=== ANALISI GUADAGNO ENERGETICO ===")
        print(f"Particelle iniziali: {len(initial_active)}")
        print(f"Particelle finali attive: {len(final_active)}")
        print(f"Efficienza di trasporto: {survival_rate:.1f}%")
        if len(final_active) > 0:
            print(f"Energia media iniziale: {initial_active['energy_keV'].mean():.1f} keV")
            print(f"Energia media finale: {final_active['energy_keV'].mean():.1f} keV")
            print(f"Guadagno energetico medio: {final_active['energy_keV'].mean() - initial_active['energy_keV'].mean():.1f} keV")
    
    def full_analysis(self):
        """
        Esegue l'analisi completa dei risultati PIC
        """
        print("\n=== ANALISI COMPLETA RISULTATI PIC DLA ===\n")
        
        # 1. Geometria
        print("1. Analisi geometria...")
        fields_df = self.load_fields(final=True)
        if fields_df is not None:
            self.plot_geometry(fields_df)
        
        # 2. Campi elettromagnetici
        print("\n2. Analisi campi elettromagnetici...")
        if fields_df is not None:
            self.plot_electric_field(fields_df)
        
        # 3. Traiettorie particelle
        print("\n3. Analisi traiettorie particelle...")
        self.plot_particle_trajectories(final=True)
        
        # 4. Analisi energetica
        print("\n4. Analisi guadagno energetico...")
        self.analyze_energy_gain()
        
        print("\n=== ANALISI COMPLETATA ===")
        print(f"Grafici salvati in: {self.data_dir}")

def main():
    """Funzione principale"""
    parser = argparse.ArgumentParser(description='Analizza risultati simulazione PIC DLA')
    parser.add_argument('--dir', '-d', default='.', 
                       help='Directory contenente i file CSV (default: directory corrente)')
    parser.add_argument('--analysis', '-a', choices=['geometry', 'fields', 'particles', 'energy', 'all'],
                       default='all', help='Tipo di analisi da eseguire')
    
    args = parser.parse_args()
    
    # Inizializza analizzatore
    analyzer = PICAnalyzer(args.dir)
    
    # Esegui analisi richiesta
    if args.analysis == 'geometry':
        fields_df = analyzer.load_fields(final=True)
        analyzer.plot_geometry(fields_df)
    elif args.analysis == 'fields':
        fields_df = analyzer.load_fields(final=True)
        analyzer.plot_electric_field(fields_df)
    elif args.analysis == 'particles':
        analyzer.plot_particle_trajectories(final=True)
    elif args.analysis == 'energy':
        analyzer.analyze_energy_gain()
    else:  # 'all'
        analyzer.full_analysis()

if __name__ == "__main__":
    main()