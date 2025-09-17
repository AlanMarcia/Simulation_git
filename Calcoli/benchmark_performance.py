#!/usr/bin/env python3
"""
Script per analizzare le prestazioni della simulazione PIC parallelizzata
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import time
import subprocess
import os
import sys

def run_benchmark(executable, num_threads_list, output_dir="benchmark_results"):
    """
    Esegue benchmark con diversi numeri di thread
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    results = []
    
    print(f"=== BENCHMARK DI {executable} ===")
    
    for num_threads in num_threads_list:
        print(f"Testando con {num_threads} thread...")
        
        # Imposta variabile di ambiente
        env = os.environ.copy()
        env['OMP_NUM_THREADS'] = str(num_threads)
        
        # Misura tempo di esecuzione
        start_time = time.time()
        
        try:
            # Esegue il programma
            result = subprocess.run(
                f'./{executable}', 
                shell=True, 
                env=env, 
                capture_output=True, 
                text=True,
                timeout=300  # timeout di 5 minuti
            )
            
            end_time = time.time()
            execution_time = end_time - start_time
            
            if result.returncode == 0:
                print(f"  Completato in {execution_time:.2f} secondi")
                results.append({
                    'threads': num_threads,
                    'time': execution_time,
                    'speedup': 0,  # Calcolato dopo
                    'efficiency': 0  # Calcolato dopo
                })
            else:
                print(f"  Errore nell'esecuzione: {result.stderr}")
                
        except subprocess.TimeoutExpired:
            print(f"  Timeout raggiunto per {num_threads} thread")
        except Exception as e:
            print(f"  Errore: {e}")
    
    # Calcola speedup ed efficienza
    if results:
        base_time = results[0]['time']  # Tempo con 1 thread
        for result in results:
            result['speedup'] = base_time / result['time']
            result['efficiency'] = result['speedup'] / result['threads']
    
    return results

def plot_performance_analysis(results, output_dir="benchmark_results"):
    """
    Crea grafici di analisi delle prestazioni
    """
    if not results:
        print("Nessun risultato da plottare")
        return
    
    df = pd.DataFrame(results)
    
    # Setup per i grafici
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 10))
    fig.suptitle('Analisi Prestazioni PIC DLA Simulator', fontsize=16)
    
    # 1. Tempo di esecuzione vs Numero di thread
    ax1.plot(df['threads'], df['time'], 'bo-', linewidth=2, markersize=8)
    ax1.set_xlabel('Numero di Thread')
    ax1.set_ylabel('Tempo di Esecuzione (s)')
    ax1.set_title('Tempo di Esecuzione')
    ax1.grid(True, alpha=0.3)
    ax1.set_yscale('log')
    
    # 2. Speedup vs Numero di thread
    ax2.plot(df['threads'], df['speedup'], 'ro-', linewidth=2, markersize=8, label='Speedup Effettivo')
    ax2.plot(df['threads'], df['threads'], 'k--', alpha=0.5, label='Speedup Ideale')
    ax2.set_xlabel('Numero di Thread')
    ax2.set_ylabel('Speedup')
    ax2.set_title('Speedup')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # 3. Efficienza vs Numero di thread
    ax3.plot(df['threads'], df['efficiency'] * 100, 'go-', linewidth=2, markersize=8)
    ax3.axhline(y=100, color='k', linestyle='--', alpha=0.5, label='Efficienza Ideale')
    ax3.set_xlabel('Numero di Thread')
    ax3.set_ylabel('Efficienza (%)')
    ax3.set_title('Efficienza di Parallelizzazione')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    ax3.set_ylim(0, 110)
    
    # 4. Tabella riassuntiva
    ax4.axis('tight')
    ax4.axis('off')
    
    table_data = []
    for _, row in df.iterrows():
        table_data.append([
            f"{int(row['threads'])}",
            f"{row['time']:.2f}s",
            f"{row['speedup']:.2f}x",
            f"{row['efficiency']*100:.1f}%"
        ])
    
    table = ax4.table(
        cellText=table_data,
        colLabels=['Thread', 'Tempo', 'Speedup', 'Efficienza'],
        loc='center',
        cellLoc='center'
    )
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 2)
    ax4.set_title('Risultati Numerici')
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/performance_analysis.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Salva risultati in CSV
    df.to_csv(f'{output_dir}/benchmark_results.csv', index=False)
    print(f"Risultati salvati in {output_dir}/benchmark_results.csv")

def analyze_parallel_efficiency(results):
    """
    Analizza l'efficienza della parallelizzazione
    """
    if not results or len(results) < 2:
        print("Dati insufficienti per l'analisi")
        return
    
    print("\n=== ANALISI EFFICIENZA PARALLELIZZAZIONE ===")
    
    df = pd.DataFrame(results)
    max_speedup = df['speedup'].max()
    max_speedup_threads = df.loc[df['speedup'].idxmax(), 'threads']
    
    print(f"Speedup massimo: {max_speedup:.2f}x con {max_speedup_threads} thread")
    
    # Trova il punto di efficienza ottimale (>80%)
    efficient_results = df[df['efficiency'] > 0.8]
    if not efficient_results.empty:
        optimal_threads = efficient_results['threads'].max()
        optimal_speedup = efficient_results.loc[efficient_results['threads'].idxmax(), 'speedup']
        print(f"Configurazione ottimale: {optimal_threads} thread (speedup: {optimal_speedup:.2f}x)")
    else:
        print("Nessuna configurazione con efficienza >80% trovata")
    
    # Analisi della scalabilità
    if len(results) >= 3:
        print("\nScalabilità:")
        for i in range(1, len(results)):
            prev_threads = results[i-1]['threads']
            curr_threads = results[i]['threads']
            prev_time = results[i-1]['time']
            curr_time = results[i]['time']
            
            theoretical_speedup = curr_threads / prev_threads
            actual_speedup = prev_time / curr_time
            scalability = actual_speedup / theoretical_speedup * 100
            
            print(f"  {prev_threads} → {curr_threads} thread: {scalability:.1f}% della scalabilità ideale")

def main():
    """
    Funzione principale per il benchmark
    """
    if len(sys.argv) > 1:
        executable = sys.argv[1]
    else:
        executable = "pic_dla_parallel"
    
    # Verifica che l'eseguibile esista
    if not os.path.exists(executable):
        print(f"Errore: {executable} non trovato!")
        print("Compila prima il programma con: make")
        return
    
    # Lista dei numeri di thread da testare
    import multiprocessing
    max_cores = multiprocessing.cpu_count()
    thread_list = [1, 2, 4]
    
    # Aggiungi più configurazioni se ci sono abbastanza core
    if max_cores >= 8:
        thread_list.extend([8])
    if max_cores >= 16:
        thread_list.extend([16])
    
    print(f"Sistema con {max_cores} core logici")
    print(f"Testando con thread: {thread_list}")
    
    # Esegue il benchmark
    results = run_benchmark(executable, thread_list)
    
    if results:
        # Analizza e visualizza i risultati
        analyze_parallel_efficiency(results)
        plot_performance_analysis(results)
        
        print(f"\n=== RACCOMANDAZIONI ===")
        df = pd.DataFrame(results)
        best_config = df.loc[df['efficiency'].idxmax()]
        print(f"Per le migliori prestazioni, usa: OMP_NUM_THREADS={int(best_config['threads'])}")
        print(f"Questo darà un speedup di {best_config['speedup']:.2f}x con efficienza {best_config['efficiency']*100:.1f}%")
        
    else:
        print("Nessun risultato ottenuto dal benchmark")

if __name__ == "__main__":
    main()