# Parallelizzazione della Simulazione PIC DLA

## Panoramica

La simulazione PIC (Particle-in-Cell) per l'acceleratore DLA è stata parallelizzata utilizzando OpenMP per sfruttare i processori multi-core moderni. Questa documentazione descrive le ottimizzazioni applicate e come utilizzare la versione parallelizzata.

## Ottimizzazioni Implementate

### 1. Parallelizzazione dei Loop sulla Griglia

**Funzioni ottimizzate:**
- `setupDielectricGeometry()` - Setup della geometria dielettrica
- `updateLaserField()` - Aggiornamento del campo laser
- `updateElectricField()` - Aggiornamento campo elettrico FDTD
- `updateMagneticField()` - Aggiornamento campo magnetico FDTD

**Tecnica utilizzata:**
```cpp
#pragma omp parallel for collapse(2) schedule(static)
for (int i = 0; i < Nx; i++) {
    for (int j = 0; j < Ny; j++) {
        // Calcoli sulla griglia
    }
}
```

**Benefici:**
- Sfrutta tutti i core disponibili per i calcoli su griglia
- `collapse(2)` distribuisce il lavoro su entrambi i loop
- `schedule(static)` garantisce bilanciamento del carico ottimale

### 2. Parallelizzazione delle Particelle

**Funzioni ottimizzate:**
- `updateParticles()` - Aggiornamento posizione e velocità particelle
- `initializeParticles()` - Inizializzazione particelle con RNG paralleli

**Tecnica utilizzata:**
```cpp
#pragma omp parallel for schedule(dynamic, 100)
for (int p = 0; p < n_particles; p++) {
    // Aggiornamento particella p
}
```

**Caratteristiche speciali:**
- `schedule(dynamic, 100)` per bilanciare il carico (particelle inattive)
- Generatori random separati per ogni thread nell'inizializzazione
- Gestione thread-safe degli aggiornamenti

### 3. Deposizione delle Correnti Thread-Safe

**Problema:** Race conditions nella deposizione su griglia
**Soluzione:** Array temporanei per thread + riduzione

```cpp
// Array temporanei per ogni thread
std::vector<std::vector<std::vector<double>>> thread_Jx(num_threads, ...);

#pragma omp parallel
{
    int thread_id = omp_get_thread_num();
    
    #pragma omp for schedule(static)
    for (int p = 0; p < n_particles; p++) {
        // Deposita su array del thread corrente
        thread_Jx[thread_id][i][j] += contribution;
    }
}

// Riduzione finale
#pragma omp parallel for collapse(2)
for (int i = 0; i < Nx; i++) {
    for (int j = 0; j < Ny; j++) {
        for (int t = 0; t < num_threads; t++) {
            Jx[i][j] += thread_Jx[t][i][j];
        }
    }
}
```

### 4. Parallelizzazione delle Operazioni di Riduzione

**Funzioni ottimizzate:**
- `countActiveParticles()` - Conteggio particelle attive
- `printStatus()` - Calcolo statistiche energetiche

**Tecnica utilizzata:**
```cpp
#pragma omp parallel for reduction(+:active_count) schedule(static)
for (size_t i = 0; i < particles.size(); i++) {
    if (particles[i].active) {
        active_count++;
    }
}
```

### 5. I/O Parallelo per File di Output

**Ottimizzazione:** Preparazione dati in parallelo, scrittura sequenziale

```cpp
std::vector<std::string> lines(data_size);

#pragma omp parallel for schedule(static)
for (int i = 0; i < data_size; i++) {
    // Prepara stringa per elemento i
    std::ostringstream oss;
    oss << formatted_data;
    lines[i] = oss.str();
}

// Scrittura sequenziale
for (const auto& line : lines) {
    file << line;
}
```

## Compilazione e Utilizzo

### Requisiti
- Compilatore con supporto OpenMP (GCC 4.9+, Clang 3.7+, MSVC 2015+)
- C++17 o superiore

### Compilazione

```bash
# Versione parallelizzata
make -f Makefile_Parallel

# Versione originale per confronto
make -f Makefile_Parallel original

# Versione debug
make -f Makefile_Parallel debug
```

### Esecuzione

```bash
# Usa tutti i core disponibili
./pic_dla_parallel

# Specifica il numero di thread
OMP_NUM_THREADS=4 ./pic_dla_parallel

# Con scheduling personalizzato
OMP_NUM_THREADS=8 OMP_SCHEDULE=dynamic ./pic_dla_parallel
```

### Benchmark delle Prestazioni

```bash
# Confronto automatico
make -f Makefile_Parallel benchmark

# Test con diversi numeri di thread
make -f Makefile_Parallel test_threads

# Analisi dettagliata con Python
python3 benchmark_performance.py pic_dla_parallel
```

## Variabili di Ambiente OpenMP Utili

| Variabile | Descrizione | Esempio |
|-----------|-------------|---------|
| `OMP_NUM_THREADS` | Numero di thread | `OMP_NUM_THREADS=8` |
| `OMP_SCHEDULE` | Tipo di scheduling | `OMP_SCHEDULE=dynamic,100` |
| `OMP_PROC_BIND` | Affinità dei thread | `OMP_PROC_BIND=true` |
| `OMP_PLACES` | Posizionamento thread | `OMP_PLACES=cores` |

## Prestazioni Attese

### Speedup Teorico
- **Setup geometria:** ~8x (perfettamente parallelizzabile)
- **Aggiornamento campi:** ~6-8x (bound da memoria)
- **Aggiornamento particelle:** ~4-6x (carico variabile)
- **Complessivo:** ~3-5x su sistemi quad-core, ~5-8x su sistemi octa-core

### Fattori Limitanti
1. **Memoria bandwidth** - I calcoli su griglia sono memory-bound
2. **Carico sbilanciato** - Particelle inattive riducono l'efficienza
3. **Sincronizzazione** - Barriere implicite nei parallel for
4. **I/O sequenziale** - Scrittura file non parallelizzabile

## Profiling e Ottimizzazione

### Compilazione con Profiling
```bash
make -f Makefile_Parallel profile
./pic_dla_parallel
gprof pic_dla_parallel gmon.out > profile_report.txt
```

### Monitoraggio Prestazioni
Il codice include timer automatici ogni 1000 step:
```
Step 1000 time: 45.2 ms
Step 2000 time: 44.8 ms
```

### Analisi Hotspot
Usa `perf` su Linux o `Instruments` su macOS:
```bash
perf record ./pic_dla_parallel
perf report
```

## Consigli per l'Ottimizzazione

### 1. Numero Ottimale di Thread
- **CPU-bound tasks:** Numero di core fisici
- **Memory-bound tasks:** Spesso meno del numero di core
- **Test empirico:** Usa `benchmark_performance.py`

### 2. Scheduling Ottimale
- **Calcoli uniformi:** `schedule(static)`
- **Carico variabile:** `schedule(dynamic, chunk_size)`
- **Particelle:** `schedule(dynamic, 100)` funziona bene

### 3. Ottimizzazioni di Memoria
- **Locality:** Mantieni dati correlati vicini
- **Alignment:** Usa `alignas(64)` per strutture grandi
- **NUMA:** Su sistemi multi-socket, considera `numactl`

### 4. Scalabilità
- **Strong scaling:** Fisso problem size, aumenta thread
- **Weak scaling:** Aumenta sia problem size che thread
- **Sweet spot:** Trova il numero ottimale di thread per il tuo hardware

## Debugging della Parallelizzazione

### Race Conditions
```bash
# Compila con thread sanitizer
g++ -fsanitize=thread -g ...

# Usa helgrind (Valgrind)
valgrind --tool=helgrind ./pic_dla_parallel
```

### Deadlock
```bash
# Debug con GDB
gdb ./pic_dla_parallel
(gdb) set environment OMP_NUM_THREADS=2
(gdb) run
# Su deadlock: Ctrl+C, poi 'thread apply all bt'
```

### Performance Issues
```bash
# Controlla affinità dei thread
OMP_DISPLAY_AFFINITY=true ./pic_dla_parallel

# Usa Intel VTune (se disponibile)
amplxe-cl -collect hotspots ./pic_dla_parallel
```

## Limitazioni e Work-Around

### 1. Memory Bandwidth
**Problema:** Calcoli limitati dalla velocità di memoria
**Soluzione:** 
- Blocking dei dati
- Uso di registri e cache L1
- Prefetching esplicito

### 2. Load Imbalancing
**Problema:** Thread con carichi di lavoro diversi
**Soluzione:**
- Dynamic scheduling
- Work stealing
- Particelle raggruppate per attività

### 3. NUMA Effects
**Problema:** Latenza memoria non uniforme
**Soluzione:**
```bash
numactl --interleave=all ./pic_dla_parallel
# oppure
numactl --cpunodebind=0 --membind=0 ./pic_dla_parallel
```

## Estensioni Future

### 1. GPU Parallelization (CUDA/OpenCL)
- Aggiornamento particelle su GPU
- FDTD solver su GPU
- Trasferimenti memoria ottimizzati

### 2. Distributed Computing (MPI)
- Domain decomposition
- Particle migration tra nodi
- Load balancing dinamico

### 3. Vectorization (AVX/SSE)
- Istruzioni SIMD per particelle
- Calcoli su griglia vettorizzati
- Compilazione con `-march=native -O3`

### 4. Algorithmic Improvements
- Tree codes per long-range forces
- Adaptive mesh refinement
- Higher-order particle shapes

## Riferimenti

1. [OpenMP 5.0 Specification](https://www.openmp.org/spec-html/5.0/openmp.html)
2. [Intel OpenMP Runtime Library](https://github.com/llvm-mirror/openmp)
3. [PIC Simulation Techniques](https://doi.org/10.1016/0021-9991(73)90157-5)
4. [Parallel PIC Algorithms](https://doi.org/10.1006/jcph.1995.1181)