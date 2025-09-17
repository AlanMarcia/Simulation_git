# Makefile per simulazioni PIC e calcoli di accelerazione
# Compilatore e flags
CXX = g++
CXXFLAGS = -std=c++17 -O3 -Wall -Wextra -fopenmp
LDFLAGS = -lm -fopenmp

# Target executables
TARGETS = PIC_DLA_Accelerator Calcoli PIC

# Default target
all: $(TARGETS)

# Compilazione PIC DLA Accelerator
PIC_DLA_Accelerator: PIC_DLA_Accelerator.cpp
	@echo "Compilando PIC DLA Accelerator..."
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)
	@echo "✓ PIC_DLA_Accelerator compilato con successo"

# Compilazione Calcoli
Calcoli: Calcoli.cpp
	@echo "Compilando Calcoli..."
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)
	@echo "✓ Calcoli compilato con successo"

# Compilazione PIC (se esiste)
PIC: PIC
	@echo "Copiando eseguibile PIC esistente..."
	@cp PIC PIC_exec
	@echo "✓ PIC copiato come PIC_exec"

# Target per compilare solo PIC
pic: PIC_DLA_Accelerator
	@echo "PIC DLA Accelerator pronto per l'esecuzione"

# Target per compilare solo Calcoli
calcoli: Calcoli
	@echo "Calcoli pronto per l'esecuzione"

# Target per eseguire PIC
run_pic: PIC_DLA_Accelerator
	@echo "Eseguendo simulazione PIC DLA..."
	./PIC_DLA_Accelerator

# Target per eseguire Calcoli
run_calcoli: Calcoli
	@echo "Eseguendo calcoli di accelerazione..."
	./Calcoli

# Target per pulire i file compilati
clean:
	@echo "Rimuovendo file compilati..."
	rm -f $(TARGETS) PIC_exec
	rm -f *.csv *.dat *.out
	@echo "✓ Pulizia completata"

# Target per pulire solo gli eseguibili
clean_exe:
	@echo "Rimuovendo solo gli eseguibili..."
	rm -f $(TARGETS) PIC_exec
	@echo "✓ Eseguibili rimossi"

# Target per pulire solo i risultati
clean_results:
	@echo "Rimuovendo file di risultati..."
	rm -f *.csv *.dat *.out
	@echo "✓ File di risultati rimossi"

# Target di debug (con simboli di debug)
debug: CXXFLAGS += -g -DDEBUG
debug: clean_exe $(TARGETS)
	@echo "✓ Versione debug compilata"

# Target per compilazione veloce (senza ottimizzazioni)
fast: CXXFLAGS = -std=c++17 -Wall -fopenmp
fast: clean_exe $(TARGETS)
	@echo "✓ Compilazione veloce completata"

# Target per informazioni
info:
	@echo "=== Makefile per Simulazioni PIC ==="
	@echo "Targets disponibili:"
	@echo "  all          - Compila tutti i programmi"
	@echo "  pic          - Compila solo PIC_DLA_Accelerator"
	@echo "  calcoli      - Compila solo Calcoli"
	@echo "  run_pic      - Compila ed esegue PIC"
	@echo "  run_calcoli  - Compila ed esegue Calcoli"
	@echo "  clean        - Rimuove tutti i file generati"
	@echo "  clean_exe    - Rimuove solo gli eseguibili"
	@echo "  clean_results- Rimuove solo i file di risultati"
	@echo "  debug        - Compila versione debug"
	@echo "  fast         - Compilazione veloce senza ottimizzazioni"
	@echo "  info         - Mostra questo messaggio"
	@echo ""
	@echo "Compilatore: $(CXX)"
	@echo "Flags: $(CXXFLAGS)"
	@echo "Link flags: $(LDFLAGS)"

# Help alias
help: info

# Prevent make from deleting intermediate files
.PRECIOUS: %.o

# Phony targets (non sono file)
.PHONY: all pic calcoli run_pic run_calcoli clean clean_exe clean_results debug fast info help