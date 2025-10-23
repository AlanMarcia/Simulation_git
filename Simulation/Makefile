# Makefile for Poisson Solver
CXX = g++
CXXFLAGS = -std=c++17 -O3 -fopenmp -Wall -Wextra
TARGET_POISSON = poisson_solver
TARGET_PROTON = proton_simulator
SOURCES_POISSON = poisson_solver.cpp geometry_definitions.cpp
SOURCES_PROTON = proton_simulator.cpp geometry_definitions.cpp
OBJECTS_POISSON = $(SOURCES_POISSON:.cpp=.o)
OBJECTS_PROTON = $(SOURCES_PROTON:.cpp=.o)

# Python command (adjust if needed)
PYTHON = python3

# Geometry types
GEOMETRIES = piana piana_rastremata piana_variabile denti_sfasati_profondi denti_uguali

# Default target
all: $(TARGET_POISSON) $(TARGET_PROTON)

# Build the main executables
$(TARGET_POISSON): $(OBJECTS_POISSON)
	$(CXX) $(CXXFLAGS) -o $(TARGET_POISSON) $(OBJECTS_POISSON)

$(TARGET_PROTON): $(OBJECTS_PROTON)
	$(CXX) $(CXXFLAGS) -o $(TARGET_PROTON) $(OBJECTS_PROTON)

# Compile object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean build files
clean:
	rm -f $(OBJECTS_POISSON) $(OBJECTS_PROTON) $(TARGET_POISSON) $(TARGET_PROTON)

# Install dependencies (run once)
install-deps:
	sudo apt update
	sudo apt install -y build-essential g++ libomp-dev python3 python3-pip
	pip3 install matplotlib numpy pandas

# ==== INDIVIDUAL GEOMETRY RUNS ====

# Run piana geometry (all steps)
run-piana: $(TARGET_POISSON) $(TARGET_PROTON)
	@echo "=== Running Piana Geometry ==="
	mkdir -p geometria_piana
	./$(TARGET_POISSON) geometria_piana piana
	./$(TARGET_PROTON) geometria_piana
	$(PYTHON) plot_results.py geometria_piana/
	$(PYTHON) plot_trajectories.py geometria_piana/
	@echo "=== Piana Geometry Complete ==="

# Run piana rastremata geometry (all steps)  
run-piana-rastremata: $(TARGET_POISSON) $(TARGET_PROTON)
	@echo "=== Running Piana Rastremata Geometry ==="
	mkdir -p geometria_piana_rastremata
	./$(TARGET_POISSON) geometria_piana_rastremata piana_rastremata
	./$(TARGET_PROTON) geometria_piana_rastremata
	$(PYTHON) plot_results.py geometria_piana_rastremata/
	$(PYTHON) plot_trajectories.py geometria_piana_rastremata/
	@echo "=== Piana Rastremata Geometry Complete ==="

# Run denti sfasati profondi geometry (all steps)
run-denti-sfasati: $(TARGET_POISSON) $(TARGET_PROTON)
	@echo "=== Running Denti Sfasati Profondi Geometry ==="
	mkdir -p geometria_Denti_sfasati_profondi
	./$(TARGET_POISSON) geometria_Denti_sfasati_profondi denti_sfasati_profondi
	./$(TARGET_PROTON) geometria_Denti_sfasati_profondi
	$(PYTHON) plot_results.py geometria_Denti_sfasati_profondi/
	$(PYTHON) plot_trajectories.py geometria_Denti_sfasati_profondi/
	@echo "=== Denti Sfasati Profondi Geometry Complete ==="

# Run piana variabile geometry (all steps)
run-piana-variabile: $(TARGET_POISSON) $(TARGET_PROTON)
	@echo "=== Running Piana Variabile Geometry ==="
	mkdir -p geometria_piana_variabile
	./$(TARGET_POISSON) geometria_piana_variabile piana_variabile
	./$(TARGET_PROTON) geometria_piana_variabile
	$(PYTHON) plot_results.py geometria_piana_variabile/
	$(PYTHON) plot_trajectories.py geometria_piana_variabile/
	@echo "=== Piana Variabile Geometry Complete ==="

# Run denti uguali geometry (all steps)
run-denti-uguali: $(TARGET_POISSON) $(TARGET_PROTON)
	@echo "=== Running Denti Uguali Geometry ==="
	mkdir -p geometria_Denti_uguali
	./$(TARGET_POISSON) geometria_Denti_uguali denti_uguali
	./$(TARGET_PROTON) geometria_Denti_uguali
	$(PYTHON) plot_results.py geometria_Denti_uguali/
	$(PYTHON) plot_trajectories.py geometria_Denti_uguali/
	@echo "=== Denti Uguali Geometry Complete ==="

# ==== BATCH OPERATIONS ====

# Run ALL geometries sequentially
run-all-geometries: run-piana run-piana-rastremata run-piana-variabile run-denti-sfasati run-denti-uguali
	@echo "=== ALL GEOMETRIES COMPLETED ==="

# Run only Poisson solver for all geometries
run-all-poisson: $(TARGET_POISSON)
	@echo "=== Running Poisson Solver for All Geometries ==="
	mkdir -p geometria_piana && ./$(TARGET_POISSON) geometria_piana piana
	mkdir -p geometria_piana_rastremata && ./$(TARGET_POISSON) geometria_piana_rastremata piana_rastremata
	mkdir -p geometria_piana_variabile && ./$(TARGET_POISSON) geometria_piana_variabile piana_variabile
	mkdir -p geometria_Denti_sfasati_profondi && ./$(TARGET_POISSON) geometria_Denti_sfasati_profondi denti_sfasati_profondi
	mkdir -p geometria_Denti_uguali && ./$(TARGET_POISSON) geometria_Denti_uguali denti_uguali
	@echo "=== All Poisson Solvers Complete ==="

# Run only Proton simulator for all geometries
run-all-proton: $(TARGET_PROTON)
	@echo "=== Running Proton Simulator for All Geometries ==="
	./$(TARGET_PROTON) geometria_piana
	./$(TARGET_PROTON) geometria_piana_rastremata
	./$(TARGET_PROTON) geometria_piana_variabile
	./$(TARGET_PROTON) geometria_Denti_sfasati_profondi
	./$(TARGET_PROTON) geometria_Denti_uguali
	@echo "=== All Proton Simulators Complete ==="

# Run only Python plotting for all geometries
run-all-plots:
	@echo "=== Running Python Plots for All Geometries ==="
	$(PYTHON) plot_results.py geometria_piana/
	$(PYTHON) plot_trajectories.py geometria_piana/
	$(PYTHON) plot_results.py geometria_piana_rastremata/
	$(PYTHON) plot_trajectories.py geometria_piana_rastremata/
	$(PYTHON) plot_results.py geometria_piana_variabile/
	$(PYTHON) plot_trajectories.py geometria_piana_variabile/
	$(PYTHON) plot_results.py geometria_Denti_sfasati_profondi/
	$(PYTHON) plot_trajectories.py geometria_Denti_sfasati_profondi/
	$(PYTHON) plot_results.py geometria_Denti_uguali/
	$(PYTHON) plot_trajectories.py geometria_Denti_uguali/
	@echo "=== All Plots Complete ==="

# Help
help:
	@echo "Available targets:"
	@echo ""
	@echo "=== BUILD TARGETS ==="
	@echo "  all                    - Build all programs (default)"
	@echo "  clean                  - Remove build files"
	@echo "  install-deps           - Install required dependencies"
	@echo ""
	@echo "=== INDIVIDUAL GEOMETRY RUNS (C++ + Python) ==="
	@echo "  run-piana             - Run complete pipeline for piana geometry"
	@echo "  run-piana-rastremata  - Run complete pipeline for piana rastremata geometry"
	@echo "  run-piana-variabile   - Run complete pipeline for piana variabile geometry"
	@echo "  run-denti-sfasati     - Run complete pipeline for denti sfasati profondi geometry"
	@echo "  run-denti-uguali      - Run complete pipeline for denti uguali geometry"
	@echo ""
	@echo "=== BATCH OPERATIONS ==="
	@echo "  run-all-geometries    - Run ALL geometries (Poisson + Proton + Plots)"
	@echo "  run-all-poisson       - Run Poisson solver for all geometries"
	@echo "  run-all-proton        - Run Proton simulator for all geometries"
	@echo "  run-all-plots         - Run Python plots for all geometries"
	@echo ""
	@echo "=== PIPELINE FOR EACH GEOMETRY ==="
	@echo "  1. Poisson Solver     - Calculates electric field"
	@echo "  2. Proton Simulator   - Simulates particle trajectories"
	@echo "  3. Python Plots       - Generates visualization plots"
	@echo ""
	@echo "  help                  - Show this help"

.PHONY: all clean install-deps run-piana run-piana-rastremata run-piana-variabile run-denti-sfasati run-denti-uguali \
        run-all-geometries run-all-poisson run-all-proton run-all-plots help
