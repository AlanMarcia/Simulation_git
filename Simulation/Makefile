# Makefile for Poisson Solver
CXX = g++
CXXFLAGS = -std=c++17 -O3 -fopenmp -Wall -Wextra
TARGET = poisson_solver
SOURCES = poisson_solver.cpp geometry_definitions.cpp
OBJECTS = $(SOURCES:.cpp=.o)

# Default target
all: $(TARGET)

# Build the main executable
$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJECTS)

# Compile object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean build files
clean:
	rm -f $(OBJECTS) $(TARGET)

# Install dependencies (run once)
install-deps:
	sudo apt update
	sudo apt install -y build-essential g++ libomp-dev

# Run examples
run-piana:
	./$(TARGET) geometria_piana piana

run-denti-sfasati:
	./$(TARGET) geometria_Denti_sfasati_profondi denti_sfasati_profondi

run-denti-uguali:
	./$(TARGET) geometria_Denti_uguali denti_uguali

# Help
help:
	@echo "Available targets:"
	@echo "  all           - Build the program (default)"
	@echo "  clean         - Remove build files"
	@echo "  install-deps  - Install required dependencies"
	@echo "  run-piana     - Run with piana geometry"
	@echo "  run-denti-sfasati - Run with denti_sfasati_profondi geometry"
	@echo "  run-denti-uguali  - Run with denti_uguali geometry"
	@echo "  help          - Show this help"

.PHONY: all clean install-deps run-piana run-denti-sfasati run-denti-uguali help
