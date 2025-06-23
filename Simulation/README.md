# Poisson Solver - Quick Start Guide

## Prerequisites

In Ubuntu bash (WSL), install the required dependencies:

```bash
sudo apt update
sudo apt install build-essential g++ libomp-dev
```

## Quick Start

### Option 1: Using the build script (Recommended)

1. Make the script executable:
```bash
chmod +x build_and_run.sh
```

2. Run with default settings (piana geometry):
```bash
./build_and_run.sh
```

3. Run with specific geometry:
```bash
./build_and_run.sh piana my_output_folder
./build_and_run.sh denti_sfasati_profondi results_deep_teeth
./build_and_run.sh denti_uguali results_equal_teeth
```

### Option 2: Using Makefile

1. Install dependencies (run once):
```bash
make install-deps
```

2. Build the program:
```bash
make
```

3. Run examples:
```bash
make run-piana
make run-denti-sfasati
make run-denti-uguali
```

### Option 3: Manual compilation

```bash
# Compile
g++ -std=c++17 -O3 -fopenmp -Wall -Wextra -o poisson_solver poisson_solver.cpp geometry_definitions.cpp

# Run
./poisson_solver output_folder_name geometry_type
```

## Available Geometry Types

- `piana` - Flat geometry (default)
- `denti_sfasati_profondi` - Deep staggered teeth geometry
- `denti_uguali` - Equal teeth geometry

## Output

The program generates CSV files in the specified output folder:
- `potential.csv` - Electric potential field
- `electric_field_x.csv` - X-component of electric field
- `electric_field_y.csv` - Y-component of electric field  
- `permittivity.csv` - Material permittivity map
- `x_coordinates.csv` - X-axis coordinates
- `y_coordinates.csv` - Y-axis coordinates
- `geometry_params.csv` - Geometry parameters used

## Visualization

Use the provided Python scripts to visualize results:
```bash
python3 plot_results.py
```

## Example Commands

```bash
# Simple run with default settings
./build_and_run.sh

# Run with specific geometry and output folder
./build_and_run.sh denti_sfasati_profondi Silicio/geometria_Denti_sfasati_test

# Clean and rebuild
make clean && make

# Get help
make help
```
