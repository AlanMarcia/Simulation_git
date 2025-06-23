#!/bin/bash

# Build and Run Script for Poisson Solver
# Usage: ./build_and_run.sh [geometry_type] [output_folder]

set -e  # Exit on any error

echo "=== Poisson Solver Build and Run Script ==="

# Check if we're in the right directory
if [ ! -f "poisson_solver.cpp" ]; then
    echo "Error: poisson_solver.cpp not found. Please run this script from the simulation directory."
    exit 1
fi

# Install dependencies if needed
if ! command -v g++ &> /dev/null; then
    echo "Installing build dependencies..."
    sudo apt update
    sudo apt install -y build-essential g++ libomp-dev
fi

# Build the program
echo "Building poisson_solver..."
make clean
make

if [ $? -eq 0 ]; then
    echo "✓ Build successful!"
else
    echo "✗ Build failed!"
    exit 1
fi

# Default parameters
GEOMETRY_TYPE=${1:-"piana"}
OUTPUT_FOLDER=${2:-"test_output"}

echo "Running simulation with:"
echo "  Geometry type: $GEOMETRY_TYPE"
echo "  Output folder: $OUTPUT_FOLDER"
echo ""

# Run the simulation
./poisson_solver "$OUTPUT_FOLDER" "$GEOMETRY_TYPE"

if [ $? -eq 0 ]; then
    echo ""
    echo "✓ Simulation completed successfully!"
    echo "Results saved in: $OUTPUT_FOLDER"
    echo ""
    echo "Generated files:"
    ls -la "$OUTPUT_FOLDER"/ 2>/dev/null || echo "Output directory not found"
else
    echo "✗ Simulation failed!"
    exit 1
fi
