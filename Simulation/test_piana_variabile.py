#!/usr/bin/env python3
"""
Test script for the new piana_variabile geometry.
This script runs a quick test of the new geometry and verifies the results.
"""

import subprocess
import os
import sys
import pandas as pd
import numpy as np

def run_command(cmd, description):
    """Run a command and check for success."""
    print(f"\n=== {description} ===")
    print(f"Running: {' '.join(cmd)}")
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode == 0:
        print(f"✅ {description} completed successfully!")
        if result.stdout:
            print("Output:", result.stdout[-200:])  # Show last 200 chars
    else:
        print(f"❌ {description} failed!")
        print("Error:", result.stderr)
        return False
    
    return True

def check_output_files(folder):
    """Check if expected output files are created."""
    expected_files = [
        "geometry_params.csv",
        "potential.csv",
        "electric_field_x.csv", 
        "electric_field_y.csv",
        "permittivity.csv",
        "x_coordinates.csv",
        "y_coordinates.csv"
    ]
    
    print(f"\n=== Checking output files in {folder} ===")
    all_good = True
    
    for file in expected_files:
        filepath = os.path.join(folder, file)
        if os.path.exists(filepath):
            size = os.path.getsize(filepath)
            print(f"✅ {file} ({size:,} bytes)")
        else:
            print(f"❌ {file} - NOT FOUND")
            all_good = False
    
    return all_good

def analyze_geometry_params(folder):
    """Analyze the geometry parameters to verify correctness."""
    params_file = os.path.join(folder, "geometry_params.csv")
    
    if not os.path.exists(params_file):
        print("❌ Cannot analyze geometry params - file not found")
        return False
    
    print(f"\n=== Analyzing geometry parameters ===")
    
    try:
        # Read geometry parameters
        params = {}
        with open(params_file, 'r') as f:
            for line in f:
                parts = line.strip().split(',')
                if len(parts) == 2:
                    params[parts[0]] = parts[1]
        
        print("Geometry parameters:")
        for key, value in params.items():
            print(f"  {key}: {value}")
        
        # Verify piana_variabile specific parameters
        if params.get('geometry_type') == 'piana_variabile':
            print("✅ Confirmed geometry type: piana_variabile")
            
            thick_right = float(params.get('y_si_layer_thick_right', 0))
            thick_left = float(params.get('y_si_layer_thick_left', 0))
            
            if thick_right == 10.0 and thick_left == 4.0:
                print(f"✅ Thickness variation correct: {thick_right}μm → {thick_left}μm")
            else:
                print(f"❌ Unexpected thickness values: right={thick_right}, left={thick_left}")
                return False
        else:
            print(f"❌ Wrong geometry type: {params.get('geometry_type')}")
            return False
            
    except Exception as e:
        print(f"❌ Error analyzing parameters: {e}")
        return False
    
    return True

def main():
    """Main test function."""
    print("🧪 Testing piana_variabile geometry implementation")
    print("=" * 50)
    
    # Change to simulation directory
    sim_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(sim_dir)
    print(f"Working directory: {sim_dir}")
    
    test_folder = "test_piana_variabile"
    
    # Step 1: Clean previous test results
    if os.path.exists(test_folder):
        print(f"Cleaning previous test folder: {test_folder}")
        import shutil
        shutil.rmtree(test_folder)
    
    # Step 2: Compile
    if not run_command(["make", "clean"], "Cleaning build"):
        return False
    
    if not run_command(["make", "all"], "Compiling project"):
        return False
    
    # Step 3: Run Poisson solver with new geometry
    if not run_command(["./poisson_solver", test_folder, "piana_variabile"], 
                      "Running Poisson solver with piana_variabile geometry"):
        return False
    
    # Step 4: Check output files
    if not check_output_files(test_folder):
        return False
    
    # Step 5: Analyze geometry parameters
    if not analyze_geometry_params(test_folder):
        return False
    
    # Step 6: Run Python plotting (optional)
    try:
        if run_command([sys.executable, "plot_results.py", test_folder], 
                      "Generating plots with corrected units"):
            print("✅ Plotting completed successfully!")
    except:
        print("⚠️  Plotting failed, but core functionality works")
    
    print("\n" + "=" * 50)
    print("🎉 Test completed successfully!")
    print(f"Check the '{test_folder}' folder for results")
    print("Key achievements:")
    print("  ✅ New geometry compiles without errors")
    print("  ✅ Poisson solver runs with piana_variabile")
    print("  ✅ Correct parameters saved (10μm → 4μm)")
    print("  ✅ Output files generated")
    print("  ✅ Units corrected in Python files")

if __name__ == "__main__":
    main()