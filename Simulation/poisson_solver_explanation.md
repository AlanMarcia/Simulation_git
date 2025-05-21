# Code Explanation: 2D Poisson Solver for Electrostatic Simulation

This C++ program simulates the electrostatic potential and electric field in a 2D domain with defined silicon structures and a vacuum gap. It uses the Finite Difference Method (FDM) to discretize Poisson's equation and the Successive Over-Relaxation (SOR) method to solve the resulting system of linear equations. The geometry includes two silicon layers with a base and periodic teeth, separated by a vacuum gap.

## 1. Header Files and Preprocessor Definitions

*   **Standard Libraries**:
    *   `<iostream>`: For console input/output (e.g., `std::cout`, `std::cerr`).
    *   `<vector>`: For using dynamic arrays (`std::vector`).
    *   `<cmath>`: For mathematical functions (e.g., `std::abs`, `fmod`, `sqrt`).
    *   `<fstream>`: For file input/output (e.g., `std::ofstream` for writing CSV files).
    *   `<iomanip>`: For output formatting (e.g., `std::fixed`, `std::setprecision`, `std::scientific`).
    *   `<string>`: For using `std::string`.
    *   `<algorithm>`: For `std::min`, `std::max`.
*   **OpenMP**:
    *   `<omp.h>`: For parallel processing capabilities to speed up SOR iterations.
*   **Platform-Specific Directory Creation**:
    *   `#if defined(_WIN32)`: Conditional compilation for Windows.
        *   `<direct.h>`: For `_mkdir`.
        *   `MKDIR(path)` macro defined as `_mkdir(path)`.
        *   `STAT_STRUCT` and `STAT_FUNC` for checking directory existence.
    *   `#else`: Conditional compilation for non-Windows (e.g., Linux, macOS).
        *   `<sys/stat.h>` and `<sys/types.h>`: For `mkdir` and `stat`.
        *   `MKDIR(path)` macro defined as `mkdir(path, 0755)` (creates directory with read/write/execute permissions for owner, and read/execute for group/others).
        *   `STAT_STRUCT` and `STAT_FUNC` for checking directory existence.

## 2. Helper Functions

These functions encapsulate common tasks, making the `main` function cleaner.

### `saveToCSV(const std::vector<std::vector<double>>& data, const std::string& filename)`

*   Saves a 2D vector (matrix) of `double` values to a CSV file.
*   Iterates through columns (y-dimension) first, then rows (x-dimension) to match typical plotting software orientation.
*   Values in each row are separated by commas.
*   Prints a success or error message to the console.

### `saveCoordinatesToCSV(const std::vector<double>& coords, const std::string& filename)`

*   Saves a 1D vector of `double` values (typically x or y coordinates) to a CSV file.
*   Each value is written on a new line (single column).
*   Prints a success or error message to the console.

### `saveGeometryParamsToCSV(...)`

*   Saves various geometry and simulation parameters to a CSV file.
*   Parameters include:
    *   `h_val`: Grid spacing.
    *   `x_fs`: X-coordinate of free space before the structure.
    *   `x_sl`: Length of the silicon structure.
    *   `y_sibh`: Height of the silicon base.
    *   `y_th`: Height of the silicon teeth.
    *   `y_vgt`: Thickness of the vacuum gap.
    *   `x_tw`: Width of the silicon teeth.
    *   `x_ts`: Spacing between silicon teeth.
    *   `H_tot`: Total height of the simulation domain.
*   Uses `std::fixed` and `std::setprecision(10)` for consistent floating-point output.
*   Prints a success message to the console.

## 3. `main()` Function

This is the entry point and core of the simulation.

### 3.1. Parameters

*   **Grid and Dimensions**:
    *   `h`: Grid spacing (0.5 µm).
    *   `L_total`: Total length of the simulation domain (320 µm).
    *   `H_total`: Total height of the simulation domain (30 µm). This is calculated based on silicon and vacuum layer thicknesses.
    *   `x_free_space`: Length of the free space region before the silicon structure (10 µm).
    *   `x_structure_len`: Length of the silicon structure (300 µm).
*   **Silicon Layer Geometry (Toothed Structure)**:
    *   `y_si_base_height`: Height of the base of the silicon layers (5 µm).
    *   `y_teeth_height`: Height of the teeth on the silicon layers (5 µm).
    *   `x_teeth_width`: Width of each tooth (10 µm).
    *   `x_teeth_spacing`: Spacing between teeth (10 µm).
    *   `tooth_period`: Calculated as `x_teeth_width + x_teeth_spacing`.
    *   `y_vacuum_gap_thick`: Thickness of the vacuum gap between the two silicon layers (10 µm).
*   **SOR Parameters**:
    *   `omega`: Relaxation factor for SOR (1.8). Typically between 1 and 2.
    *   `tolerance`: Convergence criterion for SOR (1e-4). Iteration stops when the maximum potential change in an iteration falls below this value.
    *   `max_iterations`: Maximum number of SOR iterations (500,000) to prevent infinite loops if convergence is not met.
*   **Material Properties**:
    *   `eps_si`: Relative permittivity of silicon (11.7).
    *   `eps_vac`: Relative permittivity of vacuum (1.0).
*   **Boundary Conditions**:
    *   `V_left`: Potential applied to the left side of the silicon structures (0.0 Volts).
    *   `V_right`: Potential applied to the right side of the silicon structures (-1000.0 Volts).
*   **Output Folder**:
    *   `output_folder`: Name of the directory to save results ("geometria_Denti_uguali").

### 3.2. Output Folder Creation

*   Uses `STAT_FUNC` and `MKDIR` (platform-specific macros) to check if the `output_folder` exists.
*   If it doesn't exist, it attempts to create it.
*   If it exists but is not a directory, an error is reported.
*   Error handling is included for directory creation failure.

### 3.3. Save Geometry Parameters

*   Calls `saveGeometryParamsToCSV` to store the defined geometry and key simulation parameters.

### 3.4. Grid Setup

*   `Nx`, `Ny`: Number of grid points in x and y directions, calculated based on total dimensions and grid spacing `h`.
*   `x_coords`, `y_coords`: 1D vectors storing the physical coordinates of grid points.
*   `V`: 2D vector (`Nx` x `Ny`) storing the electric potential at each grid point, initialized to 0.0.
*   `eps_r`: 2D vector (`Nx` x `Ny`) storing the relative permittivity at each grid point, initialized to `eps_vac`.
*   `fixed_potential_mask`: 2D boolean vector (`Nx` x `Ny`) marking grid points where the potential is fixed (Dirichlet boundary conditions), initialized to `false`.

### 3.5. Define Material Regions (Permittivity Map `eps_r`)

This section defines the geometry of the silicon structures by assigning `eps_si` to the appropriate grid points in the `eps_r` matrix.
*   **Index Calculations**:
    *   `idx_x_struct_start`, `idx_x_struct_end`: Integer indices corresponding to the start and end x-positions of the silicon structure.
    *   `idx_y_bot_si_base_end`: Y-index for the top of the bottom silicon base.
    *   `idx_y_bot_si_teeth_end`: Y-index for the top of the bottom silicon teeth.
    *   `idx_y_vac_end`: Y-index for the end of the vacuum gap (and start of the top silicon teeth).
    *   `idx_y_top_si_teeth_start`: Y-index for the start of the top silicon teeth (bottom edge of the top Si structure).
    *   `idx_y_top_si_base_start`: Y-index for the start of the top silicon base (bottom edge of the top Si base).
*   **Looping through Grid**: Iterates through all grid points (`i` for x, `j` for y).
    *   Defaults `eps_r[i][j]` to `eps_vac`.
    *   **Inside Structure Region (`idx_x_struct_start` to `idx_x_struct_end`)**:
        *   **Bottom Silicon Layer**:
            *   **Base**: If `j` is within the base height (`0` to `idx_y_bot_si_base_end`), `eps_r[i][j] = eps_si`.
            *   **Teeth**: If `j` is within the teeth height (above the base, up to `idx_y_bot_si_teeth_end`):
                *   `x_coord_in_structure`: Calculates the x-position relative to the start of the structure.
                *   `pos_in_period`: Uses `fmod(x_coord_in_structure, tooth_period)` to find the position within a single tooth-gap period.
                *   If `pos_in_period < x_teeth_width`, the point is within a tooth, so `eps_r[i][j] = eps_si`.
        *   **Top Silicon Layer**: (Mirrored structure, teeth pointing downwards)
            *   **Base**: If `j` is within the top base region (`idx_y_top_si_base_start` to `Ny-1`), `eps_r[i][j] = eps_si`.
            *   **Teeth**: If `j` is within the top teeth region (below the top base, from `idx_y_top_si_teeth_start` to `idx_y_top_si_base_start - 1`):
                *   Similar `x_coord_in_structure` and `pos_in_period` logic as bottom teeth.
                *   If `pos_in_period < x_teeth_width`, the point is within a tooth, so `eps_r[i][j] = eps_si`.
*   **Console Output**: Prints grid size and key y-indices for verification.

### 3.6. Boundary Conditions (Potential `V` and `fixed_potential_mask`)

*   Iterates along the y-axis at the start (`idx_x_struct_start`) and end (`idx_x_struct_end`) x-indices of the structure.
*   If a grid point `(idx_x_struct_start, j)` or `(idx_x_struct_end, j)` is identified as silicon (by checking `eps_r`), its potential is set:
    *   `V[idx_x_struct_start][j] = V_left`
    *   `V[idx_x_struct_end][j] = V_right`
*   The corresponding `fixed_potential_mask` for these points is set to `true`.

### 3.7. SOR Iteration

This is the core solver loop.
*   **Outer Loop**: Iterates from `0` to `max_iterations - 1`.
    *   `max_diff_iter`: Stores the maximum absolute change in potential during the current iteration. Reset to 0.0 at the start of each iteration.
*   **Red-Black SOR**: The grid points are updated in two phases (red and black points, like a checkerboard pattern) to allow for effective parallelization. This ordering ensures that when updating a point, its neighbors used in the stencil have values from the previous half-iteration or the current half-iteration (for already updated points of the same color in sequential processing, or previous full iteration for points of the other color).
    *   **Phase 1: Update "Red" Points (`(i+j) % 2 == 0`)**:
        *   `#pragma omp parallel for reduction(max:max_diff_iter)`: Parallelizes the loops over `i` and `j`. The `reduction(max:max_diff_iter)` clause ensures that `max_diff_iter` is correctly updated by finding the maximum value across all threads.
        *   Iterates `i` from `1` to `Nx-2` and `j` from `1` to `Ny-2` (interior points).
        *   Skips points with fixed potentials (`fixed_potential_mask[i][j]`).
        *   **Interface Permittivity**: Calculates effective permittivities at interfaces between the current cell `(i,j)` and its neighbors:
            *   `eps_e = (eps_r[i][j] + eps_r[i+1][j]) / 2.0` (east)
            *   `eps_w = (eps_r[i][j] + eps_r[i-1][j]) / 2.0` (west)
            *   `eps_n = (eps_r[i][j] + eps_r[i][j+1]) / 2.0` (north)
            *   `eps_s = (eps_r[i][j] + eps_r[i][j-1]) / 2.0` (south)
        *   `sum_eps = eps_e + eps_w + eps_n + eps_s`. If `sum_eps` is zero (can happen in unusual configurations, though unlikely here), skip update.
        *   **Gauss-Seidel Value (`val_GS`)**: Calculates the updated potential based on neighbors using the FDM stencil for Poisson's equation with variable permittivity:
            `val_GS = (eps_e * V[i+1][j] + eps_w * V[i-1][j] + eps_n * V[i][j+1] + eps_s * V[i][j-1]) / sum_eps`
        *   `v_old_ij = V[i][j]`.
        *   **SOR Update**: `V[i][j] = (1.0 - omega) * v_old_ij + omega * val_GS`.
        *   `current_point_diff = std::abs(V[i][j] - v_old_ij)`.
        *   Updates `max_diff_iter = std::max(max_diff_iter, current_point_diff)`.
    *   **Phase 2: Update "Black" Points (`(i+j) % 2 != 0`)**:
        *   Identical logic to Phase 1, but for "black" points.
*   **Neumann Boundary Conditions (Outer Simulation Box)**:
    *   Applied to the four outer boundaries of the simulation domain (i=0, i=Nx-1, j=0, j=Ny-1).
    *   `dV/dn = 0` implies the potential at the boundary is equal to the potential of its immediate interior neighbor.
    *   Example: `V[0][j_idx] = V[1][j_idx]` (left boundary), if not a fixed potential point.
    *   These loops are also parallelized with `#pragma omp parallel for`.
*   **Re-apply Fixed Potentials**:
    *   A loop ensures that points marked by `fixed_potential_mask` retain their assigned fixed values. This is crucial if any fixed potential points coincide with the Neumann boundaries (though not the case in this specific geometry's Dirichlet setup).
*   **Convergence Check & Output**:
    *   Every 100 iterations, prints the current iteration number and `max_diff_iter`.
    *   If `max_diff_iter < tolerance`, convergence is achieved, a message is printed, and the loop breaks.
    *   If `iteration == max_iterations - 1`, the maximum number of iterations is reached, a message is printed.

### 3.8. Calculate Electric Field (`Ex`, `Ey`)

*   `Ex`, `Ey`: 2D vectors to store x and y components of the electric field.
*   **Interior Points**: Calculated using the central difference formula:
    *   `Ex[i][j] = -(V[i+1][j] - V[i-1][j]) / (2 * h)`
    *   `Ey[i][j] = -(V[i][j+1] - V[i][j-1]) / (2 * h)`
*   **Boundary Points**: Calculated using forward or backward difference formulas:
    *   `Ex[0][j] = -(V[1][j] - V[0][j]) / h` (forward for left boundary)
    *   `Ex[Nx-1][j] = -(V[Nx-1][j] - V[Nx-2][j]) / h` (backward for right boundary)
    *   Similar logic for `Ey` at top/bottom boundaries.

### 3.9. Output Results to CSV

*   Calls `saveToCSV` to save:
    *   `V` (potential) to "potential.csv".
    *   `Ex` (x-component of electric field) to "electric_field_x.csv".
    *   `Ey` (y-component of electric field) to "electric_field_y.csv".
    *   `eps_r` (permittivity map) to "permittivity.csv".
*   Calls `saveCoordinatesToCSV` to save:
    *   `x_coords` to "x_coordinates.csv".
    *   `y_coords` to "y_coordinates.csv".

### 3.10. Summary Message

*   Prints a summary message indicating that results have been saved and can be plotted using external tools.

### 3.11. Return 0

*   Indicates successful program execution.

## Compilation and Execution Notes

*   To compile with OpenMP support, you typically need to add a compiler flag (e.g., `-fopenmp` for GCC/Clang, `/openmp` for MSVC).
*   The program expects to be able to create a subdirectory (e.g., "geometria_Denti_uguali") in the directory where it is run.

This detailed explanation should provide a good understanding of the code's functionality and structure.
