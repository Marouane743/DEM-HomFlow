# HomogenizationDEM

## Project Overview

HomogenizationDEM is a computational framework for applying the Goldhirsch homogenization method to Discrete Element Method (DEM) simulations. This project enables the calculation of continuum fields (density, velocity, stress tensor, and viscosity) from discrete particle data, providing a bridge between microscopic particle interactions and macroscopic material behavior.

## Author

**Marouane ELBISSSOURI**  
*Date: 2024-10-01*

## Project Architecture

The project consists of three main components:

1. **Integral Computation (`computing_X_v3.c`)**: Pre-computes integral lookup tables for efficient homogenization
2. **Homogenization Process (`homogenization_nv1.c`)**: Processes DEM simulation data to calculate continuum fields
3. **Shared Configuration (`shared_config.h`)**: Contains common parameters and constants used by both components

### Directory Structure

```
HomogenizationDEM/
├── computing_X_v3.c         # Integral computation program
├── homogenization_nv1.c     # Main homogenization program
├── shared_config.h          # Shared configuration header
├── integrals_csv/           # Directory for storing integral lookup tables
├── DEM_data/                # Default input directory for DEM simulation data
└── output_data/             # Default output directory for homogenized fields
```

## Theoretical Background

The Goldhirsch homogenization method is a coarse-graining approach that transforms discrete particle data into continuous field variables. It uses a weight function to distribute particle properties to a computational grid, allowing for the calculation of:

- **Density field**: Distribution of mass throughout the domain
- **Velocity field**: Continuous representation of particle velocities
- **Stress tensor field**: Distribution of internal forces within the material
- **Viscosity field**: Measure of the material's resistance to deformation

## Workflow

### Step 1: Integral Computation

The `computing_X_v3.c` program pre-computes integral values needed for the homogenization process:

1. Creates a discretization of angles (theta) from 0 to π
2. For each angle, computes integral values for different distances
3. Stores these values in CSV files in the `integrals_csv` directory
4. Uses OpenMP for parallel computation to speed up the process

### Step 2: Homogenization

The `homogenization_nv1.c` program processes DEM simulation data:

1. Reads particle data from CSV files in the `DEM_data` directory
2. Loads the pre-computed integral tables
3. Creates a computational grid covering the simulation domain
4. Applies the coarse-graining approach to calculate continuum fields
5. Saves the results to CSV files in the `output_data` directory

## Input Data Format

The homogenization program expects the following input files in the `DEM_data` directory:

1. **Contact_pairs_xxxx.csv**: Contains particle contact information
   - Columns: Particle A ID, Particle B ID, Force components (x,y,z)

2. **DEMdemo_output_xxxx.csv**: Contains particle positions and velocities
   - Columns: Particle ID, Position (x,y,z), Radius, Velocity (vx,vy,vz)

3. **header0000.txt**: Contains simulation parameters
   - Domain dimensions
   - Particle density
   - Maximum particle radius
   - Simulation dimension (2D or 3D)

## Output Data Format

The homogenization process generates `output_data_xxxx.csv` files in the `output_data` directory with the following columns:

1. Grid point coordinates (x,y,z)
2. Density at the grid point
3. Velocity components (vx,vy,vz)
4. Stress tensor components (σxx, σxy, σxz, σyy, σyz, σzz)
5. Viscosity at the grid point

## Configuration Parameters

Key parameters in `shared_config.h`:

- **NUM_THETA**: Number of angular discretization points (default: 100)
- **NUM_X**: Number of distance discretization points (default: 2000)
- **NUM_DJ**: Number of particle-particle interaction distance points (default: 2000)
- **INTEGRALS_DIR**: Directory for integral CSV files (default: "integrals_csv")
- **MAX_LINE_LENGTH**: Maximum length for reading lines from input files (default: 100000)

## Compilation and Execution

### Compiling the Integral Computation Program

```bash
gcc -o compute_integrals computing_X_v3.c -lm -fopenmp
```

### Running the Integral Computation

```bash
./compute_integrals
```

### Compiling the Homogenization Program

```bash
gcc -o homogenize homogenization_nv1.c -lm
```

### Running the Homogenization Process

```bash
./homogenize [start_index] [end_index]
```

Where:
- `start_index`: First DEM data file index to process (optional, default: 0)
- `end_index`: Last DEM data file index to process (optional, default: 9999)

## Performance Considerations

- The integral computation is computationally intensive but only needs to be performed once
- OpenMP parallelization significantly speeds up the integral computation
- The homogenization process scales with the number of particles and grid points
- Memory usage depends on the grid resolution and number of particles

## References

1. Goldhirsch, I. (2010). Stress, stress asymmetry and couple stress: from discrete particles to continuous fields. *Granular Matter*, 12(3), 239-252.
2. Lucy, L. B. (1977). A numerical approach to the testing of the fission hypothesis. *The Astronomical Journal*, 82, 1013-1024.
