/****************************************************************************
 * SHARED CONFIGURATION HEADER
 *
 * Author: Marouane ELBISSSOURI
 * Date: 2024-10-01
 *
 * Description:
 *   This header file contains shared constants, macros, and configuration
 *   parameters used across the DEM homogenization project. It ensures
 *   consistency between the integral computation (computing_X_v3.c) and
 *   the homogenization process (homogenization_nv1.c).
 ****************************************************************************/

#ifndef SHARED_CONFIG_H
#define SHARED_CONFIG_H

/*
 * DISCRETIZATION PARAMETERS
 * These parameters control the resolution of our numerical integration
 * and the accuracy of the homogenization process.
 */

// Number of angular discretization points (theta values from 0 to π)
// Higher values provide more accurate angular resolution for integral calculations
#ifndef NUM_THETA
  #define NUM_THETA 100
#endif

// Number of distance discretization points for x values (0 to 2*RC)
// Controls the spatial resolution of the integral lookup tables
#ifndef NUM_X
  #define NUM_X 2000
#endif

// Number of distance discretization points for dj values (0 to 2*RC)
// Controls the resolution for particle-particle interaction distances
#ifndef NUM_DJ
  #define NUM_DJ 2000
#endif

// Number of integral files to generate (one per theta value)
// Must match NUM_THETA for proper angular coverage
#define NUM_FILES NUM_THETA

// Matrix size for integral lookup tables
// Using NUM_X ensures consistent discretization between computing and lookup
#define MATRIX_SIZE NUM_X

/*
 * MATHEMATICAL CONSTANTS
 */
// Mathematical constant π (pi) for angular calculations
#ifndef M_PI
  #define M_PI 3.14159265358979323846
#endif

// Alternative definition of π for backward compatibility
#ifndef PI
  #define PI 3.14159265358979323846
#endif

/*
 * FILE AND DIRECTORY SETTINGS
 */
// Directory where integral CSV files will be stored
#define INTEGRALS_DIR "integrals_csv"

// Maximum length for reading lines from input files
// Increased to handle large CSV files with many columns
#define MAX_LINE_LENGTH 100000

/*
 * INTEGRATION PARAMETERS
 */
// Upper bound for integration parameter s
#define S_MAX 1.00

// Lower bound for integration parameter s
#define S_MIN 0.00

#endif /* SHARED_CONFIG_H */
