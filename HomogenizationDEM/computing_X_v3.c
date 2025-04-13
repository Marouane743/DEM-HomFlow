/****************************************************************************
 * INTEGRAL COMPUTATION FOR DEM HOMOGENIZATION
 *
 * Author: Marouane ELBISSSOURI
 * Date: 2024-10-01
 *
 * Description:
 *   This program pre-computes integral values needed for the Goldhirsch
 *   homogenization method. It generates lookup tables for different angles
 *   (theta) that will be used by the homogenization code to efficiently
 *   calculate stress tensors from particle interactions.
 *
 *   The program uses OpenMP for parallel computation to significantly
 *   speed up the generation of these integral tables.
 ****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>  // For parallel processing

#include "shared_config.h"

/**
 * Weight function phi(q) for the coarse-graining approach
 *
 * This function implements Lucy's weight function which is used in the
 * coarse-graining approach to distribute particle properties to the grid.
 *
 * @param q     Normalized distance (q = 2*r/RC where r is distance and RC is cutoff radius)
 * @param s_3   Normalization factor (s_3 = 1/(π*(RC/2)³))
 * @return      Weight value at distance q
 */
double phi(double q, double s_3) {
    if(q >= 0.0 && q <= 1.0){
        // Polynomial form for 0 ≤ q ≤ 1
        return s_3 * (1.0 - 1.5*q*q*(1.0 - q/2.0));
    } else if(q > 1.0 && q <= 2.0){
        // Polynomial form for 1 < q ≤ 2
        double term = 2.0 - q;
        return (s_3 / 4.0) * term*term*term;
    }
    // Zero outside the support (q > 2)
    return 0.0;
}

/**
 * Calculates the normalized distance for the integration
 *
 * This function computes the normalized distance between two points
 * based on the integration parameter s, distances x and dj, and angle theta.
 *
 * @param s      Integration parameter (0 ≤ s ≤ 1)
 * @param x      Distance from grid point to first particle
 * @param dj     Distance between particles
 * @param theta  Angle between vectors
 * @param RC     Cutoff radius
 * @return       Normalized distance value
 */
double norm_func(double s, double x, double dj, double theta, double RC) {
    // Calculate the three terms of the distance formula
    double term1 = (1.0 - s)*(1.0 - s)*x*x;                // Distance contribution from first particle
    double term2 = 2.0 * (1.0 - s)*s * dj*cos(theta)*x;    // Cross-term with angle dependency
    double term3 = s*s * dj*dj;                            // Distance contribution from second particle

    // Sum the terms and check for numerical errors
    double inside_sqrt = term1 + term2 + term3;
    if(inside_sqrt < 0.0){
        // If negative due to numerical error, return 0 instead of NaN
        return 0.0;
    }

    // Return normalized distance (divided by cutoff radius)
    return sqrt(inside_sqrt)/RC;
}

/**
 * Calculates the integrand value at a specific point
 *
 * This function computes the value of the integrand at a specific point
 * defined by the parameters s, x, dj, and theta.
 *
 * @param s      Integration parameter (0 ≤ s ≤ 1)
 * @param x      Distance from grid point to first particle
 * @param dj     Distance between particles
 * @param theta  Angle between vectors
 * @param s_3    Normalization factor for phi function
 * @param RC     Cutoff radius
 * @return       Value of the integrand at the specified point
 */
double integrand(double s, double x, double dj, double theta, double s_3, double RC){
    // Calculate normalized distance q = 2*r/RC
    double q = 2.0 * norm_func(s, x, dj, theta, RC);

    // Return the weight function value at this distance
    return phi(q, s_3);
}

/**
 * Performs numerical integration using Simpson's rule
 *
 * This function computes the integral of phi over the parameter s
 * using Simpson's 1/3 rule for numerical integration.
 *
 * @param x      Distance from grid point to first particle
 * @param dj     Distance between particles
 * @param theta  Angle between vectors
 * @param s_3    Normalization factor for phi function
 * @param RC     Cutoff radius
 * @return       Result of the numerical integration
 */
double integral_phi(double x, double dj, double theta, double s_3, double RC){
    // Number of intervals for Simpson's rule (must be even)
    int iter = 100;
    if(iter % 2 !=0){
        iter +=1;  // Ensure even number of intervals
    }

    // Step size for integration
    double h = (S_MAX - S_MIN)/(double)iter;
    double sum = 0.0;

    // Apply Simpson's 1/3 rule: f(a) + 4*f(x₁) + 2*f(x₂) + 4*f(x₃) + ... + f(b)
    for(int i=0; i<=iter; i++){
        // Calculate current s value
        double s = S_MIN + i*h;

        // Determine weight based on Simpson's rule pattern
        double weight = 1.0;
        if(i==0 || i==iter) weight = 1.0;      // First and last points
        else if(i % 2 == 0) weight=2.0;        // Even indices
        else weight=4.0;                       // Odd indices

        // Calculate integrand value at this point
        double val = integrand(s, x, dj, theta, s_3, RC);
        sum += weight * val;
    }

    // Apply Simpson's rule formula: (h/3) * [f(a) + 4*f(x₁) + 2*f(x₂) + ... + f(b)]
    return (h/3.0)*sum;
}

/**
 * Main function for integral computation
 *
 * This function handles command-line arguments, initializes parameters,
 * and coordinates the parallel computation of integral values for different
 * angles (theta). The results are saved as CSV files in the INTEGRALS_DIR directory.
 *
 * @param argc  Number of command-line arguments
 * @param argv  Array of command-line arguments
 * @return      0 on success, non-zero on error
 */
int main(int argc, char *argv[]){

    // Start timing for performance measurement
    double start_time = omp_get_wtime();
    double RC;

    // Validate command-line arguments
    if(argc < 2){
        fprintf(stderr, "Usage: %s <RC_value>\n", argv[0]);
        return 1;
    }

    // Parse and validate the cutoff radius (RC)
    RC = atof(argv[1]);
    if(RC <= 0.0){
        fprintf(stderr, "Error: RC must be a positive number.\n");
        return 1;
    }

    // Calculate normalization factor and discretization steps
    double s_3 = 1.0/(PI*pow(RC/2.0,3.0));     // Normalization factor for phi function
    double theta_step = PI/(double)NUM_THETA;   // Angular step size
    double x_step = (2.0*RC)/(double)NUM_X;     // Distance step for x
    double dj_step= (2.0*RC)/(double)NUM_DJ;    // Distance step for dj

    // Clean up previous integral files
    char command[512];
    snprintf(command, sizeof(command), "rm -f %s/*", INTEGRALS_DIR);
    if(system(command) != 0){
        fprintf(stderr, "Warning: Failed to remove some files in %s\n", INTEGRALS_DIR);
    }

    // Create the integrals directory if it doesn't exist
    system("mkdir -p " INTEGRALS_DIR);

    // Process each theta value in parallel
    #pragma omp parallel for schedule(dynamic)
    for(int theta_idx=0; theta_idx<NUM_THETA; theta_idx++){
        // Calculate current theta value
        double theta = theta_idx * theta_step;

        // Create output file for this theta value
        char filename[256];
        sprintf(filename, INTEGRALS_DIR "/int_theta_%03d.csv", theta_idx);
        FILE* fp = fopen(filename,"w");
        if(!fp){
            printf("Error opening %s\n", filename);
            continue;
        }

        // Write CSV header row with dj values
        fprintf(fp,"x");
        for(int dj_idx=0; dj_idx<NUM_DJ; dj_idx++){
            double dj = dj_idx * dj_step;
            fprintf(fp,",%.6f", dj);
        }
        fprintf(fp,"\n");

        // Allocate memory for integral values
        double* integrals = (double*)malloc(NUM_X*NUM_DJ*sizeof(double));
        if(!integrals){
            fclose(fp);
            continue;
        }

        // Calculate integral values in parallel for all x and dj combinations
        #pragma omp parallel for collapse(2)
        for(int x_idx=0; x_idx<NUM_X; x_idx++){
            for(int dj_idx=0; dj_idx<NUM_DJ; dj_idx++){
                double x = x_idx * x_step;
                double dj = dj_idx* dj_step;

                // Only calculate when dj >= x (physical constraint)
                if(dj >= x){
                    double val = integral_phi(x, dj, theta, s_3, RC);
                    integrals[x_idx*NUM_DJ + dj_idx] = val;
                } else {
                    // Mark as not applicable
                    integrals[x_idx*NUM_DJ + dj_idx] = NAN;
                }
            }
        }

        // Write calculated integral values to CSV file
        for(int x_idx=0; x_idx<NUM_X; x_idx++){
            double x = x_idx*x_step;
            fprintf(fp,"%.6f", x);

            for(int dj_idx=0; dj_idx<NUM_DJ; dj_idx++){
                double val = integrals[x_idx*NUM_DJ + dj_idx];
                double dj  = dj_idx*dj_step;

                if(dj >= x){
                    fprintf(fp,",%.6f", val);  // Write calculated value
                } else {
                    fprintf(fp,",");           // Empty cell for invalid combinations
                }
            }
            fprintf(fp,"\n");
        }

        // Clean up resources
        free(integrals);
        fclose(fp);
        printf("File %s generated.\n", filename);
    }

    // End timing and calculate elapsed time
    double end_time = omp_get_wtime();
    double elapsed_time = end_time - start_time;
    printf("Total computation time: %.2f seconds\n", elapsed_time);

    return 0;
}
