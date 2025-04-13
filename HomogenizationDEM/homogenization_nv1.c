/****************************************************************************
 * DISCRETE ELEMENT METHOD (DEM) HOMOGENIZATION
 *
 * Author: Marouane ELBISSSOURI
 * Date: 2024-10-01
 *
 * Description:
 *   This program implements the Goldhirsch homogenization method for DEM simulations.
 *   It processes multiple DEM file sets in the default directory "DEM_data".
 *   For each index, it reads:
 *     1) Contact_pairs_xxxx.csv - Contains particle contact information
 *     2) DEMdemo_output_xxxx.csv - Contains particle positions and velocities
 *     3) header0000.txt - Contains simulation parameters
 *
 *   The program then applies the Goldhirsch homogenization algorithm to compute:
 *     - Density field
 *     - Velocity field
 *     - Stress tensor field
 *     - Viscosity field
 *
 *   Results are saved as output_data_xxxx.csv files in the "output_data" directory.
 *   The program measures elapsed time for each data set and the total processing time.
 ****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>   // For access() on POSIX
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>

// -------------------- Constants or Macros --------------------
#include "shared_config.h"
#define FILENAME_FORMAT INTEGRALS_DIR "/int_theta_%03d.csv"
#define EPSILON 1e-8


// Maximum index based on your dataset
#define MAX_INDEX 9999

// Default directories
#define DEFAULT_INPUT_DIR "DEM_data"
#define DEFAULT_OUTPUT_DIR "output_data"

// -------------------- Data Structures --------------------
typedef struct {
    double x, y, z;
} Point3d;

// Symmetric 3D stress tensor
typedef struct {
    double xx, xy, xz;
    double yx, yy, yz;
    double zx, zy, zz;
} StressTensor;

// Particle structure
typedef struct {
    int id;
    double r;        // Radius
    double mass;     // Mass
    Point3d p;       // Position
    Point3d v;       // Velocity
} Particle3d;

// -------------------- Global for integral matrices --------------------
static double ***integral_matrices = NULL;

// -------------------- Function Declarations --------------------
double distance3d(Point3d p1, Point3d p2);
Point3d vectorSubtract(Point3d pA, Point3d pB);
double dotProduct(Point3d pA, Point3d pB);
int isNull(Point3d p);
double computeMass(double radius, double density);
double phi(double q, double s_3);

// Matrix loading routines
double** allocateMatrix(int size);
void freeMatrix(double **matrix, int size);
int splitCSV(const char *line, char **tokens, int max_tokens);
double** loadMatrix(const char *filename, int size);
int loadAllMatrices(int num_files, int size);
void freeAllMatrices(int num_files, int size);

// Progress bar
void display_progress(double progress);

/*****************************************************************************
 * processSingleSet()
 *   Core homogenization logic that was previously in your main().
 *   Takes filenames as arguments to process multiple "xxxx" indices.
 *****************************************************************************/
int processSingleSet(const char *headerFile,
                     const char *demFile,
                     const char *contactFile,
                     const char *outFile,
                     double Rc)
{
    // ------------------ Start timing for this set ------------------
    clock_t local_start = clock();

    // =========== SECTION 1: Validate Rc and s_3 ===========
    if (Rc <= 0.0) {
        fprintf(stderr, "Error: Rc must be a positive number.\n");
        return 1;
    }
    double s_3 = 1 / (M_PI * pow(Rc / 2, 3));
    printf("s_3 = %f\n", s_3);

    // =========== SECTION 2: Read Simulation Parameters ===========
    int Np = 0;               // Number of particles
    double Density = 0.0;     // Particle density
    double RadiusMax = 0.0;   // Maximum particle radius
    int D = 0;                // Dimension
    double xsmin = 0.0, xsmax = 0.0, ysmin = 0.0, ysmax = 0.0, zsmin = 0.0, zsmax = 0.0;

    FILE *file = fopen(headerFile, "r");
    if (!file) {
        perror("Error opening header file");
        return 1;
    }

    char line[MAX_LINE_LENGTH];
    while (fgets(line, sizeof(line), file)) {
        if (strncmp(line, "@ Np", 4) == 0) {
            if (fgets(line, sizeof(line), file))
                Np = atoi(line);
        } else if (strncmp(line, "@ RadiusMax", 10) == 0) {
            if (fgets(line, sizeof(line), file))
                RadiusMax = atof(line);
        } else if (strncmp(line, "@ Density", 9) == 0) {
            if (fgets(line, sizeof(line), file))
                Density = atof(line);
        } else if (strncmp(line, "@ Dimension", 11) == 0) {
            if (fgets(line, sizeof(line), file))
                D = atoi(line);
        } else if (D >= 1 && strncmp(line, "@ xsmin", 7) == 0) {
            if (fgets(line, sizeof(line), file))
                xsmin = atof(line);
        } else if (D >= 1 && strncmp(line, "@ xsmax", 7) == 0) {
            if (fgets(line, sizeof(line), file))
                xsmax = atof(line);
        } else if (D >= 2 && strncmp(line, "@ ysmin", 7) == 0) {
            if (fgets(line, sizeof(line), file))
                ysmin = atof(line);
        } else if (D >= 2 && strncmp(line, "@ ysmax", 7) == 0) {
            if (fgets(line, sizeof(line), file))
                ysmax = atof(line);
        } else if (D == 3 && strncmp(line, "@ zsmin", 7) == 0) {
            if (fgets(line, sizeof(line), file))
                zsmin = atof(line);
        } else if (D == 3 && strncmp(line, "@ zsmax", 7) == 0) {
            if (fgets(line, sizeof(line), file))
                zsmax = atof(line);
        } else if (D == 3 && strncmp(line, "@ zpmax", 7) == 0) {
            if (fgets(line, sizeof(line), file))
                zsmax = atof(line);
        }
    }
    fclose(file);
    printf("Density = %f\n", Density);

    if (D != 3 && D != 2) {
        printf("Warning: This code is set for 3D and 2D but D = %d.\n", D);
    }

    if (D == 3){
        if (zsmax <= 0) {
            printf("No particles left in the hopper! zsmax = %.5f\n", zsmax);
            return 1;
        }
    }

    // =========== SECTION 3: Initialize Grid ===========
    double step = 0.02;
    double xgmin = xsmin + 2 * Rc;
    double xgmax = xsmax - 2 * Rc;
    double ygmin = ysmin + 2 * Rc;
    double ygmax = ysmax - 2 * Rc;
    double zgmin = 0.0;
    double zgmax = 0.0;
    if (D == 2){
        zgmin = 0;
        zgmax = 0;
    }
    else {
        zgmin = 2 * Rc;
        zgmax = zsmax - 2 * Rc;
    }

    printf("zgmax = %f\n", zgmax);

    int Cx = (int)((xgmax - xgmin) / step);
    int Cy = (int)((ygmax - ygmin) / step);
    int Cz = (int)((zgmax - zgmin) / step);
    int N = (Cx + 1) * (Cy + 1) * (Cz + 1);
    printf("N = %d\n", N);

    Point3d *Grid = (Point3d *)malloc(N * sizeof(Point3d));
    if (!Grid) {
        perror("Error allocating Grid");
        return 1;
    }

    for (int k = 0; k <= Cz; k++) {
        for (int j = 0; j <= Cy; j++) {
            for (int i = 0; i <= Cx; i++) {
                int idx = i + j * (Cx + 1) + k * (Cx + 1) * (Cy + 1);
                Grid[idx].x = xgmin + step * i;
                Grid[idx].y = ygmin + step * j;
                Grid[idx].z = zgmin + step * k;
            }
        }
    }

    // =========== SECTION 4: Initialize Fields ===========
    int w = (int)(Rc / step);
    Point3d *Momentum       = (Point3d *)calloc(N, sizeof(Point3d));
    double  *DensityField   = (double  *)calloc(N, sizeof(double));
    Point3d *VelocityField  = (Point3d *)calloc(N, sizeof(Point3d));
    StressTensor *STF       = (StressTensor *)calloc(N, sizeof(StressTensor));
    Point3d *ParticlePos    = (Point3d *)malloc(Np * sizeof(Point3d));
    Point3d *Particleout_contact = (Point3d *)malloc(Np * sizeof(Point3d));
    /* Allocate array for viscosity (to be computed later) */
    double *ViscosityField = (double *)calloc(N, sizeof(double));

    if (!Momentum || !DensityField || !VelocityField || !STF || !ParticlePos || !ViscosityField) {
        fprintf(stderr, "Error allocating field arrays.\n");
        free(Grid);
        free(Momentum); free(DensityField);
        free(VelocityField); free(STF); free(ParticlePos); free(ViscosityField);
        return 1;
    }

    // =========== SECTION 5: Compute Density and Momentum Fields ===========
    FILE *particleFile = fopen(demFile, "r");
    if (!particleFile) {
        perror("Error opening DEM file");
        freeAllMatrices(NUM_FILES, MATRIX_SIZE);
        free(Grid);
        free(Momentum); free(DensityField);
        free(VelocityField); free(STF); free(ParticlePos); free(ViscosityField);
        return 1;
    }

    // Skip header line
    if (!fgets(line, sizeof(line), particleFile)) {
        fprintf(stderr, "Error reading header from %s\n", demFile);
        fclose(particleFile);
        freeAllMatrices(NUM_FILES, MATRIX_SIZE);
        free(Grid);
        free(Momentum); free(DensityField);
        free(VelocityField); free(STF); free(ParticlePos); free(ViscosityField);
        return 1;
    }

    int particleCount = 0;
    int particleoutgrid = 0;
    while (fgets(line, sizeof(line), particleFile)) {
        Particle3d Particle;
        double x, y, z, r, vx, vy, vz;
        if (sscanf(line, "%lf,%lf,%lf,%lf,%lf,%lf,%lf",
                   &x, &y, &z, &r, &vx, &vy, &vz) == 7)
        {
            Particle.p.x = x;
            Particle.p.y = y;
            Particle.p.z = z;
            Particle.r   = r;
            Particle.mass= computeMass(r, Density);
            Particle.v.x = vx;
            Particle.v.y = vy;
            Particle.v.z = vz;

            Particleout_contact[particleCount] = (Point3d){0,0,0};

            // Check if particle is inside the bounding box of the grid
            if (x < xgmin || x > xgmax || y < ygmin || y > ygmax || z < zgmin || z > zgmax) {
                ParticlePos[particleCount] = (Point3d){0, 0, 0};
                Particleout_contact[particleCount] = Particle.p;
                particleCount++;
                particleoutgrid++;
                continue;
            }

            ParticlePos[particleCount] = Particle.p;

            // Indices for local region
            int nxmin = (int)((x - xgmin) / step) - w;
            int nymin = (int)((y - ygmin) / step) - w;
            int nzmin = (int)((z - zgmin) / step) - w;
            if (nxmin < 0) nxmin = 0;
            if (nymin < 0) nymin = 0;
            if (nzmin < 0) nzmin = 0;

            int nxmax = (int)((x - xgmin) / step) + w + 1;
            int nymax = (int)((y - ygmin) / step) + w + 1;
            int nzmax = (int)((z - zgmin) / step) + w + 1;
            if (nxmax > Cx) nxmax = Cx;
            if (nymax > Cy) nymax = Cy;
            if (nzmax > Cz) nzmax = Cz;

            for (int c = nzmin; c <= nzmax; c++) {
                for (int b = nymin; b <= nymax; b++) {
                    for (int a = nxmin; a <= nxmax; a++) {
                        int idx = a + b * (Cx + 1) + c * (Cx + 1) * (Cy + 1);
                        Point3d rg = Grid[idx];
                        double q   = 2.0 * distance3d(rg, Particle.p) / Rc;

                        if (q >= 2.0) continue; // Outside kernel
                        double ph = phi(q, s_3);

                        DensityField[idx]   += Particle.mass * ph;
                        Momentum[idx].x     += Particle.mass * Particle.v.x * ph;
                        Momentum[idx].y     += Particle.mass * Particle.v.y * ph;
                        Momentum[idx].z     += Particle.mass * Particle.v.z * ph;

                        // Kinetic contribution
                        STF[idx].xx -= Particle.mass * Particle.v.x * Particle.v.x * ph;
                        STF[idx].xy -= Particle.mass * Particle.v.x * Particle.v.y * ph;
                        STF[idx].xz -= Particle.mass * Particle.v.x * Particle.v.z * ph;
                        STF[idx].yy -= Particle.mass * Particle.v.y * Particle.v.y * ph;
                        STF[idx].yz -= Particle.mass * Particle.v.y * Particle.v.z * ph;
                        STF[idx].zz -= Particle.mass * Particle.v.z * Particle.v.z * ph;
                    }
                }
            }

            particleCount++;
            if (particleCount >= Np) {   // Safety
               printf("Warning : particleCount is lager than the number of particles ! \n");
                break;
            }
        }
    }
    fclose(particleFile);

    int particleingrid = particleCount - particleoutgrid;
    printf("Number of particles outside the grid : %d \n", particleoutgrid);
    printf("Number of particles inside the grid : %d \n", particleingrid);

    // =========== SECTION 6: Compute Velocity Field & Kinetic Stress Tensor ===========
    for (int i = 0; i < N; i++) {
        if (DensityField[i] > EPSILON) {
            double invRho = 1.0 / DensityField[i];
            // Mean velocity
            VelocityField[i].x = Momentum[i].x * invRho;
            VelocityField[i].y = Momentum[i].y * invRho;
            VelocityField[i].z = Momentum[i].z * invRho;

            // Stress tensor correction:
            // Term (1): Contribution from the mean velocity field
            STF[i].xx += DensityField[i] * VelocityField[i].x * VelocityField[i].x;
            STF[i].xy += DensityField[i] * VelocityField[i].x * VelocityField[i].y;
            STF[i].xz += DensityField[i] * VelocityField[i].x * VelocityField[i].z;
            STF[i].yy += DensityField[i] * VelocityField[i].y * VelocityField[i].y;
            STF[i].yz += DensityField[i] * VelocityField[i].y * VelocityField[i].z;
            STF[i].zz += DensityField[i] * VelocityField[i].z * VelocityField[i].z;

            // Symmetry
            STF[i].yx = STF[i].xy;
            STF[i].zx = STF[i].xz;
            STF[i].zy = STF[i].yz;
        }
        else {
            VelocityField[i].x = 0.0;
            VelocityField[i].y = 0.0;
            VelocityField[i].z = 0.0;
        }
    }

    // =========== SECTION 7: Potential Contribution to Stress Tensor (Contacts) ===========
    FILE *cFile = fopen(contactFile, "r");
    if (!cFile) {
        perror("Error opening contact file");
        freeAllMatrices(NUM_FILES, MATRIX_SIZE);
        free(Grid);
        free(Momentum); free(DensityField);
        free(VelocityField); free(STF); free(ParticlePos); free(ViscosityField);
        return 1;
    }
    // Skip header
    if (!fgets(line, sizeof(line), cFile)) {
        fprintf(stderr, "Error reading header from %s\n", contactFile);
        fclose(cFile);
        freeAllMatrices(NUM_FILES, MATRIX_SIZE);
        free(Grid);
        free(Momentum); free(DensityField);
        free(VelocityField); free(STF); free(ParticlePos); free(ViscosityField);
        return 1;
    }

    int T = 0;
    while (fgets(line, sizeof(line), cFile)) {
        T = 1;
        char contact_type[10];
        int a_idx, b_idx;
        double f_x, f_y, f_z;
        double x_c, y_c, z_c;

        if (sscanf(line, "%9[^,],%d,%d,%lf,%lf,%lf,%lf,%lf,%lf",
                   contact_type, &a_idx, &b_idx,
                   &f_x, &f_y, &f_z,
                   &x_c, &y_c, &z_c) == 9)
        {
            if (strcmp(contact_type, "SS") != 0) {
                break; // skip non-"SS" lines
            }

            // Adjust indices if needed:
            int A = a_idx - 4;
            int B = b_idx - 4;
            if (A < 0 || A >= Np || B < 0 || B >= Np) {
                printf("Warning with Np! \n");
                continue;
            }

            // If stored position is zero, skip
            if ((isNull(ParticlePos[A]) == 1 && isNull(ParticlePos[B]) == 1)) {
                continue; // skipping the particles that are not in the grid
            }

            Point3d rpA = ParticlePos[A];
            Point3d rpB = ParticlePos[B];
            Point3d F   = { f_x, f_y, f_z };

            if (isNull(rpA) == 1) {
                rpA = Particleout_contact[A];
            }
            if (isNull(rpB) == 1) {
                rpB = Particleout_contact[B];
            }

            if (isNull(rpA) == 1 || isNull(rpB) == 1) {
                printf("Big problem !\n");
                return 2;
            }

            Point3d rpAB = vectorSubtract(rpA, rpB);
            // Rough bounding box for this contact
            double x_min = fmin(rpA.x, rpB.x) - Rc;
            double x_max = fmax(rpA.x, rpB.x) + Rc;
            double y_min = fmin(rpA.y, rpB.y) - Rc;
            double y_max = fmax(rpA.y, rpB.y) + Rc;
            double z_min = fmin(rpA.z, rpB.z) - Rc;
            double z_max = fmax(rpA.z, rpB.z) + Rc;

            int nxmin = (int)((x_min - xgmin) / step); if (nxmin < 0) nxmin = 0;
            int nymin = (int)((y_min - ygmin) / step); if (nymin < 0) nymin = 0;
            int nzmin = (int)((z_min - zgmin) / step); if (nzmin < 0) nzmin = 0;

            int nxmax = (int)((x_max - xgmin) / step); if (nxmax > Cx) nxmax = Cx;
            int nymax = (int)((y_max - ygmin) / step); if (nymax > Cy) nymax = Cy;
            int nzmax = (int)((z_max - zgmin) / step); if (nzmax > Cz) nzmax = Cz;

            // For each grid node in bounding region, do integral lookups
            for (int cz = nzmin; cz <= nzmax; cz++) {
                for (int by = nymin; by <= nymax; by++) {
                    for (int ax = nxmin; ax <= nxmax; ax++) {
                        int idx = ax + by * (Cx + 1) + cz * (Cx + 1) * (Cy + 1);
                        Point3d rg = Grid[idx];

                        // Example geometry for lookup
                        double xA  = distance3d(rg, rpA);
                        double xB  = distance3d(rg, rpB);

                        // Dot-based angle
                        double dotAB = dotProduct(vectorSubtract(rpA, rg),
                                                  vectorSubtract(rpB, rg));
                        double lenA  = xA;
                        double lenB  = xB;
                        double cosT  = 0.0;
                        if (lenA > EPSILON && lenB > EPSILON) {
                            cosT = dotAB / (lenA * lenB);
                        }
                        if (cosT > 1.0)  cosT = 1.0;
                        if (cosT < -1.0) cosT = -1.0;

                        double theta = acos(cosT);

                        // Indices for integral_matrices
                        int idxX = (int)(xA * MATRIX_SIZE / (2.0 * Rc));
                        int idxD = (int)(xB * MATRIX_SIZE / (2.0 * Rc));
                        int idxT = (int)(theta * NUM_FILES / M_PI);
                        if (idxX < 0) idxX = 0; if (idxX >= MATRIX_SIZE) idxX = MATRIX_SIZE - 1;
                        if (idxD < 0) idxD = 0; if (idxD >= MATRIX_SIZE) idxD = MATRIX_SIZE - 1;
                        if (idxT < 0) idxT = 0; if (idxT >= NUM_FILES ) idxT = NUM_FILES - 1;

                        double integral_value;
                        if (idxD > idxX) {
                            integral_value = integral_matrices[idxT][idxX][idxD];
                        } else {
                            integral_value = integral_matrices[idxT][idxD][idxX];
                        }
                        if (integral_value == 0.0) continue;

                        // Update STF with contact contribution
                        STF[idx].xx -= F.x * rpAB.x * integral_value;
                        STF[idx].xy -= F.x * rpAB.y * integral_value;
                        STF[idx].xz -= F.x * rpAB.z * integral_value;
                        STF[idx].yx -= F.y * rpAB.x * integral_value;
                        STF[idx].yy -= F.y * rpAB.y * integral_value;
                        STF[idx].yz -= F.y * rpAB.z * integral_value;
                        STF[idx].zx -= F.z * rpAB.x * integral_value;
                        STF[idx].zy -= F.z * rpAB.y * integral_value;
                        STF[idx].zz -= F.z * rpAB.z * integral_value;
                    }
                }
            }
        }
    }
    if (T == 0) {
        printf("Warning: Contact file '%s' is empty. No collision contributions will be added.\n", contactFile);
    }
    fclose(cFile);

    /* =========== SECTION 7.5: Compute Viscosity Field =========== */
    /* For each grid node, compute:
         1. La pression p = -1/3 (STF_xx + STF_yy + STF_zz)
         2. Le tenseur de cisaillement: tau = STF + p*I
         3. Le tenseur de taux de déformation D par différences finies sur le champ de vitesse
         4. La norme de tau et de D, puis μ = norm_tau / norm_D (si norm_D > EPSILON)
    */
    {
        int i, j, k;
        for (k = 0; k <= Cz; k++) {
            for (j = 0; j <= Cy; j++) {
                for (i = 0; i <= Cx; i++) {
                    int idx = i + j * (Cx + 1) + k * (Cx + 1) * (Cy + 1);
                    if (DensityField[idx] < EPSILON) {
                        ViscosityField[idx] = 0.0;
                        continue;
                    }

                    // 1. Compute pressure from STF (assume STF = sigma)
                    double p = -(STF[idx].xx + STF[idx].yy + STF[idx].zz) / 3.0;

                    // 2. Compute shear stress tensor: tau = STF + p*I
                    StressTensor tau;
                    tau.xx = STF[idx].xx + p;
                    tau.yy = STF[idx].yy + p;
                    tau.zz = STF[idx].zz + p;
                    tau.xy = STF[idx].xy;
                    tau.yx = STF[idx].yx;
                    tau.xz = STF[idx].xz;
                    tau.zx = STF[idx].zx;
                    tau.yz = STF[idx].yz;
                    tau.zy = STF[idx].zy;

                    // Now, explicitly compute the symmetric part of tau
                    double tau_sym_xx = tau.xx;
                    double tau_sym_yy = tau.yy;
                    double tau_sym_zz = tau.zz;
                    double tau_sym_xy = 0.5 * (tau.xy + tau.yx);
                    double tau_sym_xz = 0.5 * (tau.xz + tau.zx);
                    double tau_sym_yz = 0.5 * (tau.yz + tau.zy);


                    // Corrected norm of tau (Frobenius norm of the symmetric part)
                    double norm_tau = sqrt( tau_sym_xx*tau_sym_xx + tau_sym_yy*tau_sym_yy + tau_sym_zz*tau_sym_zz + 2.0*(tau_sym_xy*tau_sym_xy +
                            tau_sym_xz*tau_sym_xz +
                            tau_sym_yz*tau_sym_yz) );

                    // 3. Compute deformation tensor D using finite differences on VelocityField.
                    // Determine grid indices from idx:
                    int ii = i, jj = j, kk = k;
                    double dVx_dx, dVy_dy, dVz_dz;
                    double dVx_dy, dVy_dx, dVx_dz, dVz_dx, dVy_dz, dVz_dy;

                    // For derivative in x-direction (Vx):
                    if (ii == 0)
                        dVx_dx = (VelocityField[idx+1].x - VelocityField[idx].x) / step;
                    else if (ii == Cx)
                        dVx_dx = (VelocityField[idx].x - VelocityField[idx-1].x) / step;
                    else
                        dVx_dx = (VelocityField[idx+1].x - VelocityField[idx-1].x) / (2.0 * step);

                    // For derivative in y-direction (Vy):
                    int stride_y = Cx + 1;
                    if (jj == 0)
                        dVy_dy = (VelocityField[idx + stride_y].y - VelocityField[idx].y) / step;
                    else if (jj == Cy)
                        dVy_dy = (VelocityField[idx].y - VelocityField[idx - stride_y].y) / step;
                    else
                        dVy_dy = (VelocityField[idx + stride_y].y - VelocityField[idx - stride_y].y) / (2.0 * step);

                    // For derivative in z-direction (Vz):
                    int stride_z = (Cx + 1) * (Cy + 1);
                    if (kk == 0)
                        dVz_dz = (VelocityField[idx + stride_z].z - VelocityField[idx].z) / step;
                    else if (kk == Cz)
                        dVz_dz = (VelocityField[idx].z - VelocityField[idx - stride_z].z) / step;
                    else
                        dVz_dz = (VelocityField[idx + stride_z].z - VelocityField[idx - stride_z].z) / (2.0 * step);

                    // Mixed derivatives:
                    // dVx/dy
                    if (jj == 0)
                        dVx_dy = (VelocityField[idx + stride_y].x - VelocityField[idx].x) / step;
                    else if (jj == Cy)
                        dVx_dy = (VelocityField[idx].x - VelocityField[idx - stride_y].x) / step;
                    else
                        dVx_dy = (VelocityField[idx + stride_y].x - VelocityField[idx - stride_y].x) / (2.0 * step);

                    // dVy/dx
                    if (ii == 0)
                        dVy_dx = (VelocityField[idx+1].y - VelocityField[idx].y) / step;
                    else if (ii == Cx)
                        dVy_dx = (VelocityField[idx].y - VelocityField[idx-1].y) / step;
                    else
                        dVy_dx = (VelocityField[idx+1].y - VelocityField[idx-1].y) / (2.0 * step);

                    // dVx/dz
                    if (kk == 0)
                        dVx_dz = (VelocityField[idx + stride_z].x - VelocityField[idx].x) / step;
                    else if (kk == Cz)
                        dVx_dz = (VelocityField[idx].x - VelocityField[idx - stride_z].x) / step;
                    else
                        dVx_dz = (VelocityField[idx + stride_z].x - VelocityField[idx - stride_z].x) / (2.0 * step);

                    // dVz/dx
                    if (ii == 0)
                        dVz_dx = (VelocityField[idx+1].z - VelocityField[idx].z) / step;
                    else if (ii == Cx)
                        dVz_dx = (VelocityField[idx].z - VelocityField[idx-1].z) / step;
                    else
                        dVz_dx = (VelocityField[idx+1].z - VelocityField[idx-1].z) / (2.0 * step);

                    // dVy/dz
                    if (kk == 0)
                        dVy_dz = (VelocityField[idx + stride_z].y - VelocityField[idx].y) / step;
                    else if (kk == Cz)
                        dVy_dz = (VelocityField[idx].y - VelocityField[idx - stride_z].y) / step;
                    else
                        dVy_dz = (VelocityField[idx + stride_z].y - VelocityField[idx - stride_z].y) / (2.0 * step);

                    // dVz/dy
                    if (jj == 0)
                        dVz_dy = (VelocityField[idx + stride_y].z - VelocityField[idx].z) / step;
                    else if (jj == Cy)
                        dVz_dy = (VelocityField[idx].z - VelocityField[idx - stride_y].z) / step;
                    else
                        dVz_dy = (VelocityField[idx + stride_y].z - VelocityField[idx - stride_y].z) / (2.0 * step);

                    // 4. Compute components of deformation tensor D
                    double Dxx = dVx_dx;
                    double Dyy = dVy_dy;
                    double Dzz = dVz_dz;
                    double Dxy = 0.5 * (dVx_dy + dVy_dx);
                    double Dxz = 0.5 * (dVx_dz + dVz_dx);
                    double Dyz = 0.5 * (dVy_dz + dVz_dy);

                    // Norm of D: norm_D = sqrt(2*(Dxx^2 + Dyy^2 + Dzz^2 + Dxy^2 + Dxz^2 + Dyz^2))
                    double norm_D = sqrt(2.0 * (Dxx*Dxx + Dyy*Dyy + Dzz*Dzz + Dxy*Dxy + Dxz*Dxz + Dyz*Dyz));

                    // 5. Compute viscosity: if norm_D > EPSILON, μ = norm_tau/norm_D, else 0.
                    double mu = (norm_D > EPSILON) ? norm_tau / norm_D : 0.0;
                    ViscosityField[idx] = mu;
                }
            }
        }
    }

    // =========== SECTION 8: Write Data to CSV ===========
    FILE *outputFile = fopen(outFile, "w");
    if (!outputFile) {
        perror("Error opening output CSV file");
        freeAllMatrices(NUM_FILES, MATRIX_SIZE);
        free(Grid);
        free(Momentum); free(DensityField);
        free(VelocityField); free(STF); free(ParticlePos); free(ViscosityField);
        return 1;
    }

    /* Updated CSV header to include Viscosity */
    fprintf(outputFile, "NodeID,X,Y,Z,Density,Vx,Vy,Vz,STF_xx,STF_xy,STF_xz,"
                        "STF_yx,STF_yy,STF_yz,STF_zx,STF_zy,STF_zz,Viscosity\n");
    for (int i = 0; i < N; i++) {
        fprintf(outputFile,
            "%d,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,"
            "%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f\n",
            i,
            Grid[i].x, Grid[i].y, Grid[i].z,
            DensityField[i],
            VelocityField[i].x, VelocityField[i].y, VelocityField[i].z,
            STF[i].xx, STF[i].xy, STF[i].xz,
            STF[i].yx, STF[i].yy, STF[i].yz,
            STF[i].zx, STF[i].zy, STF[i].zz,
            ViscosityField[i]);
    }
    fclose(outputFile);

    printf("   [Set] Generated: %s\n", outFile);

    // =========== SECTION 10: Free memory ===========
    free(Grid);
    free(Momentum); free(DensityField);
    free(VelocityField); free(STF); free(ParticlePos); free(ViscosityField);

    // ----------------- End timing for this set ------------------
    clock_t local_end = clock();
    double local_elapsed = (double)(local_end - local_start) / CLOCKS_PER_SEC;
    printf("   [Set] Elapsed time (this set): %.4f seconds\n", local_elapsed);

    return 0;
}

/*****************************************************************************
 * MAIN: Processes all DEM files in "DEM_data" for indices from 0..9999
 *       Creates "output_data" directory if it doesn't exist.
 *       If "DEM_data" directory is missing, prompts the user and exits.
 *       Measures total elapsed time across all sets.
 *****************************************************************************/
int main(int argc, char *argv[])
{
    // ========== SECTION 0: Setup and Argument Parsing ==========
    // Default input and output directories
    const char *demDirectory = DEFAULT_INPUT_DIR;
    const char *outputDirectory = DEFAULT_OUTPUT_DIR;

    // Check for Rc argument
    if (argc < 2) {
        fprintf(stderr,
                "Usage: %s <Rc>\n"
                "  <Rc> = cutoff radius\n",
                argv[0]);
        return 1;
    }

    // Parse Rc
    double Rc = atof(argv[1]);
    if (Rc <= 0.0) {
        fprintf(stderr, "Error: Rc must be a positive number.\n");
        return 1;
    }

    // ========== SECTION 1: Check if DEM_data directory exists ==========
    struct stat st;
    if (stat(demDirectory, &st) != 0 || !S_ISDIR(st.st_mode)) {
        fprintf(stderr, "Error: Input directory '%s' does not exist.\n", demDirectory);
        fprintf(stderr, "Please create the '%s' directory and add the necessary files.\n", demDirectory);
        return 1;
    }

    // ========== SECTION 2: Create output_data directory if it doesn't exist ==========
    if (stat(outputDirectory, &st) != 0) {
        // Directory does not exist, attempt to create it
        if (mkdir(outputDirectory, 0755) != 0) {
            perror("Error creating output_data directory");
            return 1;
        }
        printf("Created output directory: %s\n", outputDirectory);
    }
    else if (!S_ISDIR(st.st_mode)) {
        fprintf(stderr, "Error: '%s' exists but is not a directory.\n", outputDirectory);
        return 1;
    }

    // ========== SECTION 5: Load Integral Matrices ===========
    clock_t startInt = clock();
    if (loadAllMatrices(NUM_FILES, MATRIX_SIZE) != 0) {
        fprintf(stderr, "Error loading integral matrices.\n");
        return 1;
    }
    clock_t endInt = clock();
    double cpu_time_integrals = ((double)(endInt - startInt)) / CLOCKS_PER_SEC;
    printf(" Elapsed time for loading integrals: %.4f s\n", cpu_time_integrals);

    // ========== SECTION 3: Start total timer ==========
    clock_t total_start = clock();

    // ========== SECTION 4: Loop over possible indices ==========
    for (int i = 0; i <= MAX_INDEX; i++) {
        // Build the 4-digit string, e.g. "0000", "0001", ...
        char indexStr[5];
        snprintf(indexStr, sizeof(indexStr), "%04d", i);

        // Build the expected filenames
        char headerFile[512];
        char demFile[512];
        char contactFile[512];
        char outFile[512];

        snprintf(headerFile,  512, "%s/header0000.txt", demDirectory);
        snprintf(demFile,     512, "%s/DEMdemo_output_%s.csv", demDirectory, indexStr);
        snprintf(contactFile, 512, "%s/Contact_pairs_%s.csv",  demDirectory, indexStr);
        snprintf(outFile,     512, "%s/output_data_%s.csv",    outputDirectory, indexStr);

        // Check if all three input files exist
        if (access(headerFile,  F_OK) == 0 &&
            access(demFile,     F_OK) == 0 &&
            access(contactFile, F_OK) == 0)
        {
            // We found a valid "set" => process it
            printf("-------------------------------------------------------\n");
            printf("[Index %s] Processing files:\n", indexStr);
            printf("   header   = %s\n", headerFile);
            printf("   demFile  = %s\n", demFile);
            printf("   contact  = %s\n", contactFile);

            int ret = processSingleSet(headerFile, demFile, contactFile, outFile, Rc);
            if (ret != 0) {
                fprintf(stderr, "[Index %s] Error in processSingleSet()\n", indexStr);
            }
        }
        // else: skip any index that doesn't have a complete set
    }

    freeAllMatrices(NUM_FILES, MATRIX_SIZE);

    // ========== SECTION 5: Stop total timer ==========
    clock_t total_end = clock();
    double total_elapsed = (double)(total_end - total_start) / CLOCKS_PER_SEC;
    printf("-------------------------------------------------------\n");
    printf("All possible data sets processed.\n");
    printf("Total elapsed time: %.4f seconds\n", total_elapsed);

    return 0;
}

/*****************************************************************************
 *                          Function Definitions
 *****************************************************************************/

/**
 * Calculate the Euclidean distance between two 3D points
 *
 * @param p1  First 3D point
 * @param p2  Second 3D point
 * @return    Euclidean distance between p1 and p2
 */
double distance3d(Point3d p1, Point3d p2) {
    double dx = p2.x - p1.x;
    double dy = p2.y - p1.y;
    double dz = p2.z - p1.z;
    return sqrt(dx*dx + dy*dy + dz*dz);
}

/**
 * Subtract two 3D vectors (points)
 *
 * @param pA  First 3D point
 * @param pB  Second 3D point
 * @return    Vector from pB to pA (pA - pB)
 */
Point3d vectorSubtract(Point3d pA, Point3d pB) {
    Point3d result;
    result.x = pA.x - pB.x;
    result.y = pA.y - pB.y;
    result.z = pA.z - pB.z;
    return result;
}

/**
 * Check if a 3D point is the null vector (0,0,0)
 *
 * @param p  3D point to check
 * @return   1 if point is (0,0,0), 0 otherwise
 */
int isNull(Point3d p){
    if (p.x == 0.0 && p.y == 0.0 && p.z == 0.0){
        return 1;
    }
    return 0;
}

/**
 * Calculate the dot product of two 3D vectors
 *
 * @param pA  First 3D vector
 * @param pB  Second 3D vector
 * @return    Dot product (pA·pB)
 */
double dotProduct(Point3d pA, Point3d pB) {
    return pA.x * pB.x + pA.y * pB.y + pA.z * pB.z;
}

/**
 * Calculate the mass of a spherical particle
 *
 * @param radius   Radius of the particle
 * @param density  Density of the particle material
 * @return         Mass of the particle (4/3 * π * r³ * ρ)
 */
double computeMass(double radius, double density) {
    return (4.0 / 3.0) * M_PI * pow(radius, 3) * density;
}

/**
 * Lucy's weight function for coarse-graining approach
 *
 * This function implements the weight function used to distribute
 * particle properties to the grid points in the coarse-graining approach.
 *
 * @param q     Normalized distance (q = 2*r/RC where r is distance and RC is cutoff radius)
 * @param s_3   Normalization factor (s_3 = 1/(π*(RC/2)³))
 * @return      Weight value at distance q
 */
double phi(double q, double s_3) {
    if (q >= 0 && q <= 1) {
        // Polynomial form for 0 ≤ q ≤ 1
        return s_3 * (1 - 1.5 * q * q * (1 - q / 2.0));
    } else if (q > 1 && q <= 2) {
        // Polynomial form for 1 < q ≤ 2
        return (s_3 / 4) * pow(2 - q, 3);
    }
    // Zero outside the support (q > 2)
    return 0.0;
}

/****************************************************************************
 * MATRIX UTILITIES
 ****************************************************************************/

/**
 * Allocate a square matrix of the specified size
 *
 * This function allocates memory for a square matrix and initializes all elements to zero.
 *
 * @param size  Size of the square matrix to allocate
 * @return      Pointer to the allocated matrix or NULL if allocation fails
 *
 * @author Marouane ELBISSSOURI
 * @date 2024-10-01
 */
double** allocateMatrix(int size) {
    // Allocate array of row pointers
    double **matrix = (double **)malloc(size * sizeof(double *));
    if (!matrix) {
        fprintf(stderr, "Memory allocation failed for matrix rows.\n");
        return NULL;
    }

    // Allocate each row
    for (int i = 0; i < size; i++) {
        matrix[i] = (double *)malloc(size * sizeof(double));
        if (!matrix[i]) {
            // Clean up previously allocated memory if allocation fails
            fprintf(stderr, "Memory allocation failed for row %d.\n", i);
            for (int j = 0; j < i; j++) {
                free(matrix[j]);
            }
            free(matrix);
            return NULL;
        }
        // Initialize all elements to zero
        for (int j = 0; j < size; j++) {
            matrix[i][j] = 0.0;
        }
    }
    return matrix;
}

/**
 * Free memory allocated for a square matrix
 *
 * This function properly deallocates all memory used by a matrix
 * created with allocateMatrix().
 *
 * @param matrix  Pointer to the matrix to free
 * @param size    Size of the matrix
 *
 * @author Marouane ELBISSSOURI
 * @date 2024-10-01
 */
void freeMatrix(double **matrix, int size) {
    // Check for NULL pointer
    if (!matrix) return;

    // Free each row
    for (int i = 0; i < size; i++) {
        free(matrix[i]);
    }

    // Free the array of row pointers
    free(matrix);
}

/****************************************************************************
 * CSV PARSING UTILITIES
 ****************************************************************************/

/**
 * Split a CSV line into tokens
 *
 * This function parses a CSV line and splits it into individual tokens,
 * handling empty fields correctly. It allocates memory for each non-empty token,
 * which must be freed by the caller. Empty fields are represented by NULL pointers.
 *
 * @param line        The CSV line to split
 * @param tokens      Array to store pointers to the extracted tokens
 * @param max_tokens  Maximum number of tokens to extract
 * @return            Number of tokens extracted, or -1 on error
 *
 * @author Marouane ELBISSSOURI
 * @date 2024-10-01
 */
int splitCSV(const char *line, char **tokens, int max_tokens) {
    int count = 0;
    const char *ptr = line;
    const char *start = ptr;

    // Process each character in the line
    while (*ptr != '\0' && count < max_tokens) {
        if (*ptr == ',') {
            // Found a field separator (comma)
            int len = ptr - start;

            // Remove trailing carriage return if present
            if (len > 0 && start[len - 1] == '\r') {
                len--;
            }

            if (len == 0) {
                // Handle empty field
                tokens[count++] = NULL;
            } else {
                // Allocate memory for the token
                char *token = (char *)malloc(len + 1);
                if (!token) {
                    // Clean up previously allocated tokens on error
                    for (int t = 0; t < count; t++) {
                        if (tokens[t]) free(tokens[t]);
                    }
                    return -1;
                }

                // Copy the token and null-terminate it
                strncpy(token, start, len);
                token[len] = '\0';
                tokens[count++] = token;
            }

            // Move to the next field
            ptr++;
            start = ptr;
        } else {
            ptr++;
        }
    }

    // Handle the last field (after the final comma)
    if (count < max_tokens) {
        int len = ptr - start;

        // Remove trailing newline and carriage return if present
        if (len > 0 && start[len - 1] == '\n') len--;
        if (len > 0 && start[len - 1] == '\r') len--;

        if (len == 0) {
            // Handle empty last field
            tokens[count++] = NULL;
        } else {
            // Allocate memory for the last token
            char *token = (char *)malloc(len + 1);
            if (!token) {
                // Clean up on error
                for (int t = 0; t < count; t++) {
                    if (tokens[t]) free(tokens[t]);
                }
                return -1;
            }

            // Copy the token and null-terminate it
            strncpy(token, start, len);
            token[len] = '\0';
            tokens[count++] = token;
        }
    }

    // Fill remaining slots with NULL pointers
    while (count < max_tokens) {
        tokens[count++] = NULL;
    }

    return count;
}

/**
 * Load a matrix from a CSV file
 *
 * This function reads a CSV file containing integral values and loads them
 * into a dynamically allocated matrix. The CSV file is expected to have
 * a header row and the first column contains x values.
 *
 * @param filename  Path to the CSV file to load
 * @param size      Size of the square matrix to create
 * @return          Pointer to the loaded matrix or NULL on error
 *
 * @author Marouane ELBISSSOURI
 * @date 2024-10-01
 */
double** loadMatrix(const char *filename, int size) {
    // Open the CSV file
    FILE *file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error opening: %s\n", filename);
        return NULL;
    }

    // Allocate memory for the matrix
    double **matrix = allocateMatrix(size);
    if (!matrix) {
        fclose(file);
        return NULL;
    }

    // Buffer for reading lines from the file
    char line[MAX_LINE_LENGTH];

    // Skip the header row (contains column labels)
    if (!fgets(line, sizeof(line), file)) {
        fprintf(stderr, "Error reading header from %s\n", filename);
        freeMatrix(matrix, size);
        fclose(file);
        return NULL;
    }

    // Prepare array for storing tokens from each CSV line
    int max_tokens = size + 1;  // +1 for the first column (x values)
    char **tokens = (char **)malloc(max_tokens * sizeof(char *));
    if (!tokens) {
        fprintf(stderr, "Memory allocation failed for tokens.\n");
        freeMatrix(matrix, size);
        fclose(file);
        return NULL;
    }

    // Process each row of the CSV file
    for (int i = 0; i < size; i++) {
        // Read the next line
        if (!fgets(line, sizeof(line), file)) {
            fprintf(stderr, "Error reading line %d in %s\n", i + 1, filename);
            free(tokens);
            freeMatrix(matrix, size);
            fclose(file);
            return NULL;
        }

        // Split the line into tokens
        int num_tokens = splitCSV(line, tokens, max_tokens);
        if (num_tokens < 0) {
            fprintf(stderr, "Error parsing CSV line in %s\n", filename);
            free(tokens);
            freeMatrix(matrix, size);
            fclose(file);
            return NULL;
        }

        // Skip the first column (x values) and any values in the lower triangle
        int current_token = 1;  // Start from the second column (first data column)
        for (int skip = 0; skip < i; skip++) {
            current_token++;
            if (current_token >= num_tokens) break;
        }

        // Process the upper triangle of the matrix
        for (int j = i; j < size; j++) {
            if (current_token >= num_tokens) break;

            // Convert token to double if it exists, otherwise use 0.0
            if (tokens[current_token]) {
                matrix[i][j] = atof(tokens[current_token]);
            } else {
                matrix[i][j] = 0.0;
            }
            current_token++;
        }

        // Free memory allocated for tokens
        for (int t = 0; t < num_tokens; t++) {
            if (tokens[t]) free(tokens[t]);
        }
    }

    // Clean up resources
    free(tokens);
    fclose(file);

    return matrix;
}

/**
 * Load all integral matrices from CSV files
 *
 * This function loads all the integral matrices needed for the homogenization
 * process. It creates one matrix for each theta value and stores them in the
 * global integral_matrices array.
 *
 * @param num_files  Number of files to load (one per theta value)
 * @param size       Size of each matrix
 * @return           0 on success, -1 on error
 *
 * @author Marouane ELBISSSOURI
 * @date 2024-10-01
 */
int loadAllMatrices(int num_files, int size) {
    // Allocate memory for the array of matrices
    integral_matrices = (double ***)malloc(num_files * sizeof(double **));
    if (!integral_matrices) {
        fprintf(stderr, "Error: integral_matrices allocation.\n");
        return -1;
    }

    // Initialize progress tracking variables
    long total_files = num_files;
    long current_file = 0;
    int last_percent = -1;

    printf("=== Loading Integral Matrices ===\n");

    // Load each matrix file
    for (int k = 0; k < num_files; k++) {
        // Construct the filename for this theta value
        char filename[128];
        snprintf(filename, sizeof(filename), FILENAME_FORMAT, k);

        // Load the matrix from the file
        double **mat = loadMatrix(filename, size);
        if (!mat) {
            fprintf(stderr, "Failed loading matrix %d (%s)\n", k, filename);

            // Clean up previously loaded matrices
            for (int m = 0; m < k; m++) {
                freeMatrix(integral_matrices[m], size);
            }
            free(integral_matrices);
            integral_matrices = NULL;
            return -1;
        }

        // Store the loaded matrix
        integral_matrices[k] = mat;

        // Update and display progress
        current_file++;
        int percent = (int)((double)current_file / total_files * 100);
        if (percent != last_percent) {
            display_progress((double)current_file / total_files);
            last_percent = percent;
        }
    }

    printf("\nDone loading integrals.\n");
    return 0;
}

/**
 * Free all integral matrices
 *
 * This function properly deallocates all memory used by the integral matrices
 * that were loaded with loadAllMatrices().
 *
 * @param num_files  Number of matrices to free
 * @param size       Size of each matrix
 *
 * @author Marouane ELBISSSOURI
 * @date 2024-10-01
 */
void freeAllMatrices(int num_files, int size) {
    // Check if matrices were allocated
    if (!integral_matrices) return;

    // Free each matrix
    for (int k = 0; k < num_files; k++) {
        freeMatrix(integral_matrices[k], size);
    }

    // Free the array of matrices and reset the pointer
    free(integral_matrices);
    integral_matrices = NULL;
}

/**
 * Display a progress bar in the console
 *
 * This function displays a visual progress bar to show the status
 * of a long-running operation.
 *
 * @param progress  Progress value between 0.0 and 1.0
 *
 * @author Marouane ELBISSSOURI
 * @date 2024-10-01
 */
void display_progress(double progress) {
    int barWidth = 50;  // Width of the progress bar in characters

    // Print the start of the progress bar
    printf("[");

    // Calculate the position of the progress indicator
    int pos = (int)(barWidth * progress);

    // Print the progress bar characters
    for (int i = 0; i < barWidth; i++) {
        if (i < pos) printf("=");         // Completed portion
        else if (i == pos) printf(">");   // Current position indicator
        else printf(" ");                 // Remaining portion
    }

    // Print the percentage and carriage return (\r) to update in-place
    printf("] %3d%%\r", (int)(progress * 100));
    fflush(stdout);  // Ensure output is displayed immediately
}
