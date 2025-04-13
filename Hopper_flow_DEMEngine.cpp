//  Copyright (c) 2021, SBEL GPU Development Team
//  Copyright (c) 2021, University of Wisconsin - Madison
//
//        SPDX-License-Identifier: BSD-3-Clause

// =============================================================================
// Hopper Flow DEM Simulation
// =============================================================================
// This simulation models granular material flowing through a hopper/funnel system.
// The simulation consists of several phases:
// 1. Setup of the hopper geometry using mesh objects
// 2. Creation of spherical particles with defined properties
// 3. Settling phase where particles are generated and allowed to settle
// 4. Flow phase where the gate opens and particles flow through the hopper
//
// The simulation outputs position, velocity, and contact data at specified intervals
// for post-processing and analysis using the HomogenizationDEM framework.
// =============================================================================

#include <core/ApiVersion.h>
#include <core/utils/ThreadManager.h>
#include <DEM/API.h>
#include <DEM/HostSideHelpers.hpp>
#include <DEM/utils/Samplers.hpp>

#include <cstdio>
#include <chrono>
#include <filesystem>
#include <random>


#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

void modify(const std::string& filename, double X, double Y, double Z);

using namespace deme;
using namespace std::filesystem;

int main(){

        /* Part 0: Defining the Solver and Materials */
        // This section initializes the DEM solver and defines material properties for the simulation
        // The solver uses the Hertzian contact model which accounts for elastic deformation and friction

        // Initialize the DEME solver and set up the contact model
        DEMSolver DEMSim;
        DEMSim.UseFrictionalHertzianModel();  // Use Hertzian model for frictional contact - more accurate for elastic materials

        // Configure solver verbosity and output settings
        DEMSim.SetVerbosity(INFO);
        DEMSim.SetOutputFormat(OUTPUT_FORMAT::CSV); // Set output format to CSV
        // Output XYZ and velocity data (it controls the contact of the data that will be extracted in the Part4: Flow phase)
        DEMSim.SetOutputContent(OUTPUT_CONTENT:YZ | OUTPUT_CONTENT::VEL);
        // Contact data output (Same as the previous but for contact data)
        DEMSim.SetContactOutputContent(CNT_TYPE | OWNER | FORCE | DEME_POINT);
        DEMSim.EnsureKernelErrMsgLineNum(); // Ensure error messages show line numbers

        // Set additional solver options for collecting data and error handling
        DEMSim.SetCollectAccRightAfterForceCalc(true);
        DEMSim.SetErrorOutAvgContacts(80); // Output error if average contacts exceed threshold

        // Define materials and their properties
        auto mat_type_bottom = DEMSim.LoadMaterial({{"E", 10e9}, {"nu", 0.3}, {"CoR", 0.60}});
        auto mat_type_flume = DEMSim.LoadMaterial({{"E", 10e9}, {"nu", 0.3}, {"CoR", 0.60}});
        auto mat_type_walls = DEMSim.LoadMaterial({{"E", 10e9}, {"nu", 0.3}, {"CoR", 0.60}});
        auto mat_spheres = DEMSim.LoadMaterial({{"E", 1.0e9}, {"nu", 0.35}, {"CoR", 0.85}, {"mu", 0.40}, {"Crr", 0.04}});

        // Set material interaction properties
        DEMSim.SetMaterialPropertyPair("mu", mat_type_walls, mat_spheres, 0.30);
        DEMSim.SetMaterialPropertyPair("CoR", mat_type_walls, mat_spheres, 0.7);
        DEMSim.SetMaterialPropertyPair("Crr", mat_type_walls, mat_spheres, 0.05);
        DEMSim.SetMaterialPropertyPair("mu", mat_type_flume, mat_spheres, 0.30);
        DEMSim.SetMaterialPropertyPair("CoR", mat_type_flume, mat_spheres, 0.70);
        DEMSim.SetMaterialPropertyPair("Crr", mat_type_flume, mat_spheres, 0.05);

        /* End of Part 0 */

        /*Part1: Create the hopper*/
        // This section defines the geometry of the hopper system
        // The hopper consists of two angled walls (funnel) and a movable gate at the bottom
        // The geometry is created using mesh objects loaded from OBJ files
        // The dimensions and angles can be adjusted to control the flow characteristics

        //Domain Dimensions for the simulation box:
        double X = 0.4;
        double Y = 0.4;
        double Z = 2;

        double x = X/2;
        double y = Y/2;
        double z0 = Z/2;

        DEMSim.InstructBoxDomainDimension({-x, x}, {-y, y}, {-z0, z0});
        DEMSim.InstructBoxDomainBoundingBC("top_open", mat_type_walls);
        DEMSim.SetGravitationalAcceleration(make_float3(0, 0, -9.81));


        // Define the full paths to the mesh files
        std::string baseDir = "/export/home/elbissouri/Marouane/DEM-Engine/data/mesh/Hopper/";
        std::string file1 = baseDir + "funnel_left.obj";
        std::string file2 = baseDir + "gate.obj";

        double gateWidth = 0.05;
        double angle = deme:I/3;
        double e = 0.01;  //Thickness of the gate and funnel

        double D = gateWidth + 2 * std::sin(angle) * e;
        double L = ((X + 2 * e - D) / 2) / std::cos(angle);


        // Modify the files with the given X, Y, Z values
        modify(file1, L, Y + 2*e, e);
        modify(file2, D, Y + 2*e, e);


        // Loaded meshes are by-default fixed
        int Nmesh = 0 ; // Number fo meshes

        float3 mov0 = make_float3(0,0,0);
        float4 rot0 = make_float4(0,0,0,0);
        float4 rot1 = make_float4(0.7071, 0, 0, 0.7071);

        auto fixed_left = DEMSim.AddWavefrontMeshObject(file1, mat_type_flume);

        fixed_left->Move(mov0, rot1);

        float3 move = make_float3(-(X + 2 * e + D) / 4 + (e / 2) * sin(angle), 0 , (e / 2) * cos(angle) - (L / 2) * sin(angle) - e);
        float4 rot2 = make_float4(std::cos(angle/2),0,-std::sin(angle/2),0);

        fixed_left->Move(move, rot2);

        Nmesh ++;


        auto fixed_right = DEMSim.AddWavefrontMeshObject(file1, mat_type_flume);

        fixed_right->Move(mov0, rot1);

        move = make_float3((X + 2 * e + D) / 4 - (e / 2) * sin(angle), 0 , (e / 2) * cos(angle) - (L / 2) * sin(angle) - e);
        rot2 = make_float4(std::cos(angle/2),0,std::sin(angle/2),0);

        fixed_right->Move(move, rot2);

        Nmesh ++;


        auto gate = DEMSim.AddWavefrontMeshObject(file2, mat_type_flume);

        double l = L*sin(angle);
        gate->Move(make_float3(0, 0, -l +  e*(cos(angle)/2 - 1)+0.0001), rot1);

        Nmesh ++;

        fixed_left->SetFamily(10);
        fixed_right->SetFamily(10);
        gate->SetFamily(3);

        std::string shake_pattern_xz = " 0.0 * sin( 300 * 2 * deme:I * t)";
        std::string shake_pattern_y = " 0.0 * sin( 30 * 2 * deme:I * t)";

        DEMSim.SetFamilyFixed(1);
        DEMSim.SetFamilyFixed(3);

        double gateSpeed = -3.5;
        DEMSim.SetFamilyPrescribedLinVel(4, "0", "0", to_string_with_precision(gateSpeed));
        DEMSim.SetFamilyPrescribedLinVel(10, shake_pattern_xz, shake_pattern_y, shake_pattern_xz);

        /*End of Part1*/

        /*Part2: Creating sphere templates */
        // This section defines the properties of the particles used in the simulation
        // Particles are created as spheres (or clumps of spheres for more complex shapes)
        // Properties include radius, density, and material characteristics
        // These templates are used to generate the actual particles during the simulation

        // total number of random clump templates to generate
        // In this case, we're using just one template for simplicity
        int num_template = 1;

        // data for the spheres
        double radiusMax = 0.003;
        double radiusMin = 0.003;
        double densitySph = 15000; // (kg/m³)


        std::mt19937 generator(123); // Seeded with a constant value for reproducibility

        std::uniform_real_distribution<double> distribution(radiusMin, radiusMax);

        // Make an array to store these generated clump templates
        std::vector<std::shared_ptr<DEMClumpTemplate>> clump_sphere;

        {  // initialize spheres
        for (int i = 0; i < num_template; i++) {
            std::vector<float> radii;
            std::vector<float3> relPos;
            std::vector<std::shared_ptr<DEMMaterial>> mat;
            auto tmp = make_float3(0, 0, 0);
            //double radius = distribution(generator);
        double radius = radiusMax;
            // double radiusMax = radiusSph;

            relPos.push_back(tmp);
            mat.push_back(mat_spheres);
            radii.push_back(radius);

            double c = radius;
            double b = radius;
            double a = radius;

            float mass = 4.0 / 3.0 * PI * a * b * c * densitySph;
            float3 MOI = make_float3(1.f / 5.f * mass * (b * b + c * c), 1.f / 5.f * mass * (a * a + c * c),
                                     1.f / 5.f * mass * (b * b + a * a));
            //std::cout << a << " chosen moi ..." << a / radius << std::endl;

            auto clump_ptr = DEMSim.LoadClumpType(mass, MOI, radii, relPos, mat_spheres);
            // clump_ptr->AssignName("fsfs");
            clump_sphere.push_back(clump_ptr);
        }
        }


        /*End of Part2*/



        // Define output directory path and ensure it is clean
        path out_dir = current_path();
        out_dir += "/Test_Plastic_Sphere_Cylinder/";
        out_dir += "Hopper/";

        // Remove old data and create new directories
        remove_all(out_dir);
        create_directories(out_dir);



        // Prepare simulation parameters
        float step_size = 5*1.0e-7; // Define time step size

        DEMSim.SetInitTimeStep(step_size); // Set initial time step for simulation
        DEMSim.SetMaxVelocity(25.); // Set maximum allowed velocity for particles
        DEMSim.SetInitBinSize(radiusMax * 2); // Set initial bin size for spatial sorting

        DEMSim.Initialize(); // Initialize the simulation with the defined parameters

        // Initialize frame count and file name for mesh output
        int frame = 0;
        char meshfile[200];

        // Generate initial VTK file for visualization
        sprintf(meshfile, "%s/DEMdemo_funnel_%04d.vtk", out_dir.c_str(), frame);
        DEMSim.WriteMeshFile(std::string(meshfile));


        /*Part3: Settling phase setup*/
        // This section handles the generation and initial settling of particles
        // Particles are created in batches and allowed to settle under gravity
        // This creates a realistic initial configuration before opening the gate
        // The process continues until the desired number of particles is reached
        // or until the hopper is filled to the desired level

        int totalSph = 500000; // Total number of spheres to be generated
        float z = 0;
        float timeTotal = 0.0; // Total elapsed simulation time
        float settle_frame_time = 0.005; // Time step for each simulation frame

        // Inspectors to monitor particle properties during the simulation
        auto max_z_finder = DEMSim.CreateInspector("clump_max_z");
        auto min_z_finder = DEMSim.CreateInspector("clump_min_z");
        auto total_mass_finder = DEMSim.CreateInspector("clump_mass");
        auto max_v_finder = DEMSim.CreateInspector("clump_max_absv");

        char filename[200]; // Filename buffer for output files

        // Sphere generation parameters
        float plane_bottom = 0.02f; // Initial height of the plane at the bottom of the hopper
        unsigned int actualTotalSpheres = 0; // Counter for the number of generated spheres

        {
            float shift_xyz = 1.0 * (radiusMax) * 2.0; // Shift distance for placing new spheres
            float x = 0;
            float y = 0;
            double emitterZ = 1.00; // Initial Z position for sphere generation

            bool generate = true; // Flag to control generation of spheres
            bool initialization = true; // Flag to control the initialization phase
            double consolidation = true; // Flag to manage the consolidation phase

            while (initialization) {
                std::vector<std::shared_ptr<DEMClumpTemplate>> input_pile_template_type;
                std::vector<float3> input_pile_xyz;
                PDSampler sampler(shift_xyz); // Sampler for generating sphere positions

                // Determine whether to continue generating spheres
                bool generate = (plane_bottom + shift_xyz > emitterZ) ? false : true;

                if (generate) {
                    float sizeZ = (frame == 0) ? 0.20 : 0; // Initial pile height
                    float sizeX = X - shift_xyz; // Width of the pile
                    float sizeY = Y - shift_xyz; // Depth of the pile
                    float z = plane_bottom + shift_xyz + sizeZ / 2.0; // Center height of the pile

                    float3 center_xyz = make_float3(0, 0, z); // Center position for sphere generation
                    float3 size_xyz = make_float3(sizeX / 2.0, sizeY / 2.0, sizeZ / 2.0); // Size of the generation volume

                    std::cout << "Level of particle positions: " << center_xyz.z << std::endl;

                    auto heap_particles_xyz = sampler.SampleBocenter_xyz, size_xyz); // Sample positions within a box
                    unsigned int num_clumps = heap_particles_xyz.size(); // Number of spheres generated at this level
                    std::cout << "Number of particles at this level: " << num_clumps << std::endl;

                    // Assign templates to the generated spheres
                    for (unsigned int i = actualTotalSpheres; i < actualTotalSpheres + num_clumps; i++) {
                        input_pile_template_type.push_back(clump_sphere.at(i % num_template));
                    }

                    input_pile_xyz.insert(input_pile_xyz.end(), heap_particles_xyz.begin(), heap_particles_xyz.end());

                    // Add generated spheres to the simulation
                    auto the_pile = DEMSim.AddClumps(input_pile_template_type, input_pile_xyz);
                    the_pile->SetVel(make_float3(-0.00, 0.0, -0.80)); // Set initial velocity
                    the_pile->SetFamily(99); // Assign to a specific family for control

                    DEMSim.UpdateClumps(); // Update clumps in the simulation

                    std::cout << "Total number of particles: " << (int)DEMSim.GetNumClumps() << std::endl;
                    actualTotalSpheres = (int)DEMSim.GetNumClumps(); // Update the total number of spheres
                }

                // Increment simulation time
                timeTotal += settle_frame_time;
                std::cout << "Total runtime: " << timeTotal << "s; settling for: " << settle_frame_time << std::endl;
                std::cout << "Max Z position: " << max_z_finder->GetValue() << std::endl;

                // Check if more spheres need to be generated
                initialization = (actualTotalSpheres < totalSph) ? true : false;

                // Output data to file
                if (generate) {
                    std::cout << "Frame: " << frame << std::endl;
                    sprintf(filename, "%s/DEMdemo_settling_%04d.csv", out_dir.c_str(), frame);
                    DEMSim.WriteSphereFile(std::string(filename)); // Write sphere positions to a CSV file
                    frame++;
                }

                // Perform a dynamics step in the simulation
                DEMSim.DoDynamicsThenSync(settle_frame_time);

                // Update the bottom plane position for the next generation step
                plane_bottom = max_z_finder->GetValue();

                // Consolidation phase: allow the particles to settle further after initialization
                if (!initialization) {
                    for (int i = 0; i < (int)(0.4 / settle_frame_time); i++) {
                        DEMSim.DoDynamics(settle_frame_time); // Run dynamics without synchronization
                        sprintf(filename, "%s/DEMdemo_settling_%04d.csv", out_dir.c_str(), frame++);
                        DEMSim.WriteSphereFile(std::string(filename)); // Output settling data
                        std::cout << "Consolidating for " << i * settle_frame_time << "s " << std::endl;
                    }
                }
            }
        }


        /*End of Part3: Settling phase*/

        // Perform an initial synchronization step with zero time progression
        DEMSim.DoDynamicsThenSync(0.0);

        // Reset the frame counter to start tracking frames for the flow phase
        frame = 0;

        /* Part4: Flow phase */
        // This is the main simulation phase where the gate opens and particles flow through the hopper
        // The simulation captures the dynamics of granular flow under gravity
        // Data is recorded at regular intervals for post-processing and analysis
        // This data includes particle positions, velocities, and contact forces
        // The output files are compatible with the HomogenizationDEM framework for further analysis

        // Parameters for the flow phase
        double gateOpen = z0 - l + 0.01; // The distance the gate will travel to open
        float timeStep = step_size;
        int numStep = 3.0 / timeStep; // Total number of simulation steps for the flow phase
        int timeOut = 0.1 / timeStep; // In order to capture output data each 0.01 s .
        int gateMotion = (gateOpen / gateSpeed) / timeStep; // Time steps required for the gate to complete its motion
        std::cout << "Frame: " << timeOut << std::endl;

        char cnt_filename[200]; // File name for contact data output
        bool status = true;  // Status flag to control the gate's initial movement
        bool stopGate = true; // Status flag to stop the gate after opening
        float totalRunTime = 0.0f; // Total elapsed simulation time

        baseDir = "/export/home/elbissouri/Marouane/Engine_build/bin/Test_Plastic_Sphere_Cylinder/Hopper/";
        std::string file3;

        // Main simulation loop for the flow phase
        for (int i = 0; i < numStep; i++) {
            DEMSim.DoDynamics(timeStep); // Perform a simulation step
            totalRunTime += timeStep; // Update total run time

            // Output data at specified intervals
            if (i % timeOut == 0 || i == 0)  {
                sprintf(filename, "%s/DEMdemo_output_%04d.csv", out_dir.c_str(), frame);
                sprintf(meshfile, "%s/DEMdemo_mesh_%04d.vtk", out_dir.c_str(), frame);
                sprintf(cnt_filename, "%s/Contact_pairs_%04d.csv", out_dir.c_str(), frame);

                // Write the current state to output files
                DEMSim.WriteMeshFile(std::string(meshfile));
                DEMSim.WriteSphereFile(std::string(filename));
                DEMSim.WriteContactFile(std::string(cnt_filename));

                std::cout << "Frame: " << frame << std::endl;
                std::cout << "Elapsed time: " << totalRunTime << std::endl;

                // Generate header file with simulation parameters
                file3 = baseDir + "header" + std::to_string(frame).insert(0, 4 - std::to_string(frame).length(), '0') + ".txt";
                std::ofstream file(file3);
                if (!file) {
                    std::cerr << "Error opening file!" << std::endl;
                    return 1;
                }

                // Write simulation details to the file
                file << "@ Np" << std::endl;
                file << actualTotalSpheres << std::endl;
                file << "@ RadiusMax" << std::endl;
                file << radiusMax << std::endl;
                file << "@ Density" << std::endl;
                file << densitySph << std::endl;
                file << "@ Nmesh" << std::endl;
                file << Nmesh << std::endl;
                file << "@ Niter" << std::endl;
                file << numStep << std::endl;
                file << "@ timeBetween2Snap" << std::endl;
                file << timeOut << std::endl;
                file << "@ Elapsed Time" << std::endl;
                file << totalRunTime << std::endl;
                file << "@ DeltatT" << std::endl;
                file << timeStep << std::endl;
                file << "@ Dimension" << std::endl;
                file << "3" << std::endl;
                file << "@ xsmin" << std::endl;
                file << -x << std::endl;
                file << "@ xsmax" << std::endl;
                file << x << std::endl;
                file << "@ ysmin" << std::endl;
                file << -y << std::endl;
                file << "@ ysmax" << std::endl;
                file << y << std::endl;
                file << "@ zsmin" << std::endl;
                file << 0 << std::endl;
                file << "@ zsmax" << std::endl;
                file << max_z_finder->GetValue() << std::endl;

                file.close(); // Close the file after writing
                frame++; // Increment frame counter
            }

            // Control the gate movement to initiate the flow
            if ((i > (timeOut * 2)) && status) {
                DEMSim.DoDynamicsThenSync(0); // Synchronize before changing the state
                std::cout << "gate is in motion from: " << timeStep * i << " s" << std::endl;
                DEMSim.ChangeFamily(10, 1); // Change family to fix the floors
                DEMSim.ChangeFamily(3, 4); // Allow gate to move by changing its family
                status = false; // Update status to indicate gate is in motion
            }

            // Stop the gate once it has fully opened
            if ((i >= (timeOut * (2) + gateMotion - 1)) && stopGate) {
                DEMSim.DoDynamicsThenSync(0); // Synchronize before stopping the gate
                std::cout << "gate has stopped at: " << timeStep * i << " s" << std::endl;
                DEMSim.ChangeFamily(4, 3); // Stop the gate by changing its family
                stopGate = false; // Update flag to indicate gate has stopped
            }
        }


        /* End of Part4 */

        std::cout << "The simulated time is: " << totalRunTime << " s" << std::endl;
        DEMSim.ShowTimingStats();
        DEMSim.ClearTimingStats();


        std::cout << "DEMdemo exiting..." << std::endl;
        return 0;
}




/**
 * Modifies a mesh file to adjust the dimensions of the hopper components
 *
 * This function reads a mesh file (OBJ format), modifies the vertex coordinates
 * based on the provided dimensions, and writes the updated mesh back to the file.
 * It's used to dynamically adjust the hopper geometry without creating new mesh files.
 *
 * @param filename The path to the mesh file to be modified
 * @param X The width dimension to apply
 * @param Y The depth dimension to apply
 * @param Z The height/thickness dimension to apply
 */
void modify(const std::string& filename, double X, double Y, double Z) {

    double x = X/2;
    double y = Y/2;
    double z = Z/2;
    std::string A = to_string_with_precision(x);
    std::string B = to_string_with_precision(z);
    std::string C = to_string_with_precision(y);

    // Create new vertex lines based on the provided X, Y, Z values
    std::vector<std::string> newLines = {
        "v " + A + " " + B + " -" + C,
        "v " + A + " " + B + " " + C,
        "v " + A + " -" + B + " -" + C,
        "v " + A + " -" + B + " " + C,
        "v -" + A + " " + B + " -" + C,
        "v -" + A + " -" + B + " -" + C,
        "v -" + A + " " + B + " " + C,
        "v -" + A + " -" + B + " " + C,
    };

    // Open the file for reading
    std::ifstream inFile(filename);
    if (!inFile.is_open()) {
        std::cerr << "Error opening file for reading!" << std::endl;
        return;
    }

    // Read all lines into a vector
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(inFile, line)) {
        lines.push_back(line);
    }
    inFile.close();

    // Replace the first 8 lines that start with 'v' with the new lines
    int count = 0;
    for (size_t i = 0; i < lines.size() && count < 8; ++i) {
        if (!lines[i].empty() && lines[i][0] == 'v') {
            lines[i] = newLines[count];
            ++count;
        }
    }

    // Open the file again for writing (truncating the original content)
    std::ofstream outFile(filename, std::ios::trunc);
    if (!outFile.is_open()) {
        std::cerr << "Error opening file for writing!" << std::endl;
        return;
    }

    // Write all lines back to the file
    for (const auto& l : lines) {
        outFile << l << std::endl;
    }

    outFile.close();
    std::cout << "The first 8 'v' lines have been replaced successfully." << std::endl;
}