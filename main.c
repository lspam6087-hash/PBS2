/******************************************************************************/ 
/*                                                                            */
/*  A Molecular Dynamics simulation of Lennard-Jones particles                */
/*                                                                            */
/*	This code is part of the course "Particle-based Simulations"              */
/*  taught at Eindhoven University of Technology.                             */
/*  No part of this code may be reproduced without permission of the author:  */
/*  Dr. Ir. E.A.J.F. Peters                                                   */
/*                                                                            */
/*  Dr. Ir. J.T. Padding:    version 1.1, 30/1/2013                           */
/*  Jeroen Hofman:           version 1.2, 28/7/2015                           */
/*  Dr. Ir. E.A.J.F. Peters: version 6.1, 17/9/2025                           */
/******************************************************************************/ 

/** 
 * For the 2025 PBS assignment, the code needs to be extended:
 * 
 * - Initialize vectors.type so particles get the proper type 
 * - Implement bonds in initialise_bonds in file initialise.c 
 * - Implement the bonded and non-bonded force in forces.c (Make forces type-dependent) 
 * - Change the particle position initialization to account for bond lengths, angles and dihedrals 
 * - Implement a Berendsen thermostat in dynamics.c 
 * - Implement the needed on-the-fly data analysis
 */ 

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include "constants.h" 
#include "structs.h" 
#include "setparameters.h" 
#include "initialise.h" 
#include "nbrlist.h" 
#include "forces.h" 
#include "dynamics.h" 
#include "memory.h" 
#include "fileoutput.h" 
#include "tests.h"
// #define RUN_TEST_NB
#define RUN_TEST_NVE     // ← activa el test NVE


/** 
 * @brief Main MD simulation code. After initialization, 
 * a velocity-Verlet scheme is executed for a specified number of time steps. 
 * 
 * This simulation runs for a binary mixture of methane and ethane. 
 * 
 * Tasks implemented in this assignment:
 * 
 * 1. Add support for multiple particle types (CH3 and CH2) and handle non-bonded forces accordingly.
 * 2. Modify the code to include bonded interactions for n-butane molecules.
 * 3. Change the particle position initialization to account for bond lengths, angles and dihedrals.
 * 4. Implement the Berendsen thermostat to maintain the system temperature (NVT ensemble).
 * 5. Perform on-the-fly analysis for velocity distribution, torsion angle distribution and mean-square displacement.
 * 
 * @return int 0 on success, non-zero on failure. 
 */ 
int main(void) 
{ 
    struct Vectors vectors; 
    struct Parameters parameters; 
    struct Nbrlist nbrlist; 
    size_t step; 
    double Ekin, Epot, time, T_meas, T_therm; 




    // Step 1: Set the simulation parameters from input files
    set_parameters(&parameters); 

    // Step 2: Allocate memory for particles, forces, and neighbor lists
    alloc_memory(&parameters, &vectors, &nbrlist); 

    #ifdef RUN_TEST_NB
    // Run non-bonded force verification test
    run_pair_case(0,0, &parameters, &vectors, &nbrlist); // CH3–CH3
    run_pair_case(0,1, &parameters, &vectors, &nbrlist); // CH3–CH2
    run_pair_case(1,1, &parameters, &vectors, &nbrlist); // CH2–CH2
    return 0; // exit program after the test
#endif


    // Check if a restart is required
    if (parameters.load_restart == 1) 
    { 
        load_restart(&parameters, &vectors); 
        initialise_structure(&parameters, &vectors, &nbrlist); 
        step = 0; 
        time = 0.0; 
    } 
    else 
    {   
    /// \todo Initialize particle types (CH3 and CH2) in vectors.type array
    /// \todo Implement the bonds between the UA of n-butane in initialise_bonds (initialise.c)

        initialise(&parameters, &vectors, &nbrlist, &step, &time); 
    }

    // Step 3: Build the neighbor list for non-bonded interactions
    build_nbrlist(&parameters, &vectors, &nbrlist); 

    // Step 4: Calculate initial forces (non-bonded and bonded if implemented)
    Epot = calculate_forces(&parameters, &nbrlist, &vectors); 

    // Output initial particle positions in PDB format
    record_trajectories_pdb(1, &parameters, &vectors, time); 

    // Main MD loop using velocity-Verlet integration
    while (step < parameters.num_dt_steps) 
    { 
        step++; 
        time += parameters.dt; 

        // Update velocities (half-step)
        /// \todo Implement the use of type-dependent masses
        Ekin = update_velocities_half_dt(&parameters, &nbrlist, &vectors); 

        T_meas = calc_temp(&parameters, Ekin); // Calculate instantaneous temperature LAURA B1

        /// \todo Implement and apply the Berendsen thermostat to maintain temperature (dynamics.c)
        #ifndef RUN_TEST_NVE
        T_therm = thermostat(&parameters, &vectors, Ekin, T_meas); 
        #endif

        // Update positions
        update_positions(&parameters, &nbrlist, &vectors); 

        // Apply boundary conditions
        boundary_conditions(&parameters, &vectors); 

        // Rebuild neighbor list if needed
        update_nbrlist(&parameters, &vectors, &nbrlist); 

        // Calculate forces for the current configuration (bonded forces if implemented)
        Epot = calculate_forces(&parameters, &nbrlist, &vectors); 

        // Final velocity update (half-step)
        Ekin = update_velocities_half_dt(&parameters, &nbrlist, &vectors); 


#ifdef RUN_TEST_NVE
static int have_ref = 0;
static double Etot0 = 0.0;
double Etot = Epot + Ekin;
if (!have_ref) { Etot0 = Etot; have_ref = 1; }

if (step % 1000 == 0) {
    double drift_rel = fabs(Etot - Etot0) / fmax(1.0, fabs(Etot0));
    printf("[NVE] step %lu  Etot=%.8e  drift_rel=%.3e\n",
           (unsigned long)step, Etot, drift_rel);
}
#endif
    
        // Output system state every 'num_dt_pdb' steps
        if (step % parameters.num_dt_pdb == 0) 
            record_trajectories_pdb(0, &parameters, &vectors, time); 

        // Save restart file every 'num_dt_restart' steps
        if (step % parameters.num_dt_restart == 0) 
            save_restart(&parameters, &vectors); 

        /// \todo Implement on-the-fly analysis of velocity distribution, torsion angle distribution and mean-square displacement

        // Print to the screen to monitor the progress of the simulation
        /// \todo Write the output (also) to file, and extend the output with temperature
        printf("Step %lu, Time %f, Epot %f, Ekin %f, Etot %f, Tmeas %f\n", (long unsigned)step, time, Epot, Ekin, Epot + Ekin, T_meas);

        // NEW: write diagnostics to CSV file
        record_diagnostics_csv((step == 1) ? 1 : 0, &parameters, time,
                       Ekin, Epot, T_meas); 
    } 



    // Save final state
    save_restart(&parameters, &vectors); 

    // Step 5: Free memory and clean up
    free_memory(&vectors, &nbrlist); 

    return 0; 
}

/*
 * General workflow of the Molecular Dynamics simulation:
 *
 * 1. Initialization
 *    - Read simulation parameters (dt, number of steps, masses, σ, ε, cutoff, etc.).
 *    - Allocate memory for positions, velocities, forces, and types.
 *    - Initialize particle positions and velocities (or load from restart).
 *    - Build the initial neighbor list.
 *    - Compute the initial forces and potential energy.
 *
 * 2. Main simulation loop (repeated every time step)
 *    a) Update velocities by half a time step using the current forces.
 *    b) Apply the thermostat (Berendsen) to control the temperature.
 *    c) Update particle positions by a full time step.
 *    d) Apply boundary conditions (periodic box).
 *    e) Update the neighbor list if required.
 *    f) Recompute forces and potential energy at the new positions.
 *    g) Complete the second half-step velocity update.
 *    h) Calculate instantaneous temperature from the kinetic energy.
 *    i) Output system state (trajectories, diagnostics, restart files).
 *    j) Print progress information to the screen.
 *
 * 3. Finalization
 *    - Save the final state to a restart file.
 *    - Free allocated memory.
 *
 * In summary: initialization → (velocity-Verlet integration + thermostat + I/O) → finalization.
 */

