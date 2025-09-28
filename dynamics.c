#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constants.h"
#include "structs.h"

// This function updates particle positions using their velocities.
// The positions are advanced by one full time step (dt), and displacement vectors
// (dr) for one time step are updated. The displacement since the last neighbor 
// list creation (stored in p_nbrlist->dr) is also updated for each particle.
void update_positions(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist, struct Vectors *p_vectors)
{
    struct Vec3D dr_loc;
    struct Vec3D *r = p_vectors->r;     // Particle positions
    struct Vec3D *dr = p_vectors->dr;   // Displacement in one timestep
    struct Vec3D *v = p_vectors->v;     // Particle velocities
    struct DeltaR *dr_nbrlist = p_nbrlist->dr;  // Displacement since last neighbor list creation
    size_t num_part = p_parameters->num_part;
    double dt = p_parameters->dt;

    // Loop over all particles to update their positions
    for (size_t i = 0; i < num_part; i++)
    {
        dr_loc.x = v[i].x * dt;  // Compute displacement in x-direction for one timestep
        dr_loc.y = v[i].y * dt;  // Compute displacement in y-direction for one timestep
        dr_loc.z = v[i].z * dt;  // Compute displacement in z-direction for one timestep

        dr[i] = dr_loc;          // Store the displacement for this timestep
        r[i].x += dr_loc.x;      // Update position in x-direction
        r[i].y += dr_loc.y;      // Update position in y-direction
        r[i].z += dr_loc.z;      // Update position in z-direction

        // Update the displacement since last neighbor list creation
        dr_nbrlist[i].x += dr_loc.x;
        dr_nbrlist[i].y += dr_loc.y;
        dr_nbrlist[i].z += dr_loc.z;
        dr_nbrlist[i].sq = (dr_nbrlist[i].x) * (dr_nbrlist[i].x) + 
                           (dr_nbrlist[i].y) * (dr_nbrlist[i].y) + 
                           (dr_nbrlist[i].z) * (dr_nbrlist[i].z);  // Square of total displacement since last neighbor list creation
    }
}

// This function updates particle velocities by half a time step using the current forces.
// The updated velocities are used in the velocity-Verlet integration scheme.
// The function also calculates and returns the kinetic energy of the system.
double update_velocities_half_dt(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist, struct Vectors *p_vectors)
{
    double Ekin = 0.0;  // Initialize kinetic energy
    struct Vec3D *v = p_vectors->v;  // Particle velocities
    struct Vec3D *f = p_vectors->f;  // Forces acting on particles

    // Loop over all particles y actualiza sus velocidades usando su masa
    for (size_t i = 0; i < p_parameters->num_part; i++)
    {
        double factor = 0.5 / p_parameters->mass[p_vectors->type[i]] * p_parameters->dt; // LAURA: masa por tipo
        v[i].x += factor * f[i].x;  // Update velocity in x-direction
        v[i].y += factor * f[i].y;  // Update velocity in y-direction
        v[i].z += factor * f[i].z;  // Update velocity in z-direction
        Ekin += p_parameters->mass[p_vectors->type[i]] * (v[i].x * v[i].x + v[i].y * v[i].y + v[i].z * v[i].z);  // LAURA: masa por tipo
    }

    // Final kinetic usimng specific masses
    return 0.5 * Ekin;  // (considers each particle)
}

// This function applies periodic boundary conditions to ensure particles stay inside the simulation box.
// If a particle moves beyond the box, it is wrapped around to the opposite side.
void boundary_conditions(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
    // TODO: include neighbor list and get box size from there
    struct Vec3D invL;  // Inverse of the box size
    struct Vec3D *r = p_vectors->r;  // Particle positions
    struct Vec3D L = p_parameters->L;  // Box dimensions
    size_t num_part = p_parameters->num_part;  // Number of particles

    invL.x = 1.0 / L.x;
    invL.y = 1.0 / L.y;
    invL.z = 1.0 / L.z;

    // Loop over all particles and apply periodic boundary conditions
    for (size_t i = 0; i < num_part; i++)
    {
        r[i].x -= L.x * floor(r[i].x * invL.x);  // Apply periodic boundary in x-direction
        r[i].y -= L.y * floor(r[i].y * invL.y);  // Apply periodic boundary in y-direction
        r[i].z -= L.z * floor(r[i].z * invL.z);  // Apply periodic boundary in z-direction
    }
}

double calc_temp(struct Parameters *p_parameters, double Ekin)
{
    double N = p_parameters->num_part;
    double d = 3; //Number of degrees of freedom
    double N_free = d * (N - 1);
    double T_meas = (2 * Ekin) / (N_free);
    return T_meas;
}

// This function applies a thermostat to maintain the system's temperature.
double thermostat(struct Parameters *p_parameters, struct Vectors *p_vectors, double Ekin, double T_meas)
/// \todo Change velocities by thermostatting
{
    int N = p_parameters->num_part;
    double T_therm, lambda = 0.0;
    double dt = p_parameters->dt;
    double tau = 1000; // [s]
    double T0 = 298; // [K]
    struct Vec3D *v = p_vectors->v;
    
    lambda = sqrt(1+dt/tau * (T0/T_meas - 1));

    printf("Lambda= %lf, Tmeas=%lf\n", lambda, T_therm);

    for (int i= 0; i<N; i++)
    {
        p_vectors->v[i].x = lambda*p_vectors->v[i].x;
        p_vectors->v[i].y = lambda*p_vectors->v[i].y;
        p_vectors->v[i].z = lambda*p_vectors->v[i].z;
    }
    return T_therm;
}