#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "forces.h"
#include "constants.h"
#include "structs.h"
#include "nbrlist.h"

// This function calculates all forces acting on the particles (bonded and non-bonded).
// It initializes the forces array, then calculates bond-stretch, angle-bend, dihedral-torsion,
// and non-bonded forces. The total potential energy is returned.
double calculate_forces(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist, struct Vectors *p_vectors)
{
    struct Vec3D *f = p_vectors->f;
    size_t num_part = p_parameters->num_part;

    // Initialize the forces to zero for all particles
    for (size_t i = 0; i < num_part; i++)
        f[i] = (struct Vec3D){0.0, 0.0, 0.0};

    // Calculate the forces and accumulate the potential energy from each type of interaction
    double Epot = calculate_forces_bond(p_parameters, p_vectors);
    Epot += calculate_forces_angle(p_parameters, p_vectors);
    Epot += calculate_forces_dihedral(p_parameters, p_vectors);
    Epot += calculate_forces_nb(p_parameters, p_nbrlist, p_vectors);

    return Epot;
}

// This function calculates bond-stretch forces based on the current positions of the bonded particles.
// It applies the minimum image convention to calculate the distance between bonded pairs and then
// computes the force and potential energy due to the bond interaction.
double calculate_forces_bond(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
    double Epot = 0.0;
    struct Bond *bonds = p_vectors->bonds;
    size_t num_bonds = p_vectors->num_bonds;
    struct Vec3D *f = p_vectors->f;
    struct Vec3D *r = p_vectors->r;
    struct Vec3D L = p_parameters->L;
    struct Vec3D rij;
    struct Vec3D fi = {0.0, 0.0, 0.0};

    // Loop through each bond and calculate the forces
    for (size_t q = 0; q < num_bonds; ++q)
    {
        size_t i = bonds[q].i;
        size_t j = bonds[q].j;

    // Apply the minimum image convention for calculating distances
        rij.x = r[i].x - r[j].x;
        rij.x = rij.x - L.x * floor(rij.x / L.x + 0.5);
        rij.y = r[i].y - r[j].y;
        rij.y = rij.y - L.y * floor(rij.y / L.y + 0.5);
        rij.z = r[i].z - r[j].z;
        rij.z = rij.z - L.z * floor(rij.z / L.z + 0.5);

    /// \todo Provide the bond force calculation and assign forces to particles i and j
        double rij_sq = sqrt(rij.x*rij.x + rij.y*rij.y + rij.z*rij.z);

        fi.x = -k_b * (rij_sq - r0) * (rij.x/rij_sq);
        fi.y = -k_b * (rij_sq - r0) * (rij.y/rij_sq);
        fi.z = -k_b * (rij_sq - r0) * (rij.z/rij_sq);
        
        f[i].x += fi.x;
        f[i].y += fi.y;
        f[i].z += fi.z;
        f[j].x -= fi.x;
        f[j].y -= fi.y;
        f[j].z -= fi.z;

        Epot += (k_B/2) * (rij_sq - r0)*(rij_sq - r0);

    }
}


// This function calculates angle-bend forces based on the current positions of the angle-defined particles.
// It uses the minimum image convention and computes forces due to angle interactions.
double calculate_forces_angle(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
    double Epot = 0.0;
    struct Angle *angles = p_vectors->angles;
    size_t num_angles = p_vectors->num_angles;
    struct Vec3D *f = p_vectors->f;
    struct Vec3D *r = p_vectors->r;
    struct Vec3D L = p_parameters->L;
    struct Vec3D rij, rkj;
    struct Vec3D fi = {0.0, 0.0, 0.0}, fk = {0.0, 0.0, 0.0};

    // Loop through each angle and calculate the forces
    // Note: The angles triplets ijk are computed from the bond information during initialization. This is already implemented.
    for (size_t q = 0; q < num_angles; ++q)
    {
        size_t i = angles[q].i;
        size_t j = angles[q].j;
        size_t k = angles[q].k;

    // Apply the minimum image convention for calculating distances
        rij.x = r[i].x - r[j].x;
        rij.x = rij.x - L.x * floor(rij.x / L.x + 0.5);
        rij.y = r[i].y - r[j].y;
        rij.y = rij.y - L.y * floor(rij.y / L.y + 0.5);
        rij.z = r[i].z - r[j].z;
        rij.z = rij.z - L.z * floor(rij.z / L.z + 0.5);

        rkj.x = r[k].x - r[j].x;
        rkj.x = rkj.x - L.x * floor(rkj.x / L.x + 0.5);
        rkj.y = r[k].y - r[j].y;
        rkj.y = rkj.y - L.y * floor(rkj.y / L.y + 0.5);
        rkj.z = r[k].z - r[j].z;
        rkj.z = rkj.z - L.z * floor(rkj.z / L.z + 0.5);

    /// \todo Provide the angle force calculation and assign forces to particles i, j, and k
        double r1_sq = (rij.x*rij.x + rij.y*rij.y + rij.z*rij.z);
        double r2_sq = (rkj.x*rkj.x + rkj.y*rkj.y + rkj.z*rkj.z);
        double theta = acos((rij.x*rkj.x + rij.y*rkj.y + rij.z*rkj.z) / sqrt(r1_sq*r2_sq));
        double prefac = k_theta * (theta - theta0) / (sin(theta));

        fi.x = prefac * ((rkj.x/(r1_sq*r2_sq)) - (cos(theta)*rij.x)/(r1_sq*r1_sq));
        fi.y = prefac * ((rkj.y/(r1_sq*r2_sq)) - (cos(theta)*rij.y)/(r1_sq*r1_sq));
        fi.z = prefac * ((rkj.z/(r1_sq*r2_sq)) - (cos(theta)*rij.z)/(r1_sq*r1_sq));
 
        fk.x = prefac * ((rij.x/(r1_sq*r2_sq)) - (cos(theta)*rkj.x)/(r2_sq*r2_sq));
        fk.y = prefac * ((rij.y/(r1_sq*r2_sq)) - (cos(theta)*rkj.y)/(r2_sq*r2_sq));
        fk.z = prefac * ((rij.z/(r1_sq*r2_sq)) - (cos(theta)*rkj.z)/(r2_sq*r2_sq));

        f[i].x += fi.x;
        f[i].y += fi.y;
        f[i].z += fi.z;
        f[j].x -= (fi.x + fk.x);
        f[j].y -= (fi.y + fk.y);
        f[j].z -= (fi.z + fk.z);
        f[k].x += fk.x;
        f[k].y += fk.y;
        f[k].z += fk.z;

        Epot += (k_B/2) * (theta - theta0)*(theta - theta0);
    } 
    return Epot;
}

// This function calculates dihedral-torsion forces based on the positions of four connected particles.
// Ryckaert–Bellemans: U(φ) = c0 + c1 cosφ + c2 cos^2φ + c3 cos^3φ
// Forces: F_a = - (dU/dφ) * (∂φ/∂r_a)
double calculate_forces_dihedral(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
    double Epot = 0.0;
    struct Dihedral *dihedrals = p_vectors->dihedrals;
    size_t num_dihedrals = p_vectors->num_dihedrals;
    struct Vec3D *f = p_vectors->f;
    struct Vec3D *r = p_vectors->r;
    struct Vec3D L = p_parameters->L;

    // --- RB coefficients (in units of kB*T with kB=1 here). Move to Parameters if pref.
    const double c0 = 1010.0;
    const double c1 = -2018.9;
    const double c2 = 136.4;
    const double c3 = 3165.3;

    for (size_t q = 0; q < num_dihedrals; ++q)
    {
        size_t i = dihedrals[q].i;
        size_t j = dihedrals[q].j;
        size_t k = dihedrals[q].k;
        size_t l = dihedrals[q].l;

        // --- Minimum image vectors
        struct Vec3D rij, rjk, rkl;

        // i - j
        rij.x = r[i].x - r[j].x; rij.x -= L.x * floor(rij.x / L.x + 0.5);
        rij.y = r[i].y - r[j].y; rij.y -= L.y * floor(rij.y / L.y + 0.5);
        rij.z = r[i].z - r[j].z; rij.z -= L.z * floor(rij.z / L.z + 0.5);

        // k - j
        rjk.x = r[k].x - r[j].x; rjk.x -= L.x * floor(rjk.x / L.x + 0.5);
        rjk.y = r[k].y - r[j].y; rjk.y -= L.y * floor(rjk.y / L.y + 0.5);
        rjk.z = r[k].z - r[j].z; rjk.z -= L.z * floor(rjk.z / L.z + 0.5);

        // l - k
        rkl.x = r[l].x - r[k].x; rkl.x -= L.x * floor(rkl.x / L.x + 0.5);
        rkl.y = r[l].y - r[k].y; rkl.y -= L.y * floor(rkl.y / L.y + 0.5);
        rkl.z = r[l].z - r[k].z; rkl.z -= L.z * floor(rkl.z / L.z + 0.5);

        // Build bond vectors b1 = r_j - r_i = -rij ; b2 = r_k - r_j = rjk ; b3 = r_l - r_k = rkl
        struct Vec3D b1 = { -rij.x, -rij.y, -rij.z };
        struct Vec3D b2 = rjk;
        struct Vec3D b3 = rkl;

        // Norms and dots
        double b2n2 = b2.x*b2.x + b2.y*b2.y + b2.z*b2.z;
        double b2n  = sqrt(b2n2);

        // n1 = b1 x b2 ; n2 = b3 x b2
        struct Vec3D n1 = { b1.y*b2.z - b1.z*b2.y, b1.z*b2.x - b1.x*b2.z, b1.x*b2.y - b1.y*b2.x };
        struct Vec3D n2 = { b3.y*b2.z - b3.z*b2.y, b3.z*b2.x - b3.x*b2.z, b3.x*b2.y - b3.y*b2.x };

        double n1n2 = n1.x*n2.x + n1.y*n2.y + n1.z*n2.z;
        double n1n2_abs, n1n, n2n;
        double n1n2_den;

        double n1n2a = n1.x*n1.x + n1.y*n1.y + n1.z*n1.z;
        double n2n2a = n2.x*n2.x + n2.y*n2.y + n2.z*n2.z;

        // Guard against collinear cases
        if (b2n < 1e-12 || n1n2a < 1e-24 || n2n2a < 1e-24)
            continue;

        n1n = sqrt(n1n2a);
        n2n = sqrt(n2n2a);
        n1n2_den = n1n * n2n;

        // cosφ and sinφ (atan2-free; we only need cos and sin)
        double cosphi = n1n2 / n1n2_den;
        if (cosphi >  1.0) cosphi =  1.0;
        if (cosphi < -1.0) cosphi = -1.0;

        // sinφ sign from b1·n2
        double b1_dot_n2 = b1.x*n2.x + b1.y*n2.y + b1.z*n2.z;
        double sinphi = (b2n * b1_dot_n2) / n1n2_den;

        // Energy (RB up to cos^3)
        double cos2 = cosphi * cosphi;
        double cos3 = cos2 * cosphi;
        double U = c0 + c1*cosphi + c2*cos2 + c3*cos3;

        // dU/dφ = -sinφ * (c1 + 2 c2 cosφ + 3 c3 cos^2φ)
        double dUdphi = -sinphi * (c1 + 2.0*c2*cosphi + 3.0*c3*cos2);

        // ∂φ/∂r vectors (standard dihedral gradients)
        // dphi/di = - (|b2| / |n1|^2) * n1
        // dphi/dl =   (|b2| / |n2|^2) * n2
        // dphi/dj =  ( (b1·b2)/|b2|^2 ) dphi/di - ( (b3·b2)/|b2|^2 ) dphi/dl
        // dphi/dk = -dphi/di - dphi/dj - dphi/dl
        double inv_n1n2a = 1.0 / n1n2a;
        double inv_n2n2a = 1.0 / n2n2a;
        double inv_b2n2  = 1.0 / b2n2;

        struct Vec3D dphi_di = { -b2n * n1.x * inv_n1n2a, -b2n * n1.y * inv_n1n2a, -b2n * n1.z * inv_n1n2a };
        struct Vec3D dphi_dl = {  b2n * n2.x * inv_n2n2a,  b2n * n2.y * inv_n2n2a,  b2n * n2.z * inv_n2n2a };

        double b1_dot_b2 = b1.x*b2.x + b1.y*b2.y + b1.z*b2.z;
        double b3_dot_b2 = b3.x*b2.x + b3.y*b2.y + b3.z*b2.z;

        struct Vec3D dphi_dj = {
            ( b1_dot_b2 * dphi_di.x - b3_dot_b2 * dphi_dl.x ) * inv_b2n2,
            ( b1_dot_b2 * dphi_di.y - b3_dot_b2 * dphi_dl.y ) * inv_b2n2,
            ( b1_dot_b2 * dphi_di.z - b3_dot_b2 * dphi_dl.z ) * inv_b2n2
        };
        struct Vec3D dphi_dk = {
            -(dphi_di.x + dphi_dj.x + dphi_dl.x),
            -(dphi_di.y + dphi_dj.y + dphi_dl.y),
            -(dphi_di.z + dphi_dj.z + dphi_dl.z)
        };

        // Forces: F = - dU/dφ * ∂φ/∂r
        struct Vec3D fi = { -dUdphi * dphi_di.x, -dUdphi * dphi_di.y, -dUdphi * dphi_di.z };
        struct Vec3D fj = { -dUdphi * dphi_dj.x, -dUdphi * dphi_dj.y, -dUdphi * dphi_dj.z };
        struct Vec3D fk = { -dUdphi * dphi_dk.x, -dUdphi * dphi_dk.y, -dUdphi * dphi_dk.z };
        struct Vec3D fl = { -dUdphi * dphi_dl.x, -dUdphi * dphi_dl.y, -dUdphi * dphi_dl.z };

        // Accumulate
        f[i].x += fi.x; 
        f[i].y += fi.y; 
        f[i].z += fi.z;
        f[j].x += fj.x;
        f[j].y += fj.y; 
        f[j].z += fj.z;
        f[k].x += fk.x; 
        f[k].y += fk.y; 
        f[k].z += fk.z;
        f[l].x += fl.x; 
        f[l].y += fl.y; 
        f[l].z += fl.z;

        Epot += U;
    }
    return Epot;
}


// This function calculates non-bonded forces between particles using the neighbor list.
// The potential energy and forces are calculated using the Lennard-Jones potential.
/*
 * This function computes the non-bonded interactions between particles
 * using the Lennard–Jones potential with an energy shift at the cutoff.
 * For each pair of neighbors from the neighbor list, the Lorentz–Berthelot
 * mixing rules are applied to obtain σ_ij and ε_ij. If the distance is 
 * smaller than the cutoff, the truncated and shifted Lennard–Jones potential 
 * is evaluated to accumulate the total potential energy, and the corresponding 
 * force vector is calculated as the negative gradient of the potential.
 * The force is then distributed to the two particles (equal and opposite)
 * according to Newton’s third law. In this way, the routine updates both
 * the forces and the non-bonded contribution to the energy of the system.
 */

double calculate_forces_nb(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist, struct Vectors *p_vectors)
{
    struct Vec3D df;
    double r_cutsq, sigmasq, sr2, sr6, sr12, fr;
    struct DeltaR rij;
    struct Pair *nbr = p_nbrlist->nbr;
    const size_t num_nbrs = p_nbrlist->num_nbrs;
    struct Vec3D *f = p_vectors->f;
    size_t num_part = p_parameters->num_part;

    r_cutsq = p_parameters->r_cut * p_parameters->r_cut;
    double Epot = 0.0, Epot_cutoff;
    // The energy shift (Epot_cutoff) will be calculated per pair using sigma_ij (Lorentz-Berthelot mixing rule).
    // This is correct for mixtures, since each pair can have a different sigma.

    // Loop through the neighbor list and calculate the forces for each particle pair
    for (size_t k = 0; k < num_nbrs; k++)
    {
        rij = nbr[k].rij;
        size_t i = nbr[k].i;
        size_t j = nbr[k].j;

    // Compute forces if the distance is smaller than the cutoff distance
        if (rij.sq < r_cutsq)
        {
            // Lorentz-Berthelot mixing rules for sigma and epsilon
            int type_i = p_vectors->type[i];
            int type_j = p_vectors->type[j];
            double sigma_ij = 0.5 * (p_parameters->sigma[type_i] + p_parameters->sigma[type_j]);
            double epsilon_ij = sqrt(p_parameters->epsilon[type_i] * p_parameters->epsilon[type_j]);
            sigmasq = sigma_ij * sigma_ij;

            sr2 = sigmasq / rij.sq;
            sr6 = sr2 * sr2 * sr2;
            sr12 = sr6 * sr6;

            // --- IMPORTANT: Per-pair energy shift ---
            // To ensure the truncated potential is continuous in mixtures,
            // the shift (Epot_cutoff) must be calculated using sigma_ij for each pair.
            double sr2c  = sigmasq / r_cutsq;
            double sr6c  = sr2c * sr2c * sr2c;
            double sr12c = sr6c * sr6c;
            double Epot_cutoff = sr12c - sr6c;

        
            // Calculate the potential energy
            Epot += 4.0 * epsilon_ij * (sr12 - sr6 - Epot_cutoff);

            // Compute the force and apply it to both particles
            fr = 24.0 * epsilon_ij * (2.0 * sr12 - sr6) / rij.sq;  // Force divided by distance
            df.x = fr * rij.x;
            df.y = fr * rij.y;
            df.z = fr * rij.z;

            f[i].x += df.x;
            f[i].y += df.y;
            f[i].z += df.z;
            f[j].x -= df.x;
            f[j].y -= df.y;
            f[j].z -= df.z;
        }
    }

    return Epot;  // Return the potential energy due to non-bonded interactions
}

/*
 * This routine computes BOTH the non-bonded potential energy and the forces:
 *
 * 1. Potential Energy:
 *    - For each neighbor pair (i,j) with distance < r_cut, the truncated and
 *      shifted Lennard–Jones potential is evaluated.
 *    - The contributions are accumulated into the variable Epot.
 *    - At the end, the function returns the total non-bonded potential energy.
 *
 * 2. Forces:
 *    - From the same potential, the force vector is obtained as the derivative
 *      with respect to the distance.
 *    - This vector is added to particle i and subtracted from particle j
 *      (Newton's third law: equal and opposite).
 *    - The result is stored in p_vectors->f, which contains the forces acting
 *      on all particles.
 *
 * In summary:
 * - The return value is the total non-bonded potential energy.
 * - The particle force vectors are updated as a side effect.
 */

 