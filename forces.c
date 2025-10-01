#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "forces.h"
#include "constants.h"
#include "structs.h"
#include "nbrlist.h"
#include "vec3d.h"

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
    // Epot += calculate_forces_angle(p_parameters, p_vectors);
    // Epot += calculate_forces_dihedral(p_parameters, p_vectors);
    // Epot += calculate_forces_nb(p_parameters, p_nbrlist, p_vectors);

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

    /// \todo (done) Provide the bond force calculation and assign forces to particles i and j
        double rij_sq = sqrt(rij.x*rij.x + rij.y*rij.y + rij.z*rij.z);
        if (rij_sq < 1e-8) rij_sq = 1e-8;

        fi.x = -k_b * (rij_sq - r0) * (rij.x/rij_sq);
        fi.y = -k_b * (rij_sq - r0) * (rij.y/rij_sq);
        fi.z = -k_b * (rij_sq - r0) * (rij.z/rij_sq);
        
        f[i].x += fi.x;
        f[i].y += fi.y;
        f[i].z += fi.z;
        f[j].x -= fi.x;
        f[j].y -= fi.y;
        f[j].z -= fi.z;

        Epot += (k_b/2) * (rij_sq - r0)*(rij_sq - r0);

    }
    return Epot;
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

    /// \todo (done) Provide the angle force calculation and assign forces to particles i, j, and k
        double rij_sq = (rij.x*rij.x + rij.y*rij.y + rij.z*rij.z);
        double rkj_sq = (rkj.x*rkj.x + rkj.y*rkj.y + rkj.z*rkj.z);

        if (rij_sq < 1e-8) rij_sq = 1e-8;
        if (rkj_sq < 1e-8) rkj_sq = 1e-8;

        double cos_theta = (rij.x*rkj.x + rij.y*rkj.y + rij.z*rkj.z) / sqrt(rij_sq*rkj_sq);

        if (cos_theta > 1) cos_theta = 1;
        if (cos_theta < -1) cos_theta = -1;

        double theta = acos(cos_theta);

        double dtheta = theta - theta0;
        double sin_theta = sin(theta);

        if (fabs(sin_theta) < 1e-8) 
        {
            if (sin_theta >= 0.0) 
                sin_theta = 1e-8;  
            else 
                sin_theta = -1e-8;  
        }

        // double theta = acos((rij.x*rkj.x + rij.y*rkj.y + rij.z*rkj.z) / sqrt(rij_sq*rkj_sq));
        double prefac = k_theta * dtheta / sin_theta;

        fi.x = (prefac / rij_sq) * (rkj.x - (cos_theta*rij.x));
        fi.y = (prefac / rij_sq) * (rkj.y - (cos_theta*rij.y));
        fi.z = (prefac / rij_sq) * (rkj.z - (cos_theta*rij.z));

        fk.x = (prefac / rkj_sq) * (rij.x - (cos_theta*rkj.x));
        fk.y = (prefac / rkj_sq) * (rij.y - (cos_theta*rkj.y));
        fk.z = (prefac / rkj_sq) * (rij.z - (cos_theta*rkj.z));

        f[i].x += fi.x;
        f[i].y += fi.y;
        f[i].z += fi.z;
        f[j].x -= (fi.x + fk.x);
        f[j].y -= (fi.y + fk.y);
        f[j].z -= (fi.z + fk.z);
        f[k].x += fk.x;
        f[k].y += fk.y;
        f[k].z += fk.z;
        
        Epot += (k_theta/2) * (dtheta*dtheta);
    } 
    return Epot;
}

// This function calculates dihedral-torsion forces based on the positions of four connected particles.
double calculate_forces_dihedral(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
    double Epot = 0.0;
    struct Dihedral *dihedrals = p_vectors->dihedrals;
    size_t num_dihedrals = p_vectors->num_dihedrals;
    struct Vec3D *f = p_vectors->f;
    struct Vec3D *r = p_vectors->r;
    struct Vec3D L = p_parameters->L;
    struct Vec3D rij, rkj, rkl;
    struct Vec3D fi, fj, fk, fl = {0.0, 0.0, 0.0};

    // --- RB coefficients (in units of kB*T). Move to Parameters if pref.
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

        rij.x = r[i].x - r[j].x; 
        rij.x -= L.x * floor(rij.x / L.x + 0.5);
        rij.y = r[i].y - r[j].y; 
        rij.y -= L.y * floor(rij.y / L.y + 0.5);
        rij.z = r[i].z - r[j].z; 
        rij.z -= L.z * floor(rij.z / L.z + 0.5);

        rkj.x = r[k].x - r[j].x;
        rkj.x = rkj.x - L.x * floor(rkj.x / L.x + 0.5);
        rkj.y = r[k].y - r[j].y;
        rkj.y = rkj.y - L.y * floor(rkj.y / L.y + 0.5);
        rkj.z = r[k].z - r[j].z;
        rkj.z = rkj.z - L.z * floor(rkj.z / L.z + 0.5);

        rkl.x = r[l].x - r[k].x; 
        rkl.x -= L.x * floor(rkl.x / L.x + 0.5);
        rkl.y = r[l].y - r[k].y; 
        rkl.y -= L.y * floor(rkl.y / L.y + 0.5);
        rkl.z = r[l].z - r[k].z; 
        rkl.z -= L.z * floor(rkl.z / L.z + 0.5);

        struct Vec3D nj = cross(rij, rkj);
        struct Vec3D nk = cross(rkj, rkl);
        double nj_mag = norm(nj);
        double nk_mag = norm(nk);

        if (nj_mag < 1e-2) nj_mag = 1e-2;
        if (nk_mag < 1e-2) nk_mag = 1e-2;

        struct Vec3D nj_hat = scl(1.0/ nj_mag, nj);
        struct Vec3D nk_hat = scl(1.0/ nk_mag, nk);

        double cos_phi = dot(nj_hat, nk_hat);
        if (cos_phi >  1.0) cos_phi =  1.0;
        if (cos_phi < -1.0) cos_phi = -1.0;

        double dU_dcosphi = c1 + 2 * c2 * cos_phi + 3 * c3 * cos_phi * cos_phi;

        // Gradients calculation
        double rkj_mag = norm(rkj);
        double rkj_sq = rkj_mag * rkj_mag;
        struct Vec3D grad_ij = scl(-dot(rkl, nj_hat) * rkj_sq / (nj_mag * nk_mag), nj_hat);
        struct Vec3D grad_kl = scl(-dot(rij, nk_hat) * rkj_sq / (nj_mag * nk_mag), nk_hat);

        struct Vec3D grad_kj_1 = scl(dot(rkj, rij) * dot(rkl, nj_hat) / (nj_mag * nk_mag), nj_hat);
        struct Vec3D grad_kj_2 = scl(dot(rkj, rkl) * dot(rij, nk_hat) / (nj_mag * nk_mag), nk_hat);
        struct Vec3D grad_kj = add(grad_kj_1, grad_kj_2);

        // Atom gradient calculations
        struct Vec3D grad_i = grad_ij;
        struct Vec3D grad_j = add(scl(-1.0, grad_ij), scl(-1.0, grad_kj));
        struct Vec3D grad_k = add(grad_kj, grad_kl);
        struct Vec3D grad_l = scl(-1.0, grad_kl);

        // Forces calculation
        struct Vec3D fi = scl(-dU_dcosphi, grad_i);
        struct Vec3D fj = scl(-dU_dcosphi, grad_j);
        struct Vec3D fk = scl(-dU_dcosphi, grad_k);
        struct Vec3D fl = scl(-dU_dcosphi, grad_l);

        //Add forces
        f[i] = add(f[i], fi);
        f[j] = add(f[j], fj);
        f[k] = add(f[k], fk);
        f[l] = add(f[l], fl);

        Epot += c0 + c1*cos_phi + c2*cos_phi*cos_phi + c3*cos_phi*cos_phi*cos_phi;
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
    double sr2_cut, sr6_cut, sr12_cut;
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
            // Skip the LJ for bonded interactions
            if ( (p_parameters->exclude_12_nb && is_connected_12(i,j,p_nbrlist)) || (p_parameters->exclude_13_nb && is_connected_13(i,j,p_nbrlist)) || (p_parameters->exclude_14_nb && is_connected_14(i,j,p_nbrlist)) )
                continue;

            // Lorentz-Berthelot mixing rules for sigma and epsilon
            int type_i = p_vectors->type[i];
            int type_j = p_vectors->type[j];

            double sigma_ij = 0.5 * (p_parameters->sigma[type_i] + p_parameters->sigma[type_j]);
            double sigma_ij_sq = sigma_ij * sigma_ij;
            double epsilon_ij = sqrt(p_parameters->epsilon[type_i] * p_parameters->epsilon[type_j]);
            sigmasq = sigma_ij * sigma_ij;

            //Calculate Epot_cut
            sr2_cut = sigmasq / r_cutsq;
            sr6_cut = sr2_cut * sr2_cut * sr2_cut;
            sr12_cut = sr6_cut * sr6_cut;
            Epot_cutoff = sr12_cut - sr6_cut;

            // Calculate LJ potential
            sr2 = sigmasq / rij.sq;
            sr6 = sr2 * sr2 * sr2;
            sr12 = sr6 * sr6;

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
 * - The return cos_theta is the total non-bonded potential energy.
 * - The particle force vectors are updated as a side effect.
 */

 