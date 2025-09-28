// tests.c — Finite-difference test for non-bonded LJ forces (embedded in project)
// ---------------------------------------------------------------------------------
// PURPOSE
//   Verify that the analytical non-bonded Lennard–Jones (LJ) force implemented in
//   calculate_forces(...) matches the numerical derivative of the (shifted) LJ energy
//   for a two-particle system.
//
// WHAT THIS FILE DOES
//   • Places two particles at a distance r on the x-axis (y=z=0).
//   • Calls your production force routine (calculate_forces) to get the analytical force.
//   • Computes a central finite-difference derivative of the energy: 
//       F_num ≈ -[U(r+h) - U(r-h)] / (2h)
//   • Compares F_anal vs F_num at several distances r and prints the relative error.
//
// IMPORTANT IMPLEMENTATION DETAILS
//   • We temporarily set num_part=2 and DO NOT change the box size L (no re-allocation).
//   • We temporarily DISABLE 1–2/1–3/1–4 exclusions, because the test does not create any bond list.
//   • The projection direction uses r_hat = (r_i - r_j)/|r_i - r_j| to match your
//     internal definition of rij in the force accumulation (df = fr * rij).
//
// USAGE
//   • Add this file to your build.
//   • In main.c, after alloc_memory(...), call run_pair_case(...) under a compile-time flag:
//       #define RUN_TEST_NB
//       run_pair_case(0,0, ...); run_pair_case(0,1, ...); run_pair_case(1,1, ...);
//       return 0;
//   • Comment the flag to run the normal simulation again.

#include <math.h>
#include <stdio.h>
#include "constants.h"
#include "structs.h"
#include "nbrlist.h"
#include "forces.h"

// --- Zero all forces (safety) -------------------------------------------------
// Utility to clear the force array before calling calculate_forces, so we avoid
// accidentally reusing forces from a previous configuration.
static inline void zero_forces(struct Vectors *v, size_t N){
    for (size_t i = 0; i < N; ++i){
        v->f[i].x = 0.0; 
        v->f[i].y = 0.0; 
        v->f[i].z = 0.0;
    }
}

// --- energy_and_proj_force ----------------------------------------------------
// Rebuilds the neighbor list for the current coordinates, calls the production
// calculate_forces(...) to obtain E_pot and force vectors, then projects the
// force on particle i=0 along r_hat = (r_i - r_j)/|r_i - r_j|.
//
// RETURNS
//   *Epot_out  : shifted, truncated LJ potential energy returned by calculate_forces.
//   *Fproj_out : analytical force on particle 0 projected along r_hat.
static void energy_and_proj_force(struct Parameters *p, struct Vectors *v, struct Nbrlist *nl,
                                  double *Epot_out, double *Fproj_out)
{
    // Build the current neighbor list (with present positions and box)
    build_nbrlist(p, v, nl);

    // Ensure forces are zeroed before accumulation
    zero_forces(v, p->num_part);

    // Call your production routine: fills v->f[*] and returns total E_pot
    double Epot = calculate_forces(p, nl, v);

    // r_hat from j -> i (consistent with rij = r_i - r_j used in your force code)
    struct Vec3D dr = { v->r[0].x - v->r[1].x, 
                        v->r[0].y - v->r[1].y, 
                        v->r[0].z - v->r[1].z };
    double r  = sqrt(dr.x*dr.x + dr.y*dr.y + dr.z*dr.z);
    double rx = dr.x / r, ry = dr.y / r, rz = dr.z / r;

    // Project the analytical force on particle 0 onto r_hat
    double Fproj = v->f[0].x * rx + v->f[0].y * ry + v->f[0].z * rz;

    *Epot_out  = Epot;
    *Fproj_out = Fproj;
}

/*
 * run_pair_case(type_i, type_j)
 * -----------------------------------------------------------------------------
 * Runs a one-dimensional finite-difference force test for a given pair of types.
 *
 * INPUT
 *   type_i, type_j : particle types (e.g., 0=CH3, 1=CH2) to be assigned to the
 *                    two test particles.
 *
 * SIDE-EFFECTS (temporary, restored at the end)
 *   - p->num_part set to 2: we use only the first two entries of the arrays.
 *   - Exclusions (1–2/1–3/1–4) disabled: no bond list exists in this synthetic test.
 *
 * ALGORITHM
 *   1) Place particles i=0 and j=1 on the x-axis separated by r.
 *   2) Compute F_anal by calling calculate_forces(...) and projecting onto r_hat.
 *   3) Compute U(r+h) and U(r-h), then F_num = -[U(r+h) - U(r-h)] / (2h).
 *   4) Print r, F_anal, F_num and the relative error.
 */
void run_pair_case(int type_i, int type_j, struct Parameters *p, struct Vectors *v, struct Nbrlist *nl)
{
    // --- Backup fields we will modify and must restore later
    size_t num_part_backup = p->num_part;
    int save_ex12 = p->exclude_12_nb;
    int save_ex13 = p->exclude_13_nb;
    int save_ex14 = p->exclude_14_nb;

    // --- Use only two particles (arrays were already allocated larger)
    p->num_part = 2;

    // --- Disable bonded exclusions for this bondless test
    //     (avoids dereferencing head12/head13/head14 which are NULL here)
    p->exclude_12_nb = 0;
    p->exclude_13_nb = 0;
    p->exclude_14_nb = 0;

    // --- Assign types and set positions along x (keep y=z=0)
    v->type[0] = type_i;
    v->type[1] = type_j;
    v->r[0].x = 0.0; v->r[0].y = 0.0; v->r[0].z = 0.0;
    v->r[1].y = 0.0; v->r[1].z = 0.0;

    // --- Lorentz–Berthelot mixing for sigma_ij (used for step size h and r values)
    double sigma_ij = 0.5*(p->sigma[type_i] + p->sigma[type_j]);

    // Finite-difference displacement (tune 1e-5..1e-3 * sigma_ij if needed for stability)
    double h = 1e-4 * sigma_ij;

    // Radii to probe (avoid exactly r_cut to not hit the derivative discontinuity)
    double rlist[] = { 0.9*sigma_ij, 1.0*sigma_ij, 1.2*sigma_ij, 2.0*sigma_ij, p->r_cut - 1e-3 };
    int nr = (int)(sizeof(rlist)/sizeof(rlist[0]));

    printf("\n== Pair test (%d-%d)  sigma_ij=%.3f  r_cut=%.3f  h=%.1e ==\n",
           type_i, type_j, sigma_ij, p->r_cut, h);
    printf("%10s %18s %18s %14s\n", "r (A)", "F_anal", "F_num", "rel_err");

    for (int k = 0; k < nr; ++k){
        // --- Set interparticle distance r along +x
        double r = rlist[k];
        v->r[1].x = r;

        // Analytical force (via your production code) and energy at r
        double U0, F_anal;
        energy_and_proj_force(p, v, nl, &U0, &F_anal);

        // Central finite-difference on energy: U(r+h) and U(r-h)
        v->r[1].x = r + h; double Uplus, dummy;
        energy_and_proj_force(p, v, nl, &Uplus, &dummy);

        v->r[1].x = r - h; double Uminus;
        energy_and_proj_force(p, v, nl, &Uminus, &dummy);

        // Numerical force estimate (negative energy gradient)
        double F_num  = -(Uplus - Uminus) / (2.0*h);

        // Relative error (use max(1, |F_anal|) to avoid exploding ratios near zero)
        double relerr = fabs(F_anal - F_num) / fmax(1.0, fabs(F_anal));

        printf("%10.4f %18.8e %18.8e %14.3e\n", r, F_anal, F_num, relerr);
    }

    // --- Restore original exclusions and number of particles
    p->exclude_12_nb = save_ex12;
    p->exclude_13_nb = save_ex13;
    p->exclude_14_nb = save_ex14;
    p->num_part      = num_part_backup;
}
