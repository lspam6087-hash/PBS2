#include <math.h>
#include <stdio.h>
#include "constants.h"
#include "structs.h"
#include "nbrlist.h"
#include "forces.h"

//// Initialize the forces to zero for all particles

static inline void zero_forces(struct Vectors *v, size_t n){
    for (size_t i = 0; i < n; ++i)
        v->f[i] = (struct Vec3D){0.0, 0.0, 0.0};
}

// Calculate the analytical energy and force 
static void calc_energy_force_anal(struct Parameters *p, struct Vectors *v, struct Nbrlist *nl,
                                  double *Epot_ana, double *F_ana)
{
    
    build_nbrlist(p, v, nl);
    zero_forces(v, p->num_part);
    double Epot = calculate_forces(p, nl, v);

    //Compute the force along the x-direction
    double Fprojx = -v->f[0].x;

    *Epot_ana = Epot;
    *F_ana = Fprojx;
}

//System of two particles
void run_pair_case(int type_i, int type_j, struct Parameters *p_global,
                   struct Vectors *v, struct Nbrlist *nl)
{
    // Parameters for the test (just two particles and just LN forces)
    struct Parameters ptest = *p_global;
    ptest.num_part      = 2;
    ptest.exclude_12_nb = 0;
    ptest.exclude_13_nb = 0;
    ptest.exclude_14_nb = 0;

    // Particle type and x axis
    v->type[0] = type_i;
    v->type[1] = type_j;

    v->r[0] = (struct Vec3D){0.0, 0.0, 0.0};
    v->r[1] = (struct Vec3D){0.0, 0.0, 0.0};

    double r_test = 0.5 * p_global->r_cut;   
    double h      = 1e-4 * r_test;

    // Analytical force
    v->r[1].x = r_test;           
    double U0, F_anal;
    calc_energy_force_anal(&ptest, v, nl, &U0, &F_anal);

    // Numerical force 
    v->r[1].x = r_test + h; double Uplus, dummy;
    calc_energy_force_anal(&ptest, v, nl, &Uplus, &dummy);

    v->r[1].x = r_test - h; double Uminus;
    calc_energy_force_anal(&ptest, v, nl, &Uminus, &dummy);

    double F_num  = -(Uplus - Uminus) / (2.0*h);
    double relerr = fabs(F_anal - F_num) / fmax(1.0, fabs(F_anal));

    // Results
    printf("\n== Simple pair test (%d-%d)  r=%.4f  h=%.1e ==\n",
       type_i, type_j, r_test, h);
    printf("%10s %18s %18s %14s\n", "r (A)", "F_anal", "F_num", "rel_err");
    printf("%10.4f %18.8e %18.8e %14.3e\n", r_test, F_anal, F_num, relerr);
}