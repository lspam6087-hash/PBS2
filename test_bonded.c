#include <math.h>
#include <stdio.h>
#include "constants.h"
#include "structs.h"
#include "forces.h"
#include "initialise.h"
#include "test_bonded.h"


static inline void zero_forces(struct Vectors *v, size_t n){
    for (size_t i = 0; i < n; ++i)
        v->f[i] = (struct Vec3D){0.0, 0.0, 0.0};
}

//Analytical calcultaion nbrlist=NULL to have only bonded forces
static void calc_energy_force_bonded(struct Parameters *p, struct Vectors *v,
                                     double *Epot_ana, struct Vec3D *F_ana, int atom_id)
{
    zero_forces(v, p->num_part);
    *Epot_ana = calculate_forces(p, NULL, v);   
    *F_ana = v->f[atom_id];                     
}

//Set up one molecule only
static void setup_single_molecule(struct Parameters *p, struct Vectors *v)
{
    p->num_part = 4;

    //Only bonded
    for (int t = 0; t < 2; ++t) p->epsilon[t] = 0.0;
    p->r_cut = 1e-6;
    p->exclude_12_nb = 1;
    p->exclude_13_nb = 1;
    p->exclude_14_nb = 1;

    v->type[0] = 0;
    v->type[1] = 1;
    v->type[2] = 1;
    v->type[3] = 0;

    //Initial geometry to avoid 0 forces by symetry
    v->r[0] = (struct Vec3D){0.0,     0.000,  0.000};
    v->r[1] = (struct Vec3D){r0,      0.003,  0.000};
    v->r[2] = (struct Vec3D){2.0*r0, -0.004,  0.001};
    v->r[3] = (struct Vec3D){3.0*r0,  0.000,  0.000};


    initialise_bonds(p, v);   
}

//Test

void run_bonded_fd_test(struct Parameters *p_global, struct Vectors *v)
{
    struct Parameters p = *p_global;   
    setup_single_molecule(&p, v);

    const double h = 1e-4 * r0;

    for (int atom_id = 0; atom_id < 4; ++atom_id) {
        struct Vec3D r0v = v->r[atom_id];

        // x axis
        double U0, Up, Um;
        struct Vec3D Fa, Fdummy;

        calc_energy_force_bonded(&p, v, &U0, &Fa, atom_id);
        v->r[atom_id].x = r0v.x + h;  calc_energy_force_bonded(&p, v, &Up, &Fdummy, atom_id);
        v->r[atom_id].x = r0v.x - h;  calc_energy_force_bonded(&p, v, &Um, &Fdummy, atom_id);
        double Fx_num = -(Up - Um) / (2.0 * h);
        double errx = fabs(Fa.x - Fx_num) / fmax(1.0, fabs(Fa.x));

        // y axis
        v->r[atom_id] = r0v;
        calc_energy_force_bonded(&p, v, &U0, &Fa, atom_id);
        v->r[atom_id].y = r0v.y + h;  calc_energy_force_bonded(&p, v, &Up, &Fdummy, atom_id);
        v->r[atom_id].y = r0v.y - h;  calc_energy_force_bonded(&p, v, &Um, &Fdummy, atom_id);
        double Fy_num = -(Up - Um) / (2.0 * h);
        double erry = fabs(Fa.y - Fy_num) / fmax(1.0, fabs(Fa.y));

        // z axis
        v->r[atom_id] = r0v;
        calc_energy_force_bonded(&p, v, &U0, &Fa, atom_id);
        v->r[atom_id].z = r0v.z + h;  calc_energy_force_bonded(&p, v, &Up, &Fdummy, atom_id);
        v->r[atom_id].z = r0v.z - h;  calc_energy_force_bonded(&p, v, &Um, &Fdummy, atom_id);
        double Fz_num = -(Up - Um) / (2.0 * h);
        double errz = fabs(Fa.z - Fz_num) / fmax(1.0, fabs(Fa.z));


        printf("Atom %d | Fx ana=% .6e num=% .6e err=%.3e | "
               "Fy ana=% .6e num=% .6e err=%.3e | "
               "Fz ana=% .6e num=% .6e err=%.3e\n",
               atom_id, Fa.x, Fx_num, errx, Fa.y, Fy_num, erry, Fa.z, Fz_num, errz);
    }
}