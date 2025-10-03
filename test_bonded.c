#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "constants.h"
#include "structs.h"
#include "forces.h"

static void setup_single_molecule(struct Parameters *p, struct Vectors *v)
{

    p->num_part = 4;
    p->L = (struct Vec3D){50.0, 50.0, 50.0};

    // Avoid non bonded forces
    for (int t = 0; t < 2; ++t) p->epsilon[t] = 0.0;
    p->r_cut = 1e-6;
    p->exclude_12_nb = p->exclude_13_nb = p->exclude_14_nb = 1;

    // Types
    v->type[0]=0; v->type[1]=1; v->type[2]=1; v->type[3]=0;

    
    v->r[0] = (struct Vec3D){0.0,     0.000,  0.000};
    v->r[1] = (struct Vec3D){r0,      0.003,  0.000};
    v->r[2] = (struct Vec3D){2.0*r0, -0.004,  0.001};
    v->r[3] = (struct Vec3D){3.0*r0,  0.000,  0.000};


    v->num_bonds = 3;
    v->bonds = (struct Bond*)malloc(3*sizeof(struct Bond));
    v->bonds[0] = (struct Bond){0,1};
    v->bonds[1] = (struct Bond){1,2};
    v->bonds[2] = (struct Bond){2,3};

    v->num_angles = 2;
    v->angles = (struct Angle*)malloc(2*sizeof(struct Angle));
    v->angles[0] = (struct Angle){0,1,2};
    v->angles[1] = (struct Angle){1,2,3};

    v->num_dihedrals = 1;
    v->dihedrals = (struct Dihedral*)malloc(sizeof(struct Dihedral));
    v->dihedrals[0] = (struct Dihedral){0,1,2,3};
}

static void free_single_molecule(struct Vectors *v)
{
    free(v->bonds);     v->bonds = NULL;     v->num_bonds = 0;
    free(v->angles);    v->angles = NULL;    v->num_angles = 0;
    free(v->dihedrals); v->dihedrals = NULL; v->num_dihedrals = 0;
}

static void force_only_bonded(struct Parameters *p, struct Vectors *v,
                              int atom_id, double *U, struct Vec3D *F)
{
    
    for (size_t i = 0; i < p->num_part; ++i) v->f[i] = (struct Vec3D){0,0,0};
    *U = calculate_forces(p, NULL, v);
    *F = v->f[atom_id];
}

void run_bonded_fd_test(struct Parameters *p_global, struct Vectors *v)
{
    struct Parameters p = *p_global;  
    setup_single_molecule(&p, v);

    const double h = 1e-4 * r0;

    for (int a = 0; a < 4; ++a) {
        struct Vec3D r0v = v->r[a];
        double U0, Up, Um; struct Vec3D Fa, dum;

        // X
        force_only_bonded(&p, v, a, &U0, &Fa);
        v->r[a].x = r0v.x + h; force_only_bonded(&p, v, a, &Up, &dum);
        v->r[a].x = r0v.x - h; force_only_bonded(&p, v, a, &Um, &dum);
        double Fx_num = -(Up - Um) / (2*h);
        double errx = fabs(Fa.x - Fx_num) / fmax(1.0, fabs(Fa.x));

        // Y
        v->r[a] = r0v;         force_only_bonded(&p, v, a, &U0, &Fa);
        v->r[a].y = r0v.y + h; force_only_bonded(&p, v, a, &Up, &dum);
        v->r[a].y = r0v.y - h; force_only_bonded(&p, v, a, &Um, &dum);
        double Fy_num = -(Up - Um) / (2*h);
        double erry = fabs(Fa.y - Fy_num) / fmax(1.0, fabs(Fa.y));

        // Z
        v->r[a] = r0v;         force_only_bonded(&p, v, a, &U0, &Fa);
        v->r[a].z = r0v.z + h; force_only_bonded(&p, v, a, &Up, &dum);
        v->r[a].z = r0v.z - h; force_only_bonded(&p, v, a, &Um, &dum);
        double Fz_num = -(Up - Um) / (2*h);
        double errz = fabs(Fa.z - Fz_num) / fmax(1.0, fabs(Fa.z));

        printf("Atom %d | Fx ana=% .6e num=% .6e err=%.3e | "
               "Fy ana=% .6e num=% .6e err=%.3e | "
               "Fz ana=% .6e num=% .6e err=%.3e\n",
               a, Fa.x, Fx_num, errx, Fa.y, Fy_num, erry, Fa.z, Fz_num, errz);

               
    }


    free_single_molecule(v);
}
