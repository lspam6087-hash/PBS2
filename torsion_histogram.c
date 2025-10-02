#include "torsion_histogram.h"
#include "constants.h"   // 'pi' en minúscula
#include <stdio.h>
#include <math.h>


static inline double wrap_min_image(double dx, double L){
    if (dx >  0.5*L) dx -= L;
    if (dx < -0.5*L) dx += L;
    return dx;
}
static inline struct Vec3D min_image_vec(struct Vec3D a, struct Vec3D b, struct Vec3D L){
    struct Vec3D d = {a.x - b.x, a.y - b.y, a.z - b.z};
    d.x = wrap_min_image(d.x, L.x);
    d.y = wrap_min_image(d.y, L.y);
    d.z = wrap_min_image(d.z, L.z);
    return d;
}
static inline double dot3(struct Vec3D a, struct Vec3D b){
    return a.x*b.x + a.y*b.y + a.z*b.z;
}
static inline struct Vec3D cross3(struct Vec3D a, struct Vec3D b){
    return (struct Vec3D){ a.y*b.z - a.z*b.y,
                           a.z*b.x - a.x*b.z,
                           a.x*b.y - a.y*b.x };
}
static inline double norm3(struct Vec3D a){ return sqrt(dot3(a,a)); }

/* ϕ en rad (-π, π] */
double dihedral_phi_rad(struct Vec3D r1, struct Vec3D r2,
                        struct Vec3D r3, struct Vec3D r4,
                        struct Vec3D L)
{
    struct Vec3D b1 = min_image_vec(r2, r1, L);
    struct Vec3D b2 = min_image_vec(r3, r2, L);
    struct Vec3D b3 = min_image_vec(r4, r3, L);

    struct Vec3D n1 = cross3(b1, b2);
    struct Vec3D n2 = cross3(b2, b3);

    double n1n = norm3(n1), n2n = norm3(n2), b2n = norm3(b2);
    if (n1n < 1e-14 || n2n < 1e-14 || b2n < 1e-14) return 0.0;

    double x = dot3(n1, n2);
    double y = dot3(cross3(n1, b2), n2) / b2n;

    return atan2(y, x);
}


void write_torsion_samples(struct Parameters *P,
                                struct Vectors *V,
                                size_t step)
{
    if (P->sample_interval_dih <= 0) return;
    if (step % (size_t)P->sample_interval_dih != 0) return;

    FILE *fp = (step == (size_t)P->sample_interval_dih)
             ? fopen(P->filename_dih_samples, "w")
             : fopen(P->filename_dih_samples, "a");
    if (!fp) return;

    size_t Nmol = P->num_part / 4;   
    struct Vec3D L = P->L;

    for (size_t m = 0; m < Nmol; ++m) {
        size_t i = 4*m + 0, j = 4*m + 1, k = 4*m + 2, l = 4*m + 3;
        double phi_rad = dihedral_phi_rad(V->r[i], V->r[j], V->r[k], V->r[l], L);
        double phi_deg = phi_rad * (180.0 / PI);
        fprintf(fp, "%.12f\n", phi_deg);
    }
    fclose(fp);
}
