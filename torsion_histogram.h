#ifndef TORSION_HISTOGRAM_H
#define TORSION_HISTOGRAM_H

#include <stddef.h>
#include "structs.h"   

#ifdef __cplusplus
extern "C" {
#endif


typedef struct {
    size_t nbins;
    double amin_deg;
    double amax_deg;
    double bin_width;
    unsigned long *counts;
    unsigned long total_counts;
} PhiHist;

PhiHist *ph_create(size_t nbins, double amin_deg, double amax_deg);
void     ph_destroy(PhiHist *h);
void     ph_reset(PhiHist *h);
void     ph_add_sample(PhiHist *h, double phi_deg);
void     ph_write_ascii(PhiHist *h, const char *filename);


double dihedral_phi_rad(struct Vec3D r1, struct Vec3D r2,
                        struct Vec3D r3, struct Vec3D r4,
                        struct Vec3D L);


void write_torsion_samples(struct Parameters *p_parameters,
                           struct Vectors *p_vectors,
                           size_t step);

#ifdef __cplusplus
}
#endif
#endif /* TORSION_HISTOGRAM_H */
