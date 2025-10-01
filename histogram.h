/* velocity_histogram.h */
#ifndef VELOCITY_HISTOGRAM_H
#define VELOCITY_HISTOGRAM_H

#include <stddef.h>
#include "structs.h"

VelHist *vh_create(size_t nbins, double vmin, double vmax);
void vh_destroy(VelHist *h);
void vh_reset(VelHist *h);
void vh_add_sample(VelHist *h, double speed);
void vh_write_ascii(VelHist *h, const char *filename); /* writes bin_center count */

/**
 * @brief Include the velocity of the particles in the histogram.
 * 
 * @param p_parameters used members: nbins, hist_vmax, num_part, L
 * @param p_vectors used members: v
 * @param step current step
 */
void write_hist(struct Parameters *p_parameters, struct Vectors *p_vectors, size_t step);
#endif