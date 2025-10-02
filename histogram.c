/* velocity_histogram.c */
#include "histogram.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

// Setting the parameters
VelHist *vh_create(size_t nbins, double vmin, double vmax) {
    if (nbins == 0 || vmax <= vmin) return NULL;
    VelHist *h = (VelHist*) malloc(sizeof(VelHist));
    if (!h) return NULL;
    h->nbins = nbins;
    h->vmin = vmin;
    h->vmax = vmax;
    h->bin_width = (vmax - vmin) / (double)nbins;
    h->counts = (unsigned long*) calloc(nbins, sizeof(unsigned long));
    h->total_counts = 0;
    if (!h->counts) { free(h); return NULL; }
    return h;
}

void write_hist(struct Parameters *p_parameters, struct Vectors *p_vectors, size_t step){
    size_t nbins = p_parameters->nbins;
    double vmin = 0.0;
    double vmax = p_parameters->hist_vmax; 
    int N = p_parameters->num_part;
    struct Vec3D *v = p_vectors->v;
    VelHist *vh = vh_create(nbins, vmin, vmax);

    // sampling parameters
    size_t sample_interval = p_parameters->sample_interval; 

    if (step % sample_interval != 0) return;

    if (step == sample_interval){
        FILE *fp = fopen("speeds.dat", "w");
        // loop over particles and add speeds
        for (size_t i = 0; i < N; ++i) {
            double speed = sqrt(v[i].x*v[i].x + v[i].y*v[i].y + v[i].z*v[i].z);
            fprintf(fp, "%.12f\n", speed);
        }
        fclose(fp);
    } else {
        FILE *fp = fopen("speeds.dat", "a");
        // loop over particles and add speeds 
        for (size_t i = 0; i < N; ++i) {
        double speed = sqrt(v[i].x*v[i].x + v[i].y*v[i].y + v[i].z*v[i].z);
        fprintf(fp, "%.12f\n", speed);
        }
        fclose(fp);
    }
}
