#include <stdio.h>
#include <stdlib.h>
#include "structs.h"

// void initialize_msd(struct Parameters *p_parameters, struct Vectors *p_vectors){
//     size_t N = p_parameters->num_part;
//     size_t N_frames = p_parameters->num_dt_steps;
//     size_t frame = 0;
//     struct Vec3D *r = p_vectors->r;
//     double *cor = malloc((size_t)N_frames * sizeof(double));
//     int *count = malloc((size_t)N_frames * sizeof(double));
//     double *store = malloc((size_t)N * 3 * N_frames * sizeof(double));

//     // access element (i, d, f)
//     #define STORE(i,d,f) store[((i)*3 + (d)) * N_frames + (f)]

//     for (size_t i = 0; i < N; i++) {
//         for (int d = 0; d < 3; d++) {
//             for (size_t f = 0; f < N_frames; f++) {
//                 STORE(i,d,f) = 0.0;
//             }
//         }
//     }
        
//     // reset variables
//     frame = 0;
//     for (int lag = 0; lag < N_frames; lag++) {
//         cor[lag] = 0.0;
//         count[lag] = 0;
//     }

//     for (int i = 0; i < N; i++) {
//         for (int d = 0; d < 3; d++) {
//             for (int f = 0; f < N_frames; f++) {
//                 STORE[i][d][f] = 0.0;
//             }
//         }
//     }
// }

// void update_msd(struct Parameters *p_parameters, struct Vectors *p_vectors){
//     size_t N = p_parameters->num_part;
//     size_t N_frames = p_parameters->num_dt_steps;
//     size_t frame = 0;
//     struct Vec3D *r = p_vectors->r;
//     double cor[N_frames];
//     int count[N_frames];
//     double store[N][3][N_frames];

//     // store current positions into buffer
//     for (int i = 0; i < N; i++) {
//         store[i][0][frame] = r[i].x; // x
//         store[i][1][frame] = r[i].y; // y
//         store[i][2][frame] = r[i].x; // z
//     }

//     frame += 1; // increment frame counter
// }