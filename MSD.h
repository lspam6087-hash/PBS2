/* MSD.h */
#ifndef MSD_H
#define MSD_H

#include <stddef.h>
#include "structs.h"

/**
 * @brief Initializes the MSD parameters.
 * 
 * @param p_parameters used members: num_part, num_dt_steps
 * @param p_vectors used members: r
 */
void initialize_msd(struct Parameters *p_parameters, struct Vectors *p_vectors);

#endif