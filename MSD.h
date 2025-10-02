/* MSD.h */
#ifndef MSD_H
#define MSD_H

#include <stddef.h>
#include "structs.h"

/**
 * @brief Initializes the MSD parameters
 * 
 * @param p_parameters used members: num_part, num_dt_steps
 * @param p_vectors used members: r
 */
void initialise_msd(struct Parameters *p_parameters, struct Vectors *p_vectors, struct MSD *p_msd);

/**
 * @brief Calculates the center of mass of the molecule
 * 
 * @param p_parameters used members: L, mass
 * @param p_vectors used members: r, type
 * @param imolecule index of the butane molecule
 */
struct Vec3D com_molecule(struct Parameters *p_parameters, struct Vectors *p_vectors, int imolecule);

/**
 * @brief updates the MSD parameters
 * 
 * @param p_parameters used members: num_part, ncor, L
 * @param p_vectors used members: r
 * @param p_msd used members: frame, prev, store, cor, frame
 */
void update_msd(struct Parameters *p_parameters, struct Vectors *p_vectors, struct MSD *p_msd);


#endif