#ifndef DYNAMICS_H_
#define DYNAMICS_H_

/**
 * @brief Update particle positions by using velocities.
 * @param[in] p_parameters member: dt
 * @param[out] p_nbrlist used members: dr
 * @param[in,out] p_vectors members r, dr, v
 */
void update_positions(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist, struct Vectors *p_vectors);

/**
 * @brief Update velocities for half a time step using forces.
 * @param[in] p_parameters used members: mass, dt
 * @param[in] p_nbrlist
 * @param[in, out] p_vectors used members: v, f
 * @return double kinetic energy
 */
double update_velocities_half_dt(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist, struct Vectors *p_vectors);

/**
 * @brief Apply boundary conditions: particles folded back in periodic box.
 * 
 * Ensures that particles moving beyond the simulation box
 * boundaries are wrapped around to the opposite side.
 * This maintains the periodic nature of the simulation box.
 * 
 * @param[in] p_parameters used members: L
 * @param[in, out] p_vectors used members: r
 */
void boundary_conditions(struct Parameters *p_parameters, struct Vectors *p_vectors);

/**
 * @brief Apply thermostat by manipulating particle velocities
 * @param[in] p_parameters parameters for the thermostat
 * @param[in, out] p_vectors used members: v
 * @param[in] Ekin current kinetic energy
 * @param[in] T_meas Measured temperature
 */
double thermostat(struct Parameters *p_parameters, struct Vectors *p_vectors, double Ekin, double T_meas);

/**
 * @brief Calculate the instantaneous measured temperature from the kinetic energy.
 *
 * Converts the total kinetic energy into temperature using the equipartition theorem:
 * 
 *     T = (2 * Ekin) / (kB * dof)
 * 
 * where:
 *   - dof = 3 * (N - 1), i.e., subtracting 3 degrees of freedom for the center of mass.
 *   - kB = 1 in reduced units.
 *
 * @param[in] p_parameters pointer to Parameters struct (provides number of particles)
 * @param[in] Ekin total kinetic energy at this time step
 * @return Instantaneous temperature (double)
 */
double calc_temp(struct Parameters *p_parameters, double Ekin);

#endif /* DYNAMICS_H_ */