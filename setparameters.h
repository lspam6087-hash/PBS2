#ifndef SETPARAMETERS_H_
#define SETPARAMETERS_H_

/* function  parameters */
/**
 * @brief Set the parameters struct used in the simulation.
 * 
 * @param[out] p_parameters Parameters of the simulation.
 */
void set_parameters(struct Parameters * p_parameters);

/**
 * @brief Set the parameters struct used in the MSD simulation.
 * 
 * @param[out] p_msd Parameters of the MSD simulation.
 */
void set_parameters_MSD(struct MSD * p_MSD);

#endif /* SETPARAMETERS_H_ */
