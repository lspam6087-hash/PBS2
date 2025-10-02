#ifndef FILEOUTPUT_H_
#define FILEOUTPUT_H_

/**
 * @brief Output particle positions to a pdb file for visualization in tools such as OVITO.
 * @param reset 1: write new file or overwrite existing file, 0: append data
 * @param p_parameters used members: filename_pdb
 * @param p_vectors used members: r
 * @param time time stamp
 */
void record_trajectories_pdb(int reset, struct Parameters * p_parameters, struct Vectors * p_vectors, double time);

/**
 * @brief Output particle positions to xyz file
 * @param reset 1: write new file or overwrite existing file, 0: append data
 * @param p_parameters used members: filename_xyz
 * @param p_vectors used members: r
 * @param time time stamp
 */
void record_trajectories_xyz(int reset, struct Parameters * p_parameters, struct Vectors * p_vectors, double time);

/**
 * @brief Save a restart file
 * 
 * @param p_parameters 
 * @param p_vectors 
 */
void save_restart(struct Parameters *p_parameters, struct Vectors *p_vectors);

/**
 * @brief Load a restart file
 * 
 * @param p_parameters 
 * @param p_vectors 
 */
void load_restart(struct Parameters *p_parameters, struct Vectors *p_vectors, struct MSD *p_msd);

/** 
 * @brief Output diagnostic data (energies, temperature) to a CSV file.
 * 
 * Writes one line per timestep containing time, kinetic energy, potential energy,
 * total energy, and instantaneous temperature.
 * 
 * @param reset 1: write new file (and CSV header), 0: append data
 * @param p_parameters used members: filename_diag
 * @param time simulation time stamp
 * @param kin_energy kinetic energy at this step
 * @param pot_energy potential energy at this step
 * @param temperature instantaneous temperature at this step
 */
void record_diagnostics_csv(int reset, struct Parameters *p_parameters, double time, double kin_energy, double pot_energy, double temperature);

/** 
 * @brief Output MSD data to a CSV file.
 * 
 * @param p_parameters used members: MSD_filename, num_dtau, dt
 * @param P_msd used members: count, cor
 */
void record_msd_csv(struct Parameters *p_parameters, struct MSD *p_msd);

#endif /* FILEOUTPUT_H_ */