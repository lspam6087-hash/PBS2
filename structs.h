#ifndef TYPES_MD_H_
#define TYPES_MD_H_
#include <stddef.h>

/* This header file contains definitions of struct types used in the molecular dynamics code */

/**
 * @file structs.h
 * @brief Defines the core data structures used in the molecular dynamics simulation.
 *
 * This header provides type definitions (structs) that organize and store both
 * simulation parameters and dynamic particle data:
 *
 * - struct Vec3D:
 *   Represents a 3D vector with x, y, z components.
 *   Used for particle positions, velocities, and forces.
 *
 * - struct Parameters:
 *   Contains global simulation parameters such as:
 *     • number of particles
 *     • number of time steps
 *     • integration time step (dt)
 *     • simulation box dimensions
 *     • file names for output (PDB, XYZ, restart, diagnostics, ...)
 *   These values define how the simulation runs and how data is saved.
 *
 * - struct Vectors:
 *   Holds arrays of particle data:
 *     • r[] : positions
 *     • v[] : velocities
 *     • f[] : forces
 *   These are updated at every time step and represent the current state
 *   of the simulated system.
 *
 * In summary, structs.h provides the central data structures that link
 * together the simulation code: Parameters define the setup, while Vectors
 * store the evolving system state.
 */


/**
 * @brief Struct to store x, y, and z component of a 3D vector.
 * 
 */
struct Vec3D
{
    double x, y, z; //!< Three three coordinates of a 3D vector
};

/**
 * @brief Struct to store all parameters. These parameters are set by the function @ref set_parameters.
 * 
 */
struct Parameters
{
    size_t num_part;         //!< Number of particles
    size_t num_dt_steps;     //!< Number of time steps
    double dt;               //!< integration time step
    struct Vec3D L;          //!< Box sizes in 3 direction
    int exclude_12_nb;       //!< If true (=1) 1-2 connected atoms are exluded from non-bonded interactions 
    int exclude_13_nb;       //!< If true (-1) 1-3 connected atoms are exluded from non-bonded interactions 
    int exclude_14_nb;       //!< If true (=1) 1-4 connected atoms are exluded from non-bonded interactions    
    double kT;               //!< Thermal energy
    double mass[2];          //!< Mass per type: [0] = CH3, [1] = CH2
    double sigma[2];         //!< σ per type (Å)
    double epsilon[2];       //!< ε per type (K, dividido ya por kB)
    double r_cut;            //!< Cut-off distance for LJ interaction
    double r_shell;          //!< Shell thickness for neighbor list
    size_t num_dt_pdb;       //!< Number of time steps between pdb saves
    double rescale_output;   //!< Rescale factor for outputting positions, typically used to rescale lenghts to a magnitude required by specific visualization software.
    char filename_pdb[1024]; //!< filename (without extension) for pdb file
    char filename_xyz[1024]; //!< filename (without extension) for pdb file
    char load_restart;       //!< if equal 1 restart file is loaded
    char filename_diag[256]; //!< base filename for diagnostics CSV output LAURA B1
    int sample_interval;     //!< sample the speeds on every interval
    int write_interval;      //!< write histogram to disk every X steps
    size_t nbins;            //!< number of bins
    double hist_vmax;        //!< maximum velocity included in the histogram
    char filename_hist[256]; //!< filename for the .csv file
    size_t nbins_dih;           //!< number of bins for dihedral histogram
    size_t sample_interval_dih;                 
    size_t write_interval_dih;                
    char filename_dih_samples[1024];
    char filename_dih_hist[256];//!< filename for the .csv file for the histogram for dihedrals
    size_t num_dt_restart;   //!< Number of time steps between saves of restart file
    char restart_in_filename[1024];  //!< filename for loaded restart file
    char restart_out_filename[1024]; //!< filename for saved restart file
    int ncor;
    char MSD_filename[1024];
    int num_dtau;
};

/**
 * @brief Struct to store a 3D vector and its square length. This is expecially useful for connecting vectors in e.g. neighbor lists.
 * 
 */
struct DeltaR
/* Structure to store a 3D vector and its square length. */
{
    double x, y, z; //!< x, y and z coordinates 
    double sq;      //!< square length 
};

/**
 * @brief Struct to store i, j, k indices of a 3D grid
 * 
 */
struct Index3D
{
    size_t i, j, k; //!< 3 indices: i,j, k
};

/**
 * @brief Struct to store indices of bonded particles i-j
 * 
 */
struct Bond
{
    size_t i,j;
};

/**
 * @brief Struct to store indices of particles in an angle i-j-k
 * 
 */
struct Angle
{
    size_t i,j,k;
};

/**
 * @brief Struct to store indices of particles in a dihedral i-j-k-l
 * 
 */
struct Dihedral
{
    size_t i,j,k,l;
};

/**
 * @brief Struct to store a pair of particles: its indices and connecting vector
 * 
 */
struct Pair
{
    size_t i, j; //!< indices of the two particles forming a pair
    struct DeltaR rij; //!< The connecting vector between the pairs rij = r[i]-r[j] corrected for periodicity
};

/**
 * @brief Struct with pointers to all particle arrays relevant for a MD simulation
 * 
 */
struct Vectors
{
    size_t size;                //!< size of particle arrays (can be > num_part)
    size_t num_bonds;           //!< number of bonds
    size_t num_angles;          //!< number of angles 
    size_t num_dihedrals;       //!< number of dihedrals
    struct Vec3D *r;            //!< positions
    struct Vec3D *dr;           //!< displacements
    struct Vec3D *v;            //!< velocities
    struct Vec3D *f;            //!< forces
    struct Bond *bonds;         //!< bonds
    struct Angle *angles;       //!< angles
    struct Dihedral *dihedrals; //!< dihedrals
    int *type;         // type of particle [N] (TYPE_CH3 o TYPE_CH2)
};

/**
 * @brief Struct used to store a cell-linked-list
 * 
 */
struct Celllist
{
    size_t *head; //!< head[icell] provides the head the list for cell icell 
    size_t *list; //!< list[i] provides the next particle index in the cell-linked-list. list[i]==SIZE_MAX encodes the end of the list.
    size_t *particle2cell; //!< provides the cell index for a particle
    size_t num_cells, num_cells_max, num_part_max; //!< number of cells used and number of cells and particles allocated for
    struct Index3D size_grid; //!< number of cells in each direction
};

/**
 * @brief Struct to store a neighbor list
 * 
 */
struct Nbrlist
{
    struct Celllist *p_celllist;   //!< pointer to celllist used to create the neighbor list
    size_t num_nbrs, num_nbrs_max; //!< number of neighbors and maximum number allocated
    struct Pair *nbr;              //!< list of non-bounded neighbor pairs
    struct DeltaR *dr;             //!< displacements particles with respect to nbrlist creation time
    size_t *head12, *pairs12;          //!< list of 12 bonded pairs
    size_t *head13, *pairs13;          //!< list of 13 bonded pairs
    size_t *head14, *pairs14;          //!< list of 14 bonded pairs
};

/**
 * @brief Struct to store data for a histogram
 * 
 */
typedef struct {
    size_t nbins;               //!< Number of bins
    double vmin, vmax;          //!< Range in speed units
    double bin_width;           //!< Width of the bins
    unsigned long *counts;      //!< Counts per bin
    unsigned long total_counts; //!< Total samples added
} VelHist;

/**
 * @brief MSD struct
 */
struct MSD{
    size_t frame;
    double *cor;
    size_t *count;
    struct Vec3D *prev;
    struct Vec3D *store;
};


#endif /* TYPES_MD_H_ */