#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "constants.h"
#include "structs.h"

// Set the parameters of this simulation
void set_parameters(struct Parameters *p_parameters)
{
p_parameters->kT = 298.0; // si sigues usando energía reducida

// Mass of both types of particles
p_parameters->mass[0] = 15.0;  // CH3
p_parameters->mass[1] = 14.0;  // CH2

// Parameters LJ table
p_parameters->sigma[0]   = 3.75; // CH3
p_parameters->sigma[1]   = 3.95; // CH2
p_parameters->epsilon[0] = 98.0; // CH3 (K, en ε/kB)
p_parameters->epsilon[1] = 46.0; // CH2 (K, en ε/kB)

// The parameters below control core functionalities of the code, but many values will need to be changed
  p_parameters->num_part = 560;                            //number of particles
  p_parameters->num_dt_steps = 20000;                        //number of time steps
  p_parameters->exclude_12_nb = 1;                          // 1-2 connected atoms exluded from non-bonded interactions 
  p_parameters->exclude_13_nb = 1;                          // 1-3 connected atoms exluded from non-bonded interactions 
  p_parameters->exclude_14_nb = 1;                          // 1-4 connected atoms exluded from non-bonded interactions                                 
  p_parameters->dt = 0.00001;                               //integration time step
  p_parameters->L = (struct Vec3D){28.2, 28.2, 28.2};       //box size
  p_parameters->r_cut = 14;                                 //cut-off distance used for neigbor list
  p_parameters->r_shell = 1.5;                              //shell thickness for neighbor list
  p_parameters->num_dt_pdb = 500;                           //number of time steps in between pdb outputs
  strcpy(p_parameters->filename_pdb, "trajectories");       //filename (without extension) for pdb file
  p_parameters->rescale_output = 1;                         //factor used to rescale output lengthscale (Most visualisation programs identify bonds based on distances of order 1)
  p_parameters->load_restart = 0;                           //if equal 1 restart file is loaded
  strcpy(p_parameters->restart_in_filename, "restart.dat"); //filename for loaded restart file
  p_parameters->num_dt_restart = 1000;                      // number of time steps between saves
  strcpy(p_parameters->restart_out_filename, "restart.dat");//filename for saved restart file
  strcpy(p_parameters->filename_diag, "diagnostics");       // filename (without extension) for diagnostics CSV LAURA
  p_parameters->nbins = 200;                                // set the number of bins
  p_parameters->sample_interval = 10;                       // sample the speeds on every interval
  p_parameters->write_interval = 200;                       // write histogram to disk every X stepsp_parameters->nbins = 200;
  p_parameters->hist_vmax = 20;                             //max velocity of the histogram
  strcpy(p_parameters->filename_hist, "velocity-hist.csv"); //filename of the histogram
  p_parameters->nbins_dih = 75;                             // number of bins for the dihedral histogram
  strcpy(p_parameters->filename_dih_hist, "velocity-dih-hist.csv"); //filename of the dihedral histogram
  p_parameters->nbins_dih           = 72;                    // number of bins dihedral histogram
  p_parameters->sample_interval_dih = 10;                   //sample interval dihedral
  p_parameters->write_interval_dih  = 500;                  // write interval dih
  strcpy(p_parameters->filename_dih_samples, "torsion_phi.dat");  // filename dihedral
  strcpy(p_parameters->filename_dih_hist,    "torsion_hist.csv"); // filename dihedral data
  strcpy(p_parameters->MSD_filename,         "MSD.csv");      // filename MSD data

  p_parameters->ncor = 2000; //number of samples
  p_parameters->num_dtau = 10; // interval for data storage MSD

  if (p_parameters->r_cut > p_parameters->L.x / 2.0)
    fprintf(stderr, "Warning! r_cut > Lx/2");
  if (p_parameters->r_cut > p_parameters->L.y / 2.0)
    fprintf(stderr, "Warning! r_cut > Ly/2");
  if (p_parameters->r_cut > p_parameters->L.z / 2.0)
    fprintf(stderr, "Warning! r_cut > Lz/2");
}
