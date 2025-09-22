#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "constants.h"
#include "structs.h"

// Set the parameters of this simulation
void set_parameters(struct Parameters *p_parameters)
{
// The parameters first 5 parameters are only used for demonstration puprposes
  p_parameters->kT = 1.0;                                   //thermal energy
  p_parameters->mass = 1.0;                                 //mass of a particle
  p_parameters->epsilon = 1.0;                              //LJ interaction strength
  p_parameters->sigma = 1.0;                                //LJ particle diameter

// The parameters below control core functionalities of the code, but many values will need to be changed
  p_parameters->num_part = 2000;                            //number of particles
  p_parameters->num_dt_steps = 2000;                        //number of time steps
  p_parameters->exclude_12_nb = 1;                          // 1-2 connected atoms exluded from non-bonded interactions 
  p_parameters->exclude_13_nb = 1;                          // 1-3 connected atoms exluded from non-bonded interactions    
  p_parameters->dt = 0.01;                                  //integration time step
  p_parameters->L = (struct Vec3D){14.938, 14.938, 14.938}; //box size
    p_parameters->r_cut = 2.5;                              //cut-off distance used for neigbor list
  p_parameters->r_shell = 0.4;                              //shell thickness for neighbor list
  p_parameters->num_dt_pdb = 500;                           //number of time steps in between pdb outputs
  strcpy(p_parameters->filename_pdb, "trajectories");       //filename (without extension) for pdb file
  p_parameters->rescale_output = 1;                         //factor used to rescale output lengthscale (Most visualisation programs identify bonds based on distances of order 1)
  p_parameters->load_restart = 0;                           //if equal 1 restart file is loaded
  strcpy(p_parameters->restart_in_filename, "restart.dat"); //filename for loaded restart file
  p_parameters->num_dt_restart = 1000;                      // number of time steps between saves
  strcpy(p_parameters->restart_out_filename, "restart.dat");//filename for saved restart file


  if (p_parameters->r_cut > p_parameters->L.x / 2.0)
    fprintf(stderr, "Warning! r_cut > Lx/2");
  if (p_parameters->r_cut > p_parameters->L.y / 2.0)
    fprintf(stderr, "Warning! r_cut > Ly/2");
  if (p_parameters->r_cut > p_parameters->L.z / 2.0)
    fprintf(stderr, "Warning! r_cut > Lz/2");
}
