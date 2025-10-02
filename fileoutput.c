#include <stdio.h>
#include <stdlib.h>
#include "constants.h"
#include "memory.h"
#include "structs.h"
#include "fileoutput.h"

// Write the particle positions to a pdf file
// The filename (without extension) is given by p_parameters->filename_pdb.
// If reset = 1 the data is written to the file deleting data it possibly contained.
// If reset = 0 the data is appended.
void record_trajectories_pdb(int reset, struct Parameters *p_parameters, struct Vectors *p_vectors, double time)
{
  FILE *fp_traj;
  char filename[1024];
  double rs = p_parameters->rescale_output;

  snprintf(filename, 1024, "%s%s", p_parameters->filename_pdb, ".pdb");
  if (reset == 1)
  {
    fp_traj = fopen(filename, "w");
  }
  else
  {
    fp_traj = fopen(filename, "a");
  }

  fprintf(fp_traj, "MODEL\n");
  fprintf(fp_traj, "REMARK TIME = %f\n", time);
  fprintf(fp_traj, "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-10s%-3s\n", rs*p_parameters->L.x, rs*p_parameters->L.y, rs*p_parameters->L.z, 90.0, 90.0, 90.0, "P 1", "1");
  for (size_t i = 0; i < p_parameters->num_part; i++)
  {
    fprintf(fp_traj, "HETATM%5u  C   UNK A   1    %8.3f%8.3f%8.3f  1.00  0.00           C\n", (unsigned int)i % 100000, rs*p_vectors->r[i].x, rs*p_vectors->r[i].y, rs*p_vectors->r[i].z);
  }
  fprintf(fp_traj, "ENDMDL\n");

  fclose(fp_traj);
}

// Write the particle positions to a xyz file
// The filename (without extension) is given by p_parameters->filename_xyz.
// If reset = 1 the data is written to the file deleting data it possibly contained.
// If reset = 0 the data is appended.
void record_trajectories_xyz(int reset, struct Parameters *p_parameters, struct Vectors *p_vectors, double time)
{
  FILE *fp_traj;
  char filename[1024];
  double rs = p_parameters->rescale_output;

  snprintf(filename, 1024, "%s%s", p_parameters->filename_xyz, ".xyz");
  if (reset == 1)
  {
    fp_traj = fopen(filename, "w");
  }
  else
  {
    fp_traj = fopen(filename, "a");
  }

  fprintf(fp_traj, "%lu\n", p_parameters->num_part);
  fprintf(fp_traj, "time = %f\n", time);
  struct Vec3D *r = p_vectors->r;
  for (size_t i = 0; i < p_parameters->num_part; i++)
  {
    fprintf(fp_traj, "  C        %10.5f %10.5f %10.5f\n", rs*r[i].x, rs*r[i].y, rs*r[i].z);
  }

  fclose(fp_traj);
}

// save arrays in vectors to binary file
void save_restart(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
  FILE* p_file = fopen( p_parameters->restart_out_filename, "wb");
  size_t num_part = p_parameters->num_part;
  size_t sz = num_part*sizeof(struct Vec3D);

  fwrite(&num_part, sizeof(size_t), 1, p_file);
  fwrite(p_vectors->r, sz, 1, p_file);
  fwrite(p_vectors->v ,sz, 1, p_file);
  fwrite(p_vectors->f, sz, 1, p_file);
  fclose(p_file);
}

// load arrays in vectors from binary file
void load_restart(struct Parameters *p_parameters, struct Vectors *p_vectors, struct MSD *p_msd)
{
  FILE* p_file = fopen( p_parameters->restart_in_filename, "rb" );
  size_t num_part;
  fread(&num_part, sizeof(size_t), 1, p_file);
  size_t sz = num_part*sizeof(struct Vec3D);
  alloc_vectors(p_vectors, p_parameters, p_msd);
  p_parameters->num_part = num_part;
  fread(p_vectors->r, sz, 1, p_file);
  fread(p_vectors->v, sz, 1, p_file);
  fread(p_vectors->f, sz, 1, p_file);
  fclose(p_file);
}

// Write diagnostic data (energies, temperature) to CSV file LAURA B1
void record_diagnostics_csv(int reset, struct Parameters *p_parameters, double time, double kin_energy, double pot_energy, double temperature)
{
    FILE *fp_diag;
    char filename[1024];
    snprintf(filename, 1024, "%s%s", p_parameters->filename_diag, ".csv");

    if (reset == 1) {
      // Print the header
      fp_diag = fopen(filename, "w");
      fprintf(fp_diag, "time,kinetic_energy,potential_energy,total_energy,temperature\n");
    } else {
      // Print data to the file
      fp_diag = fopen(filename, "a");
    }

    fprintf(fp_diag, "%f,%f,%f,%f,%f\n",
            time, kin_energy, pot_energy, kin_energy + pot_energy, temperature);

    fclose(fp_diag);
}

void record_msd_csv(struct Parameters *p_parameters, struct MSD *p_msd)
{
    FILE *fp_msd;
    char filename[1024];

    snprintf(filename, 1024, "%s", p_parameters->MSD_filename);
    fp_msd = fopen(filename, "w");
    if (!fp_msd) {
        printf("Error opening MSD file %s!\n", filename);
        return;
    }

    //print header
    fprintf(fp_msd, "# num_dtau=%zu dt=%.6e ncor=%zu\n", p_parameters->num_dtau, p_parameters->dt, p_parameters->ncor);
    fprintf(fp_msd, "time_lag,MSD,count\n");

    //output finilased msd
    for (size_t icor = 0; icor < p_parameters->ncor; icor++) {
        if (p_msd->count[icor] > 0) {
            double tau = icor * p_parameters->num_dtau * p_parameters->dt;
            double msd = p_msd->cor[icor] / (double)p_msd->count[icor];
            fprintf(fp_msd, "%.6e,%.6e,%zu\n", tau, msd, p_msd->count[icor]);
        }
    }
    fclose(fp_msd);
}

