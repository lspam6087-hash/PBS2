#include <stdio.h>
#include <stdlib.h>
#include "structs.h"
#include "nbrlist.h"

// Allocate the arrays in 'vectors' needed to store information of all particles
void alloc_vectors(struct Vectors *p_vectors, struct Parameters *p_parameters, struct MSD *p_msd)
{
    size_t sz = p_parameters->num_part;
    p_vectors->size = sz;
    p_vectors->type = (int *)malloc(sz * sizeof(int));
    p_vectors->r = (struct Vec3D *)malloc(sz * sizeof(struct Vec3D));
    p_vectors->dr = (struct Vec3D *)malloc(sz * sizeof(struct Vec3D));
    p_vectors->v = (struct Vec3D *)malloc(sz * sizeof(struct Vec3D));
    p_vectors->f = (struct Vec3D *)malloc(sz * sizeof(struct Vec3D));
    p_vectors->bonds = NULL;
    p_vectors->num_bonds = 0;
    p_vectors->angles = NULL;
    p_vectors->num_angles = 0;
    p_vectors->dihedrals = NULL;
    p_vectors->num_dihedrals = 0;

    int molecules = sz / 4;

    p_msd->cor   = malloc(p_parameters->ncor * sizeof(double));
    p_msd->count = malloc(p_parameters->ncor * sizeof(size_t));
    p_msd->prev  = malloc(molecules * sizeof(struct Vec3D));
    p_msd->store = malloc(molecules * p_parameters->ncor * sizeof(struct Vec3D));
}

// Free the arrays in 'vectors'
void free_vectors(struct Vectors *p_vectors, struct MSD *p_msd)
{
    free(p_vectors->type);
    p_vectors->type = NULL;
    free(p_vectors->r);
    p_vectors->r = NULL;
    free(p_vectors->dr);
    p_vectors->dr = NULL;
    free(p_vectors->v);
    p_vectors->v = NULL;
    free(p_vectors->f);
    p_vectors->f = NULL;
    free(p_vectors->bonds);
    p_vectors->bonds = NULL;
    free(p_vectors->angles);
    p_vectors->angles = NULL;
    free(p_vectors->dihedrals);
    p_vectors->dihedrals = NULL;
    p_vectors->size = 0;
    p_vectors->num_bonds = 0;
    p_vectors->num_angles = 0;
    p_vectors->num_dihedrals = 0;

    free(p_msd->cor);
    free(p_msd->count);
    free(p_msd->prev);
    free(p_msd->store);
}

// Allocate all variables needed in the MD simulation
void alloc_memory(struct Parameters *p_parameters, struct Vectors *p_vectors, struct Nbrlist *p_nbrlist, struct MSD *p_msd)
{    
    alloc_vectors(p_vectors, p_parameters, p_msd);
    alloc_nbrlist(p_parameters, p_nbrlist);
}

// Free the memory allocated by alloc_memory
void free_memory(struct Vectors *p_vectors, struct Nbrlist *p_nbrlist, struct MSD *p_msd)
{
    free_vectors(p_vectors, p_msd);
    free_nbrlist(p_nbrlist);
}
