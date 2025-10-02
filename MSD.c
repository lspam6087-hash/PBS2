#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "structs.h"
#include "Vec3D.h"
#include "MSD.h"

//Calculate COM
struct Vec3D com_molecule(struct Parameters *p_parameters, struct Vectors *p_vectors, int imolecule)
{
    struct Vec3D L = p_parameters->L;
    int molecule_index = 4 * imolecule;
    struct Vec3D com = v3(0.0, 0.0, 0.0);
    struct Vec3D *r = p_vectors->r;
    int type = p_vectors->type;
    double *mass = p_parameters->mass;
    double mtot = 0.0;
    
    //Look at first particle in molecule
    struct Vec3D r0 = r[molecule_index];
    double m0 = mass[p_vectors->type[molecule_index]];
    com = add(com, scl(m0, r0));
    mtot += m0;

    //Look at other atoms while applying minimum image convention
    for (int i = 1; i < 4; i++) 
    {
        int j = molecule_index + i;   
        int type = p_vectors->type[j];
        double mass1 = mass[type];

        struct Vec3D rij;
        rij.x = p_vectors->r[j].x - r0.x;
        rij.x = rij.x - L.x * floor(rij.x / L.x + 0.5);
        rij.y = p_vectors->r[j].y - r0.y;
        rij.y = rij.y - L.y * floor(rij.y / L.y + 0.5);
        rij.z = p_vectors->r[j].z - r0.z;
        rij.z = rij.z - L.z * floor(rij.z / L.z + 0.5);

        struct Vec3D rj = add(r0, rij);

        com = add(com, scl(mass1, rj));
        mtot += mass1;
    }
    return scl(1.0/mtot, com);
}


//Initialise MSD
void initialise_msd(struct Parameters *p_parameters, struct Vectors *p_vectors, struct MSD *p_msd)
{
    int molecules = p_parameters->num_part / 4;
    p_msd->frame = 0;

    //zero cor and count array
    for (size_t i = 0; i < p_parameters->ncor; i++) {
        p_msd->cor[i]   = 0.0;
        p_msd->count[i] = 0;
    }

    //store first COM positions
    for (int i = 0; i < molecules; i++) {
        struct Vec3D com = com_molecule(p_parameters, p_vectors, i);
        p_msd->prev[i] = com;
        p_msd->store[i] = com;
    }

    p_msd->frame++;
}

//Update MSD
void update_msd(struct Parameters *p_parameters, struct Vectors *p_vectors, struct MSD *p_msd)
{
    int molecules = p_parameters->num_part / 4;
    size_t maxcor = p_parameters->ncor - 1;
        if (p_msd->frame < p_parameters->ncor)
            maxcor = p_msd->frame;

    //current frame
    size_t curframe = p_msd->frame % p_parameters->ncor;

    //loop over all molecules
    for (int m = 0; m < molecules; m++) 
    {
        struct Vec3D com = com_molecule(p_parameters, p_vectors, m);

        //displacement of molecule with minimum image convention
        struct Vec3D d = sub(com, p_msd->prev[m]);
        d.x = d.x - p_parameters->L.x * floor(d.x / p_parameters->L.x + 0.5);
        d.y = d.y - p_parameters->L.y * floor(d.y / p_parameters->L.y + 0.5);
        d.z = d.z - p_parameters->L.z * floor(d.z / p_parameters->L.z + 0.5);

        //store COM
        struct Vec3D com_unfold = add(p_msd->prev[m], d);
        p_msd->store[curframe * molecules + m] = com_unfold;
        p_msd->prev[m] = com_unfold; //store previous for next step

        //correlate with older frames
        for (size_t icor = 0; icor <= maxcor; icor++) 
        {
            size_t icorframe = (curframe + p_parameters->ncor - icor) % p_parameters->ncor;
            struct Vec3D old = p_msd->store[icorframe * molecules + m];

            double dx = com_unfold.x - old.x;
            double dy = com_unfold.y - old.y;
            double dz = com_unfold.z - old.z;
            double dr2 = dx*dx + dy*dy + dz*dz;

            p_msd->cor[icor]   += dr2;
            p_msd->count[icor] += 1;
        }

    }

    // advance frame counter
    p_msd->frame++;
}