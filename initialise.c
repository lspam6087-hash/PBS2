#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constants.h"
#include "structs.h"
#include "random.h"
#include "initialise.h"

// This function initializes the particle types.
// Particle type 0 is assigned to all particles for now, but this could be modified
// later to initialize different types (e.g., methane and ethane in a binary mixture).
void initialise_types(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
    /// \todo (done) Initialize particle types in the vectors.type array
    for (size_t i = 0; i < p_parameters->num_part; i+=4)
    {
        p_vectors->type[i]   = 0; // Specify particle type (currently only one type)
        p_vectors->type[i+1] = 1; // Specify particle type (currently only one type)
        p_vectors->type[i+2] = 1; // Specify particle type (currently only one type)
        p_vectors->type[i+3] = 0; // Specify particle type (currently only one type)
    }
}

// This function initializes the bond connectivity between particles.
/// \todo Students need to specify bonds between pairs of particles, which are stored as bond[i].i and bond[i].j.
// This will be important for handling bonded interactions in the simulation.
void initialise_bond_connectivity(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
    size_t num_bonds = p_parameters->num_part * 3/4;  // Currently, no bonds are set up.
    int N = p_parameters->num_part;
    struct Bond *bonds = (struct Bond *)malloc(num_bonds * sizeof(struct Bond));

    /// \todo (urgent) Specify bonds between particles, i.e., bonds[i].i and bonds[i].j for bonded particle pairs.
    size_t x = 0;
    for (int n = 0; n < N; n++)
    {
        if ((n % 4) != 3)
        {
            bonds[x].i = n;
            bonds[x].j = n + 1;
            x++;
        }
        
    }
    
    p_vectors->num_bonds = num_bonds;
    p_vectors->bonds = bonds;
}

// This function initializes the molecular structure, including bonds, angles, and dihedrals,
// and computes 1-2, 1-3, and 1-4 bonded interactions for the neighbor list used in force calculations.
void initialise_structure(struct Parameters *p_parameters, struct Vectors *p_vectors, struct Nbrlist *p_nbrlist)
{
    initialise_bond_connectivity(p_parameters, p_vectors); // Initialize bonds
    
    struct Bond *bonds = p_vectors->bonds;
    size_t num_bonds = p_vectors->num_bonds;
    size_t num_part = p_parameters->num_part;
    size_t *cnt = (size_t *)calloc(num_part + 1, sizeof(size_t));
    size_t cnt_tot = 0;
    for (size_t i = 0; i < num_bonds; ++i)
    {
        ++cnt[bonds[i].i + 1];
        ++cnt[bonds[i].j + 1];
    }
    size_t *head12 = (size_t *)malloc((num_part + 1) * sizeof(size_t));
    head12[0] = 0;
    for (size_t i = 1; i <= num_part; ++i)
    {
        head12[i] = cnt[i] + head12[i - 1];
        cnt[i] = head12[i];
    }
    size_t *pairs12 = (size_t *)malloc(cnt[num_part] * sizeof(size_t));
    for (size_t i = 0; i < num_bonds; ++i)
    {
        pairs12[cnt[bonds[i].i]++] = bonds[i].j;
        pairs12[cnt[bonds[i].j]++] = bonds[i].i;
    }
    size_t num_angles = 0;
    for (size_t i = 1; i <= num_part; ++i)
    {
        size_t num_pairs = head12[i] - head12[i - 1];
        num_angles += (num_pairs * (num_pairs - 1)) / 2;
    }
    struct Angle *angles = (struct Angle *)malloc(num_angles * sizeof(struct Angle));
    size_t m = 0;
    for (size_t i = 0; i < num_part; ++i)
        for (size_t j = head12[i]; j <= head12[i + 1]; ++j)
            for (size_t k = j + 1; k < head12[i + 1]; ++k)
            {
                angles[m].i = pairs12[j];
                angles[m].j = i;
                angles[m].k = pairs12[k];
                ++m;
            }
    size_t num_dihdr = 0;
    for (size_t i = 0; i < num_bonds; ++i)
    {
        num_dihdr += (head12[bonds[i].i + 1] - head12[bonds[i].i] - 1) * (head12[bonds[i].j + 1] - head12[bonds[i].j] - 1);
    }
    struct Dihedral *dihedrals = (struct Dihedral *) malloc(num_dihdr * sizeof(struct Dihedral));
    size_t k = 0;
    for (size_t i = 0; i < num_bonds; ++i)
    {
        struct Dihedral dihdr;
        dihdr.j = bonds[i].i;
        dihdr.k = bonds[i].j;
        for (size_t j = head12[dihdr.j]; j < head12[dihdr.j + 1]; ++j)
        {
            if (pairs12[j] == dihdr.k)
                continue;
            dihdr.i = pairs12[j];
            for (size_t j = head12[dihdr.k]; j < head12[dihdr.k + 1]; ++j)
            {
                if (pairs12[j] != dihdr.j)
                {
                    dihdr.l = pairs12[j];
                    dihedrals[k++] = dihdr;
                }
            }
        }
    }
    if (p_parameters->exclude_13_nb)
    {
        for (size_t i = 0; i <= num_part; ++i)
            cnt[i] = 0;
        for (size_t i = 0; i < num_angles; ++i)
        {
            ++cnt[angles[i].i + 1];
            ++cnt[angles[i].k + 1];
        }
        size_t *head13 = (size_t *)malloc((num_part + 1) * sizeof(size_t));
        head13[0] = 0;
        for (size_t i = 1; i <= num_part; ++i)
        {
            head13[i] = cnt[i] + head13[i - 1];
            cnt[i] = head13[i];
        }
        size_t *pairs13 = (size_t *)malloc(cnt[num_part] * sizeof(size_t));
        for (size_t i = 0; i < num_angles; ++i)
        {
            pairs13[cnt[angles[i].i]++] = angles[i].k;
            pairs13[cnt[angles[i].k]++] = angles[i].i;
        }
        p_nbrlist->head13 = head13;
        p_nbrlist->pairs13 = pairs13;
    }

    for (size_t i = 0; i <= num_part; ++i)
        cnt[i] = 0;
    for (size_t i = 0; i < num_dihdr; ++i)
    {
        ++cnt[dihedrals[i].i + 1];
        ++cnt[dihedrals[i].l + 1];
    }
    size_t *head14 = (size_t *)malloc((num_part + 1) * sizeof(size_t));
    head14[0] = 0;
    for (size_t i = 1; i <= num_part; ++i)
    {
        head14[i] = cnt[i] + head14[i - 1];
        cnt[i] = head14[i];
    }
    size_t *pairs14 = (size_t *)malloc(cnt[num_part] * sizeof(size_t));
    for (size_t i = 0; i < num_dihdr; ++i)
    {
        pairs14[cnt[dihedrals[i].i]++] = dihedrals[i].l;
        pairs14[cnt[dihedrals[i].l]++] = dihedrals[i].i;
    }
    p_nbrlist->head14 = head14;
    p_nbrlist->pairs14 = pairs14;
    p_vectors->num_angles = num_angles;
    p_vectors->angles = angles;
    p_vectors->num_dihedrals = num_dihdr;
    p_vectors->dihedrals = dihedrals;
    if (p_parameters->exclude_12_nb)
    {
        p_nbrlist->head12 = head12;
        p_nbrlist->pairs12 = pairs12;
    }
    else
    {
        free(head12);
        free(pairs12);
    }
    free(cnt);
}


// This function initializes the simulation by calling subroutines to initialize
// particle types, positions, velocities, and bond connectivity.
void initialise(struct Parameters *p_parameters, struct Vectors *p_vectors, struct Nbrlist *p_nbrlist, size_t *p_step, double *p_time)
{
    initialise_types(p_parameters, p_vectors);  // Initialize particle types
    initialise_structure(p_parameters, p_vectors, p_nbrlist);  // Initialize structure (bonds, angles, dihedrals)
    srand(SEED);  // Seed random number generator
    initialise_positions(p_parameters, p_vectors);  // Initialize particle positions
    initialise_velocities(p_parameters, p_vectors);  // Initialize particle velocities
    *p_step = 0;   // Initialize step to zero
    *p_time = 0.0; // Initialize time to zero
    return;
}


// This function initializes particle positions on a cubic lattice.
// Particles are placed in a grid with spacing based on the number of particles and the box dimensions.
//Build a n-butane molecule with 4 particles in already the right geometry
void build_molecule(struct Parameters *p_parameters, struct Vec3D molecule_center, struct Vectors *p_vectors, int ipart)
{
    //place particle 1 (CH3)
    p_vectors->r[ipart].x = molecule_center.x - r0;
    p_vectors->r[ipart].y = molecule_center.y;
    p_vectors->r[ipart].z = molecule_center.z;

    //place particle 2 (CH2)
    p_vectors->r[ipart + 1].x = molecule_center.x;
    p_vectors->r[ipart + 1].y = molecule_center.y;
    p_vectors->r[ipart + 1].z = molecule_center.z;

    //place particle 3 (CH2)
    double dx = r0 * cos(theta0);  
    double dy = r0 * sin(theta0);       
    p_vectors->r[ipart + 2].x = p_vectors->r[ipart + 1].x - dx;
    p_vectors->r[ipart + 2].y = p_vectors->r[ipart + 1].y + dy;
    p_vectors->r[ipart + 2].z = molecule_center.z;

    //place particle 4 (CH3)
    p_vectors->r[ipart + 3].x = p_vectors->r[ipart + 2].x + r0;
    p_vectors->r[ipart + 3].y = p_vectors->r[ipart + 2].y;
    p_vectors->r[ipart + 3].z = molecule_center.z;
}


// This function initializes particle positions on a cubic lattice.
// Particles are placed in a grid with spacing based on the number of particles and the box dimensions.
void initialise_positions(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
    struct Vec3D dr;  // Displacement vector for positioning particles
    struct Index3D n; // Number of grid cells along each axis
    double dl;        // Lattice spacing
    int ipart = 0;    // Particle index
    int imol = 0;     //molecule index

    int num_mol = p_parameters->num_part / 4;
    double bondlength = r0 + 0.0;

    //If there is only 1 molecule: center in box
    if (num_mol == 1) 
    {
        // centre of the box
        struct Vec3D molecule_center;
        molecule_center.x = 0.5 * p_parameters->L.x;
        molecule_center.y = 0.5 * p_parameters->L.y;
        molecule_center.z = 0.5 * p_parameters->L.z;

        // place 4 atoms in a straight x-direction chain centred at box centre
        build_molecule(p_parameters, molecule_center, p_vectors, 0);
        return;
    }

    // Calculate lattice spacing based on particle number and box dimensions
    //dl = pow(p_parameters->L.x * p_parameters->L.y * p_parameters->L.z / ((double)num_mol), 1.0 / 3.0);
    //n.i = (int)ceil(p_parameters->L.x / dl);
    //n.j = (int)ceil(p_parameters->L.y / dl);
    //n.k = (int)ceil(p_parameters->L.z / dl);

    n.i = (int) floor(p_parameters->L.x / 7.0);
    if (n.i < 1) 
        n.i = 1;

    int mol_per_layer = (int) ceil((double) num_mol / n.i);
    n.j = (int) floor(sqrt(mol_per_layer));
    if (n.j < 1) 
        n.j = 1;
    n.k = (int) ceil((double) mol_per_layer / n.j);

    dr.x = p_parameters->L.x / (double)n.i;
    dr.y = p_parameters->L.y / (double)n.j;
    dr.z = p_parameters->L.z / (double)n.k;
    ipart = 0;
    imol = 0;
    for (size_t i = 0; i < n.i; ++i)
    {
        for (size_t j = 0; j < n.j; ++j)
        {
            for (size_t k = 0; k < n.k; ++k)
            {
                if (imol >= num_mol)
                    break;
                
                //Molecule center
                struct Vec3D molecule_center;
                molecule_center.x = (i + 0.5) * dr.x;
                molecule_center.y = (j + 0.5) * dr.y;
                molecule_center.z = (k + 0.5) * dr.z;
                
                build_molecule(p_parameters, molecule_center, p_vectors, ipart);

                ipart += 4;
                imol++;
            }
        if (imol >= num_mol)
            break;
        }
    if (imol >= num_mol)
        break;
    }
}

// This function initializes the velocities of particles based on the Maxwell-Boltzmann distribution.
// The total momentum is also removed to ensure zero total momentum (important for stability).
void initialise_velocities(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
    struct Vec3D sumv = {0.0, 0.0, 0.0};  // Total velocity (to remove later)

    // Assign random velocities to each particle
    for (size_t i = 0; i < p_parameters->num_part; i++)
    {
        //Use specific mass per particleelocit
        double mass = p_parameters->mass[p_vectors->type[i]];
        double sqrtktm = sqrt(p_parameters->kT / mass);

        p_vectors->v[i].x = sqrtktm * gauss();
        p_vectors->v[i].y = sqrtktm * gauss();
        p_vectors->v[i].z = sqrtktm * gauss();
        sumv.x += p_vectors->v[i].x;
        sumv.y += p_vectors->v[i].y;
        sumv.z += p_vectors->v[i].z;
    }

    //calculate center of mass velocity
    struct Vec3D momentum = {0.0, 0.0, 0.0};
    double mtot = 0.0;

    for (size_t i = 0; i < p_parameters->num_part; i++)
    {
        double m = p_parameters->mass[p_vectors->type[i]];
        momentum.x += m * p_vectors->v[i].x;
        momentum.y += m * p_vectors->v[i].y;
        momentum.z += m * p_vectors->v[i].z;
        mtot += m;
    }

    struct Vec3D Vcm = {momentum.x / mtot, momentum.y / mtot, momentum.z / mtot};

    //remove Vcm from each particle
    for (size_t i = 0; i < p_parameters->num_part; i++)
    {
        p_vectors->v[i].x -= Vcm.x;
        p_vectors->v[i].y -= Vcm.y;
        p_vectors->v[i].z -= Vcm.z;
    }
}

double box_length_from_density_A(size_t Nmol, double rho_kg_m3)
{
    const double NA = 6.02214076e23; // 1/mol
    const double M_g_mol = 58.123;   // mass molar butane g/mol
    double mass_total_kg = (Nmol * (M_g_mol / 1000.0)) / NA;
    double V_m3 = mass_total_kg / rho_kg_m3;
    double L_m  = cbrt(V_m3) * 1e10;
    return L_m;
}

void initialise_bonds(struct Parameters *p, struct Vectors *v)
{
    struct Nbrlist dummy = (struct Nbrlist){0}; // NB off en el test
    initialise_structure(p, v, &dummy);
}