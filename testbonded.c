#ifdef RUN_TEST_BONDED_FD
#include <math.h>
#include <stdio.h>
#include <string.h>

// Utilidad: evalúa SOLO energía enlazada (sin acumular fuerzas si do_forces=0)
static double eval_bonded(struct Parameters *P, struct Vectors *V, int do_forces)
{
    double E = 0.0;
    if (do_forces) {
        // limpia fuerzas
        memset(V->f, 0, P->num_part * sizeof(*V->f));
    }
    // IMPORTANTE: tus funciones deberían devolver la energía y
    // ACUMULAR fuerzas en V->f si do_forces=1. Si no, adapta.
    E += calculate_forces_bond(P, V);      // ya tienes implementado
    E += calculate_forces_angle(P, V);     // ya tienes el esqueleto
    E += calculate_forces_dihedral(P, V);  // lo estás añadiendo
    if (!do_forces) {
        // Si tus funciones siempre escriben fuerzas, guarda y restaura:
        // (opcional si ya controlas do_forces dentro)
    }
    return E;
}

static void make_single_molecule_ijkl(struct Parameters *P, struct Vectors *V)
{
    // 4 átomos: i(0), j(1), k(2), l(3)
    P->num_part = 4;
    P->L.x = P->L.y = P->L.z = 1000.0; // grande para evitar MIC accidental
    // Colocación no degenerada (activa bond, angle y dihedral):
    V->r[0] = (struct Vec3D){ 0.0, 0.0, 0.0};   // i
    V->r[1] = (struct Vec3D){ 1.5, 0.0, 0.0};   // j
    V->r[2] = (struct Vec3D){ 2.5, 0.7, 0.0};   // k
    V->r[3] = (struct Vec3D){ 3.2, 0.8, 1.0};   // l

    // Define bonds: (i-j), (j-k), (k-l)
    V->num_bonds = 3;
    V->bonds[0] = (struct Bond){ .i=0, .j=1, .k= P->kbond, .r0=P->r0bond_ij }; // adapta a tus campos
    V->bonds[1] = (struct Bond){ .i=1, .j=2, .k= P->kbond, .r0=P->r0bond_jk };
    V->bonds[2] = (struct Bond){ .i=2, .j=3, .k= P->kbond, .r0=P->r0bond_kl };

    // Angle: (i-j-k) y (j-k-l)
    V->num_angles = 2;
    V->angles[0] = (struct Angle){ .i=0, .j=1, .k=2, .k_th=P->kangle, .th0=P->th0_ijk };
    V->angles[1] = (struct Angle){ .i=1, .j=2, .k=3, .k_th=P->kangle, .th0=P->th0_jkl };

    // Dihedral: (i-j-k-l)
    V->num_dihedrals = 1;
    V->dihedrals[0] = (struct Dihedral){
        .i=0, .j=1, .k=2, .l=3,
        // ejemplo cos-series: V(φ)=∑_n k_n (1+cos(nφ-δ_n))
        .n=3, .k1=P->kdih1, .k2=P->kdih2, .k3=P->kdih3,
        .delta1=P->dih_delta1, .delta2=P->dih_delta2, .delta3=P->dih_delta3
    };

    // Desactiva no-enlazadas:
    // Si tienes lista de pares, pon num_pairs=0 o un flag para omitir no-bonded.
}

void run_test_bonded_fd(struct Parameters *P, struct Vectors *V)
{
    make_single_molecule_ijkl(P, V);

    // Fuerza analítica
    memset(V->f, 0, P->num_part * sizeof(*V->f));
    double E0 = eval_bonded(P, V, /*do_forces=*/1);

    const double h = 1e-6; // paso FD
    double max_abs = 0.0, max_rel = 0.0;
    size_t maxi=0; int maxc=0;

    // backup posiciones
    struct Vec3D r_backup[4];
    for (size_t a=0; a<4; ++a) r_backup[a] = V->r[a];

    for (size_t a=0; a<4; ++a) {
        for (int c=0; c<3; ++c) {
            // r+ = r; r+[a].c += h
            V->r[a] = r_backup[a];
            if (c==0) V->r[a].x += h; else if (c==1) V->r[a].y += h; else V->r[a].z += h;
            double Ep = eval_bonded(P, V, /*do_forces=*/0);

            // r- = r; r-[a].c -= h
            V->r[a] = r_backup[a];
            if (c==0) V->r[a].x -= h; else if (c==1) V->r[a].y -= h; else V->r[a].z -= h;
            double Em = eval_bonded(P, V, /*do_forces=*/0);

            // restaura
            V->r[a] = r_backup[a];

            double Fnum = -(Ep - Em) / (2.0*h);
            double Fana = (c==0)? V->f[a].x : (c==1? V->f[a].y : V->f[a].z);
            double abs_err = fabs(Fana - Fnum);
            double rel_err = abs_err / fmax(1.0, fabs(Fnum));

            if (abs_err > max_abs) { max_abs = abs_err; maxi=a; maxc=c; }
            if (rel_err > max_rel) { max_rel = rel_err; }
        }
    }

    printf("[BONDED-FD] E=%.15e  max_abs_err=%.3e  max_rel_err=%.3e  (atom %zu, comp %d)\n",
           E0, max_abs, max_rel, maxi, maxc);
}
#endif
