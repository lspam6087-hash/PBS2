// test_forces.c
// Test for non-bonded force: compare analytical and finite-difference force for a particle pair
// LAURA

#include <stdio.h>
#include <math.h>
#include "structs.h"
#include "forces.h"

// Simple Lennard-Jones potential for two particles at distance r
double lj_potential(double epsilon, double sigma, double r) {
    double sr6 = pow(sigma / r, 6);
    double sr12 = sr6 * sr6;
    return 4.0 * epsilon * (sr12 - sr6);
}

// Analytical force for two particles at distance r
double lj_force(double epsilon, double sigma, double r) {
    double sr6 = pow(sigma / r, 6);
    double sr12 = sr6 * sr6;
    return 24.0 * epsilon * (2.0 * sr12 - sr6) / r;
}

int main() {
    // Example parameters for CH3-CH2 interaction (replace with your actual values)
    // Values from the table (epsilon in K, sigma in Angstrom)
    double epsilon_CH3 = 98.0;
    double epsilon_CH2 = 46.0;
    double sigma_CH3 = 3.75;
    double sigma_CH2 = 3.95;

    // Lorentz-Berthelot mixing rule
    double epsilon = sqrt(epsilon_CH3 * epsilon_CH2); // Berthelot (geometric mean)
    double sigma = 0.5 * (sigma_CH3 + sigma_CH2);     // Lorentz (arithmetic mean)

    double r = 1.2; // Distance between particles
    double delta = 1e-5;

    // Potential at r and r+delta
    double U1 = lj_potential(epsilon, sigma, r);
    double U2 = lj_potential(epsilon, sigma, r + delta);

    // Numerical force (finite difference)
    double F_num = -(U2 - U1) / delta;
    // Analytical force
    double F_analytical = lj_force(epsilon, sigma, r);

    printf("LJ test for CH3-CH2 pair at r = %.5f\n", r);
    printf("Numerical force:   % .8f\n", F_num);
    printf("Analytical force:  % .8f\n", F_analytical);
    printf("Relative error:    %.2e\n", fabs(F_num - F_analytical) / fabs(F_analytical));

    return 0;
}
