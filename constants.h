#ifndef CONSTANTS_H_
#define CONSTANTS_H_

/* Define constants for the entire program here */

#define PI 3.141592653589
#define SEED 13

// Internal MD units (Kelvin, Å, amu)
#define k_b 3.19e5    // bond spring constant, in K/Å^2
#define r0  1.54      // bond equilibrium length, in Å

// Optional: keep SI constants for conversions
#define k_B 1.380649e-23 // J/K
#define k_b_SI 440.427031   // N/m   (same constant in SI units)

//#define k_B 1.380649 //Boltzmann constant in J/K
//#define k_b 440.427031 //Bond force constant in N/m

#define k_theta 8.63e-19 //Angle force constant in J/rad^2
#define theta0 1.989675347 //Equilibrium angle in rad (114 degrees)

#endif /* CONSTANTS_H_ */
