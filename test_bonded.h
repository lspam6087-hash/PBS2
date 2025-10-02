/* test_bonded.h
 * ----------------------------------------------------------------------------
 * Finite-difference test for bonded forces (single molecule).
 * Declares the public test function to be called from main().
 * ----------------------------------------------------------------------------
 */

#ifndef TEST_BONDED_H_
#define TEST_BONDED_H_

#include "structs.h"  // struct Parameters, struct Vectors

/* 
 * Finite-difference test for bonded forces (single molecule).
 * Declares the public test function to be called from main().
 */
void run_bonded_fd_test(struct Parameters *p_global, struct Vectors *v);

#endif /* TEST_BONDED_H_ */