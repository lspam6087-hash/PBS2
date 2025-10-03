#ifndef TESTS_H
#define TESTS_H

#include "structs.h"
#include "nbrlist.h"

void run_pair_case(int type_i, int type_j, struct Parameters *p_global,
                   struct Vectors *v, struct Nbrlist *nl);
void run_bonded_fd_tests(struct Parameters *p_global,
                         struct Vectors *v,
                         struct Nbrlist *nl);

#endif