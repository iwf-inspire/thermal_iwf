#include "../methods/all_methods.h"
#include "../benchmarks/ex1_ex2_reconstruction.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#ifndef RECONSTRUCTION_H_
#define RECONSTRUCTION_H_

void run_reconstrunction_trigonometric_3D(METHOD method, unsigned int N, bool RANDOM);
void run_reconstrunction_polynomial_3D(METHOD method, unsigned int N, bool RANDOM);

#endif /* RECONSTRUCTION_H_ */
