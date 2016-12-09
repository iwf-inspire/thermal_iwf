#include "../structures/singleton_geometry.h"
#include "../structures/particle.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../methods/cspm.h"
#include "../methods/laplacian.h"
#include "../methods/nmfs.h"
#include "../methods/pse.h"
#include "../methods/rkpm.h"
#include "../methods/sph.h"

#ifndef STEADY_H_
#define STEADY_H_

double solution_steady(double Tamp, double Tbnd, double x, double y, double z);

void steady_print_errors3D(particle* particles, unsigned int N, unsigned int Nbnd);
void steady_heat_perform_test_cpu3D(particle* particles, laplacian lap, double dt, unsigned int step, unsigned int N);
void steady_heat_perform_cpu3D     (particle* particles, laplacian lap, double dt, unsigned int step, unsigned int N, unsigned int Nbnd);
void steady_heat_perform_rk4_cpu3D (particle* particles, laplacian lap, double dt, unsigned int step, unsigned int N, unsigned int Nbnd);
particle *steady_init3D(unsigned int N, unsigned int Nbnd, double hdx, bool randomized);
particle* steady_test_init3D(unsigned int N, double hdx, bool randomized);

#endif /* STEADY_H_ */
