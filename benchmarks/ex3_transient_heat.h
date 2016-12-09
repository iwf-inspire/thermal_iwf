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

#ifndef TRANSIENT_H_
#define TRANSIENT_H_

double solution_transient(double r2, double t, double alpha);
void transient_print_errors3D(particle* particles, double t, unsigned int Ntot);
void transient_heat_perform_cpu3D    (particle* particles, laplacian lap, double dt, unsigned int step, unsigned int Ntot);
void transient_heat_perform_rk4_cpu3D(particle* particles, laplacian lap, double dt, unsigned int step, unsigned int Ntot);

particle* transient_init3D(unsigned int N, double hdx, bool randomized);

#endif /* TRANSIENT_H_ */
