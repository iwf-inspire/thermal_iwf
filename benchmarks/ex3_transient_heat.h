#include "../structures/singleton_geometry.h"
#include "../structures/particle.h"
#include "../methods/all_methods.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef TRANSIENT_H_
#define TRANSIENT_H_

void transient_print_errors_3D(particle* particles, double t, unsigned int Ntot);
void perform_transient_heat_euler_3D(particle* particles, METHOD method, double dt, unsigned int step, unsigned int Ntot);
void perform_transient_heat_rk4_3D  (particle* particles, METHOD method, double dt, unsigned int step, unsigned int Ntot);

particle* transient_init3D(unsigned int N, double hdx, bool RANDOM);

#endif /* TRANSIENT_H_ */
