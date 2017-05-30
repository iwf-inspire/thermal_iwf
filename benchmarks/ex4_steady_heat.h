#include "../structures/singleton_geometry.h"
#include "../structures/particle.h"

#include "../methods/sph.h"
#include "../methods/pse.h"
#include "../methods/rkpm.h"
#include "../methods/cspm.h"
#include "../methods/nmfs.h"
#include "../methods/all_methods.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef STEADY_H_
#define STEADY_H_

void steady_print_errors_3D(particle* particles, unsigned int N, unsigned int Nbnd);
void perform_steady_heat_euler_3D(particle* particles, METHOD method, double dt, unsigned int step, unsigned int N, unsigned int Nbnd);
void perform_steady_heat_rk4_3D  (particle* particles, METHOD method, double dt, unsigned int step, unsigned int N, unsigned int Nbnd);

particle *steady_init3D(unsigned int N, unsigned int Nbnd, double hdx, bool RANDOM);

#endif /* STEADY_H_ */
