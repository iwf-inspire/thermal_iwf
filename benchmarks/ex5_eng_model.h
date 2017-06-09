#include "../structures/singleton_geometry.h"
#include "../structures/singleton_physics.h"
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

#ifndef ENG_THERMAL_H_
#define ENG_THERMAL_H_

double perform_eng_heat_euler_3D(particle* particles, METHOD method, double dt, unsigned int step);

particle* eng_init3D(double hdy, unsigned int Nbnd, double *spacing);

#endif /* ENG_THERMAL_H_ */
