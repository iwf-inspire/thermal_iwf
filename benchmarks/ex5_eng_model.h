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

// density (kg/m3)
#define RHO_STL 7850.0
#define RHO_ALU 2770.0
#define RHO_PUR 2170.0
// heat conductivity (W/m.k)
#define KA_STL  175.0
#define KA_ALU  150.0
#define KA_PUR  0.100
// specific heat capacity (J/kg.k)
#define CP_STL  434.0
#define CP_ALU  875.0
#define CP_PUR  1010.0

#ifndef ENG_THERMAL_H_
#define ENG_THERMAL_H_

double perform_eng_heat_euler_3D(particle* particles, METHOD method, double dt, unsigned int step);
double perform_eng_heat_rk4_3D  (particle* particles, METHOD method, double dt, unsigned int step);

particle* eng_init3D(double hdy, unsigned int Nbnd, double *spacing);

#endif /* ENG_THERMAL_H_ */
