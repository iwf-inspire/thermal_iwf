#include "../methods/sph.h"
#include "../methods/pse.h"
#include "../methods/rkpm.h"
#include "../methods/cspm.h"
#include "../methods/nmfs.h"
#include "../methods/all_methods.h"

#include "../benchmarks/ex1_ex2_reconstruction.h"
#include "../benchmarks/ex3_transient_heat.h"
#include "../benchmarks/ex4_steady_heat.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#ifndef THERMAL_TEST_H_
#define THERMAL_TEST_H_

void run_transient_heat_3D(METHOD method, unsigned int N, bool RANDOM);
void run_steady_heat_3D   (METHOD method, unsigned int N, bool RANDOM);

#endif /* THERMAL_TEST_H_ */
