#include "../methods/sph.h"
#include "../methods/pse.h"
#include "../methods/rkpm.h"
#include "../methods/cspm.h"
#include "../methods/nmfs.h"
#include "../methods/all_methods.h"

#include "../benchmarks/ex1_ex2_reconstruction.h"
#include "../benchmarks/ex3_transient_heat.h"
#include "../benchmarks/ex4_steady_heat.h"
#include "../benchmarks/ex5_eng_model.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#ifndef ENG_HEAT_MODEL_H_
#define ENG_HEAT_MODEL_H_

void run_eng_model_3D(METHOD method);

#endif /* ENG_HEAT_MODEL_H_ */
