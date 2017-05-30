#include "../structures/particle.h"
#include "../neighboring/cell_list.h"
#include "../base/kernels.h"

#include <float.h>
#include <assert.h>
#include <stdbool.h>
#include <stdio.h>

#ifndef ALL_METHODS_H_
#define ALL_METHODS_H_

typedef enum METHOD {
	SPH,
	PSE,
	RKPM,
	CSPM,
	NMFS
} METHOD;

void initialize_meshfree_method(particle* particles, METHOD method, kernel kernel_func, unsigned int Ntot);
void correction_terms_meshfree_method(particle* particles, METHOD method);
void laplacian_meshfree_method(particle* particles, METHOD method);

void methods_wipe_out();

#endif /* ALL_METHODS_H_ */
