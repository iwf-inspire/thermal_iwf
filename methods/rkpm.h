#include "../neighboring/cell_list.h"
#include "../structures/singleton_geometry.h"
#include "../structures/particle.h"

#include <string.h>
#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../base/kernels.h"
#include "../base/utils.h"

#ifndef RKPM_H_
#define RKPM_H_

void rkpm_compute_correction_terms(particle* particles);
void perform_rkpm(particle* particles);
void rkpm_set_kernel(kernel kernel_fun);
void rkpm_init(unsigned int num_part, kernel kernel_fun);
void rkpm_wipe_out();

#endif /* RKPM_H_ */
