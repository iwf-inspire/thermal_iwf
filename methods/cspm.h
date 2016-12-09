#include "../neighboring/cell_list.h"
#include "../structures/particle.h"
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

#ifndef CSPM_H_
#define CSPM_H_

void cspm_compute_correction_terms(particle* particles);
void perform_cspm(particle* particles);
void cspm_set_kernel(kernel kernel_fun);
void cspm_init(unsigned int num_part, kernel kernel_fun);
void cspm_wipe_out();

#endif /* CSPM_H_ */
