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

#ifndef NMFS_H_
#define NMFS_H_

void nmfs_compute_correction_terms(particle* particles);
void perform_nmfs(particle* particles);
void nmfs_init(unsigned int num_part, kernel kernel_fun);
void nmfs_wipe_out();

#endif /* NMFS_H_ */
