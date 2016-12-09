#include "../neighboring/cell_list.h"
#include "../structures/particle.h"
#include "../structures/singleton_geometry.h"

#include <string.h>
#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../base/kernels.h"

#ifndef SPH_H_
#define SPH_H_

void perform_sph(particle *particles);
void perform_sph_brookshaw(particle *particles);
void sph_set_kernel(kernel kernel_fun);
void sph_init(kernel kernel_fun);

#endif /* SPH_H_ */
