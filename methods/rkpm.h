//Copyright ETH Zurich, IWF

//This file is part of iwf_thermal.

//iwf_thermal is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.

//iwf_thermal is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.

//You should have received a copy of the GNU General Public License
//along with iwf_thermal.  If not, see <http://www.gnu.org/licenses/>.

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
void rkpm_init(unsigned int num_part, kernel kernel_fun);
void rkpm_wipe_out();

#endif /* RKPM_H_ */
