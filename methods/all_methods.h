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
