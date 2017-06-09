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

#ifndef STEADY_H_
#define STEADY_H_

void steady_print_errors_3D(particle* particles, unsigned int N, unsigned int Nbnd);
void perform_steady_heat_euler_3D(particle* particles, METHOD method, double dt, unsigned int step, unsigned int N, unsigned int Nbnd);
void perform_steady_heat_rk4_3D  (particle* particles, METHOD method, double dt, unsigned int step, unsigned int N, unsigned int Nbnd);

particle *steady_init3D(unsigned int N, unsigned int Nbnd, double hdx, bool RANDOM);

#endif /* STEADY_H_ */
