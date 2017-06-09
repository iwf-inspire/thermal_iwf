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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <time.h>

#ifndef BENCHMARKS_H_
#define BENCHMARKS_H_

particle* setup_reconstruction_polynomial_3D(unsigned int N, double hdx, bool RANDOM);
particle* setup_reconstruction_trigonometric_3D(unsigned int N, double hdx, bool RANDOM);

#endif /* BENCHMARKS_H_ */
