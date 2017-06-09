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

#include "../methods/all_methods.h"
#include "../benchmarks/ex1_ex2_reconstruction.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#ifndef RECONSTRUCTION_H_
#define RECONSTRUCTION_H_

void run_reconstrunction_trigonometric_3D(METHOD method, unsigned int N, bool RANDOM);
void run_reconstrunction_polynomial_3D(METHOD method, unsigned int N, bool RANDOM);

#endif /* RECONSTRUCTION_H_ */
