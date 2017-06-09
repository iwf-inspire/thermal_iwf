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

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include "../structures/singleton_geometry.h"
#include "../structures/particle.h"

#ifndef CELL_LIST_H_
#define CELL_LIST_H_

void grid_update_grid_geom(particle* particles);		//call this before grid sort!
void grid_sort(particle* particles);
void grid_restore(particle* particles);
void grid_unhash(unsigned int idx, int *i, int *j, int *k);

#endif /* CELL_LIST_H_ */
