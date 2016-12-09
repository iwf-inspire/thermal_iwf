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
