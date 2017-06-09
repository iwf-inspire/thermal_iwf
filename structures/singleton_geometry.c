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

#include "singleton_geometry.h"

singleton_geometry * get_singleton_geometry() {
    static singleton_geometry *instance = NULL;
    if(instance == NULL) {
    	instance = (singleton_geometry*) malloc(sizeof(*instance));

    	instance->bbmin_grid_x = 0.;
    	instance->bbmin_grid_y = 0.;
    	instance->bbmin_grid_z = 0.;
    	instance->bbmax_grid_x = 0.;
    	instance->bbmax_grid_y = 0.;
    	instance->bbmax_grid_z = 0.;

    	instance->bbmin_part_x = 0.;
    	instance->bbmin_part_y = 0.;
    	instance->bbmin_part_z = 0.;
    	instance->bbmax_part_x = 0.;
    	instance->bbmax_part_y = 0.;
    	instance->bbmax_part_z = 0.;

    	instance->cell_start = NULL;
    	instance->cell_end   = NULL;

    	instance->dx = 0.;

    	instance->lx = 0.;
    	instance->ly = 0.;
    	instance->lz = 0.;

    	instance->num_cell = 0;

    	instance->n  = 0;
    	instance->nx = 0;
    	instance->ny = 0;
    	instance->nz = 0;
    }
    return instance;
}
