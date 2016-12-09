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
