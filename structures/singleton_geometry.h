#include <stdlib.h>

#ifndef SINGLETON_GEOMETRY_H_
#define SINGLETON_GEOMETRY_H_

typedef struct singleton_geometry {
	unsigned int n;		// number of particles

	unsigned int num_cell;		// number of cells
	unsigned int nx; 			// number of cells in x direction
	unsigned int ny;			// number of cells in y direction
	unsigned int nz;			// number of cells in z direction
	int *cell_start;			// points into particle vector for cell starts
	int *cell_end;				// points into particle vector for cell ends
	double dx;					// grid cell size
	double lx;					// grid extents
	double ly;
	double lz;

	double bbmin_part_x;		// bounding box for particles (tight bb)
	double bbmin_part_y;
	double bbmin_part_z;
	double bbmax_part_x;
	double bbmax_part_y;
	double bbmax_part_z;

	double bbmin_grid_x;		// bounding box of the grid (loose bb)
	double bbmin_grid_y;
	double bbmin_grid_z;
	double bbmax_grid_x;
	double bbmax_grid_y;
	double bbmax_grid_z;

	double h_max;		// maximum smoothing length
} singleton_geometry;

singleton_geometry * get_singleton_geometry();

#endif /* SINGLETON_GEOMETRY_H_ */
