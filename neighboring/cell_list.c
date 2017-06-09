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

#include "cell_list.h"

void grid_hash(particle* particles) {
	singleton_geometry *geometry = get_singleton_geometry();
	unsigned int n = geometry->n;

	//spatial hashing of particles
	for (unsigned int i = 0; i < n; i++) {
		unsigned int ix = (unsigned int) ((particles[i].px - geometry->bbmin_part_x)/geometry->dx);
		unsigned int iy = (unsigned int) ((particles[i].py - geometry->bbmin_part_y)/geometry->dx);
		unsigned int iz = (unsigned int) ((particles[i].pz - geometry->bbmin_part_z)/geometry->dx);

		particles[i].hash = ix*(geometry->ny)*(geometry->nz) + iy*(geometry->nz) + iz;
	}
}

void grid_unhash(unsigned int id, int *i, int *j, int *k) {
	singleton_geometry *geometry = get_singleton_geometry();

	*i = id / ((geometry->nz)*(geometry->ny));
	*j = (id - (*i)*(geometry->ny)*(geometry->nz)) / (geometry->nz);
	*k = id % (geometry->nz);

	assert(id == (*i)*(geometry->ny)*(geometry->nz) + (*j)*(geometry->nz) + (*k));
	assert((*i)>=0 && (*i)<(geometry->nx));
	assert((*j)>=0 && (*j)<(geometry->ny));
	assert((*k)>=0 && (*k)<(geometry->nz));
}

void grid_update_grid_geom(particle* particles) {
	singleton_geometry *geometry = get_singleton_geometry();
	unsigned int n = geometry->n;

	double bbmin_part_x = +FLT_MAX; double bbmax_part_x = -FLT_MAX;
	double bbmin_part_y = +FLT_MAX; double bbmax_part_y = -FLT_MAX;
	double bbmin_part_z = +FLT_MAX; double bbmax_part_z = -FLT_MAX;

	double h_max = -FLT_MAX;
	double lx,ly,lz;
	unsigned int nx,ny,nz;
	double dx;

	for (unsigned int i = 0; i < n; i++) {
		h_max = fmax(particles[i].h,h_max);

		bbmin_part_x = fmin(particles[i].px, bbmin_part_x);
		bbmin_part_y = fmin(particles[i].py, bbmin_part_y);
		bbmin_part_z = fmin(particles[i].pz, bbmin_part_z);

		bbmax_part_x = fmax(particles[i].px, bbmax_part_x);
		bbmax_part_y = fmax(particles[i].py, bbmax_part_y);
		bbmax_part_z = fmax(particles[i].pz, bbmax_part_z);
	}

	// prevent round-off error
	bbmin_part_x -= 1e-12;
	bbmin_part_y -= 1e-12;
	bbmin_part_z -= 1e-12;
	bbmax_part_x += 1e-12;
	bbmax_part_y += 1e-12;
	bbmax_part_z += 1e-12;


	lx = (bbmax_part_x - bbmin_part_x);
	ly = (bbmax_part_y - bbmin_part_y);
	lz = (bbmax_part_z - bbmin_part_z);

	dx = 2.f*h_max;         //ToDo: this needs to be dependent on kernel width

	nx = ceil(lx/dx);		//ToDo: this needs to be dependent on kernel width
	ny = ceil(ly/dx);
	nz = ceil(lz/dx);

	if (nx*ny*nz > geometry->num_cell) {
		geometry->cell_start = (int*) realloc (geometry->cell_start,sizeof(int)*nx*ny*nz);
		geometry->cell_end   = (int*) realloc (geometry->cell_end  ,sizeof(int)*nx*ny*nz);
	}

	geometry->dx = dx;

	geometry->nx = nx;
	geometry->ny = ny;
	geometry->nz = nz;

	geometry->num_cell = nx*ny*nz;

	geometry->lx = lx;
	geometry->ly = ly;
	geometry->lz = lz;

	geometry->bbmin_part_x = bbmin_part_x;
	geometry->bbmin_part_y = bbmin_part_y;
	geometry->bbmin_part_z = bbmin_part_z;

	geometry->bbmax_part_x = bbmax_part_x;
	geometry->bbmax_part_y = bbmax_part_y;
	geometry->bbmax_part_z = bbmax_part_z;

	geometry->bbmin_grid_x = bbmin_part_x;
	geometry->bbmin_grid_y = bbmin_part_y;
	geometry->bbmin_grid_z = bbmin_part_z;

	geometry->bbmax_grid_x = bbmin_part_x;
	geometry->bbmax_grid_y = bbmin_part_y;
	geometry->bbmax_grid_z = bbmin_part_z;
}

static int comp_hash(const void *elem1, const void *elem2) {
	particle p1 = *((particle*)elem1);
	particle p2 = *((particle*)elem2);
	if (p1.hash > p2.hash) return  1;
	if (p1.hash < p2.hash) return -1;
	return 0;
}

static int comp_idx(const void *elem1, const void *elem2) {
	particle p1 = *((particle*)elem1);
	particle p2 = *((particle*)elem2);
	if (p1.idx > p2.idx) return  1;
	if (p1.idx < p2.idx) return -1;
	return 0;
}

void grid_restore(particle* particles) {
	singleton_geometry *geometry = get_singleton_geometry();
	unsigned int n = geometry->n;

	qsort(particles,n,sizeof(particle),comp_idx);
}

void grid_sort(particle* particles) {
	singleton_geometry *geometry = get_singleton_geometry();
	unsigned int n = geometry->n;

	// spatial hashing of particles
	grid_hash(particles);

	// sort particles according to spatial hash
	qsort(particles,n,sizeof(particle),comp_hash);

	memset(geometry->cell_start,-1,sizeof(int)*geometry->num_cell);
	memset(geometry->cell_end  ,-1,sizeof(int)*geometry->num_cell);

	int current_box = particles[0].hash;
	geometry->cell_start[current_box] = 0;
	for (unsigned int i = 0; i < n-1; i++) {
		if (particles[i].hash != particles[i+1].hash) {
			geometry->cell_end[current_box] = i+1;
			current_box = particles[i+1].hash;
			assert(current_box >= 0 && current_box < geometry->num_cell);
			geometry->cell_start[current_box] = i+1;
		}
	}

	geometry->cell_end[current_box] = n;
}
