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

#ifndef PARTICLE_H_
#define PARTICLE_H_

#include <stdbool.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <string.h>

typedef enum labels {
	// general particles
	NORMAL,
	// ENG example - physical particles
	TOP,
	CONNECTOR_TOP,
	INLET,
	CONNECTOR_INLET,
	CONNECTOR_BAR_CENTER,
	BAR,
	CROSS_1,
	CROSS_2,
	CONNECTOR_HEAT_SINK,
	HEAT_SINK,
	// ENG example - nonphysical particles
	GHOST_BAR,
	GHOST_TOP,
	GHOST_CORE
} labels;

typedef struct particle {
	double rho,m;
	double f,f_x,f_y,f_z;
	double f_xy,f_xz,f_yz;
	double f_xx,f_yy,f_zz;
	double LaplF;
	double anal;
	double F;
	double q;
	unsigned int *neighbors;

	int idx,idy,idz;
	double px,py,pz;
	double vx,vy,vz;
	double derx,dery,derz;
	double derxy;
	double T,LaplT;
	double phi;
	double StrengthBC;
	double h;

	// cell-list stuff
	unsigned int part_id;
	unsigned int hash,nnbh;
	// diffusion coefficient
	double alpha;
	// flag for boundary
	bool bnd;
	// distortion
	double Ht;
	// engineering model
	labels label;

	// flags for surface particles
	bool air_surf;  // bar/top in contact with air
	bool con_surf;  // ball-bearing sides of connector
	bool bar_surf;  // ball-bearing sides of bar
	bool top_inlet; // interface band of connector_top
	bool bot_inlet; // interface band of connector_inl

	bool top_inl_surf; // contact surface particles of CONNECOTR_TOP & INLET
	bool inl_top_surf; // contact surface particles of INLET & CONNECOTR_TOP
	bool bot_inl_surf; // contact surface particles of CONNECOTR_INLET & INLET
	bool inl_bot_surf; // contact surface particles of INLET & CONNECOTR_INLET

} particle;

particle make_particle(unsigned int part_id, double px, double py, double pz);

#endif /* PARTICLE_H_ */
