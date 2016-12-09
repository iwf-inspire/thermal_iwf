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
	// physical particles
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
	// non-physical particles
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

	unsigned int part_id;
	unsigned int hash,nnbh;	 // cell hash

	double rho0;		// initial density (needed if initially particles are not spaced regularly)
	double rho_a;

	double vx_x;		// velocity gradient
	double vx_y;
	double vx_z;
	double vy_x;
	double vy_y;
	double vy_z;
	double vz_x;
	double vz_y;
	double vz_z;

	double alpha;       // diffusion coefficient

	bool bnd,blank;		// flags for boundary and blanking
	// bnd: means that the particle does not act as destination particle
	// blank: means that the particle does neither act as a source nor destination term
	// this means the particle is effectively removed from the simulation

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

particle make_particle(unsigned int idx, double px, double py, double pz, unsigned int N);

#endif /* PARTICLE_H_ */
