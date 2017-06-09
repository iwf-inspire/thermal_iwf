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

#include "eng_heat_model.h"

#define MAX_STEPS 50000

static void plot_eng(particle* particles, unsigned int Ntot, unsigned int step) {
	char fname[256];
	sprintf(fname,"./images/heat_step_%05d.txt",step);
	FILE *fp = fopen(fname,"w+");

	for (unsigned int i = 0; i < Ntot; i++) {
		if(!particles[i].bnd)
			fprintf(fp,"%e %e %e %e %d\n",particles[i].px, particles[i].py, particles[i].pz, particles[i].f, particles[i].label);
	}
	fclose(fp);
}

void run_eng_model_3D(METHOD method) {
	singleton_geometry *geometry = get_singleton_geometry();
	singleton_physics  *physics  = get_singleton_physics();

	double RHO_STL = physics->RHO_STL;
	double RHO_ALU = physics->RHO_ALU;
	double RHO_PUR = physics->RHO_PUR;
	double KA_STL  = physics->KA_STL;
	double KA_ALU  = physics->KA_ALU;
	double KA_PUR  = physics->KA_PUR;
	double CP_STL  = physics->CP_STL;
	double CP_ALU  = physics->CP_ALU;
	double CP_PUR  = physics->CP_PUR;

	// parameters *****************************
	double dx;
	double hdx = 1.5;
	unsigned int Nbnd = 2;
	double alpha_STL = KA_STL/(RHO_STL*CP_STL);
	double alpha_ALU = KA_ALU/(RHO_ALU*CP_ALU);
	double alpha_PUR = KA_PUR/(RHO_PUR*CP_PUR);
	double alpha_max = fmax(fmax(alpha_STL,alpha_ALU),fmax(alpha_STL,alpha_PUR));

	// choose your smoothing kernel
	kernel kernel_func  = &wendland_quintic3D;
	particle* particles = eng_init3D(hdx,Nbnd,&dx);
	unsigned int Ntot = geometry->n;

	// construct the cell list
	grid_update_grid_geom(particles);
	grid_sort(particles);

	// mesh-free approximation
	initialize_meshfree_method(particles,method,kernel_func,Ntot);
	correction_terms_meshfree_method(particles,method);

	double dt = 0.25*dx*dx/alpha_max;
	unsigned int step = 0;
	double totalT_old = 0;
	double totalT = 0;
	double eps = 5e-3;
	double change = 1;

	while(step<MAX_STEPS && change>eps) {
		totalT_old = totalT;

		if((step)%100==0) {
			printf("total change=%0.3f\n",change);
			plot_eng(particles,Ntot,step);
		}

		step++;

		totalT = perform_eng_heat_euler_3D(particles,method,dt,step);

		change = totalT - totalT_old;
	}

	grid_restore(particles);
	plot_eng(particles,Ntot,step);
	methods_wipe_out();
}
