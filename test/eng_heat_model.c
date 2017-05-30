#include "eng_heat_model.h"

#define MAX_STEPS 50000
static bool rk4 = false;

static void plot_eng(particle* particles, unsigned int Ntot, unsigned int step) {
	char fname[256];
	sprintf(fname,"step_%05d.txt",step);
	FILE *fp = fopen(fname,"w+");
	for (unsigned int i = 0; i < Ntot; i++) {
		if(!particles[i].bnd)
			fprintf(fp,"%e %e %e %e %d\n",particles[i].px, particles[i].py, particles[i].pz, particles[i].f, particles[i].label);
	}
	fclose(fp);
}

void run_eng_model_3D(METHOD method) {
	singleton_geometry *geometry = get_singleton_geometry();

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

	printf("CORRECTION TERMS COMPUTED SUCCESSFULLY....SIMULATION WITH IS GETTING STARTED\n\n");

	double dt = 0.5*dx*dx/alpha_max;
	double t_init = 0.;
	unsigned int step = 0;
	double totalT_old;
	double totalT = 0;
	double eps = 5e-3;
	double change = 1;

	while(step<MAX_STEPS && change>eps) {
		totalT_old = totalT;

		if(rk4) totalT = perform_eng_heat_rk4_3D  (particles,method,dt,step);
		else    totalT = perform_eng_heat_euler_3D(particles,method,dt,step);

		step++;
		change = totalT - totalT_old;

		if((step)%100==0) {
			printf("total change=%0.3f\n",change);
			plot_eng(particles,Ntot,step);
		}
	}

	grid_restore(particles);
	plot_eng(particles,Ntot,step);
	methods_wipe_out();
}
