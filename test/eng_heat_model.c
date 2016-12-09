#define MAX_STEPS 50000

#include "eng_heat_model.h"

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

void run_eng_model3D() {
	singleton_geometry *geometry = get_singleton_geometry();

	// parameters *************************************************
	double dx;
	double hdx = 1.5;
	unsigned int Nbnd = 2;
	double alpha_STL = KA_STL/(RHO_STL*CP_STL);
	double alpha_ALU = KA_ALU/(RHO_ALU*CP_ALU);
	double alpha_PUR = KA_PUR/(RHO_PUR*CP_PUR);
	double alpha_max = fmax(fmax(alpha_STL,alpha_ALU),fmax(alpha_STL,alpha_PUR));

	particle* particles = eng_init3D(hdx,Nbnd,&dx);
	laplacian lap = {&perform_nmfs};
	kernel   kern = {&wendland_quintic3D};

	double dt = 0.5*dx*dx/alpha_max;
	double t_init = 0.0;
	unsigned int step = 0;

	grid_update_grid_geom(particles);
	grid_sort(particles);

	rkpm_init(geometry->n,NULL);
	cspm_init(geometry->n,NULL);
	nmfs_init(geometry->n,NULL);

	/*
	 * as long as particles do NOT move, compute the correction terms only once
	 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	 */
	sph_set_kernel(kern);
	cspm_set_kernel(kern);
	rkpm_set_kernel(kern);
	nmfs_set_kernel(kern);

	//rkpm_compute_correction_terms(particles);
	//cspm_compute_correction_terms(particles);
	nmfs_compute_correction_terms(particles);

	const char *lapl_names = "NMFS";
	printf("CORRECTION TERMS COMPUTED SUCCESSFULLY....SIMULATION WITH %s IS GETTING STARTED\n\n",lapl_names);
	plot_eng(particles,geometry->n,step);

	double totalT_old;
	double totalT = 0;
	double eps = 5e-3;
	double change = 1;

	clock_t start_heat, end_heat;
	start_heat = clock();//===================================
	while(step<MAX_STEPS && change>eps) {
		step++;
		totalT_old = totalT;

		if(rk4) totalT = eng_heat_perform_rk4_cpu3D(particles,lap,2.0*dt,step);
		else    totalT = eng_heat_perform_cpu3D    (particles,lap,dt,step);

		change = totalT - totalT_old;

		if((step)%100==0) {
			printf("total change=%0.3f\n",change);
			plot_eng(particles,geometry->n,step);
		}
	}
	end_heat   = clock();//===================================

	cspm_wipe_out();
	rkpm_wipe_out();
	nmfs_wipe_out();

	// final values
	grid_restore(particles);
	plot_eng(particles,geometry->n,step);

	double time_heat = (end_heat - start_heat)/((double) CLOCKS_PER_SEC);
	printf("CPU took %e seconds\n",time_heat);
};

void run_eng_model3D_dbg() {
	singleton_geometry *geometry = get_singleton_geometry();

	// parameters *************************************************
	double alpha_STL = KA_STL/(RHO_STL*CP_STL);
	double alpha_ALU = KA_ALU/(RHO_ALU*CP_ALU);
	double alpha_PUR = KA_PUR/(RHO_PUR*CP_PUR);
	double alpha_max = fmax(fmax(alpha_STL,alpha_ALU),fmax(alpha_STL,alpha_PUR));
	double hdx = 1.5;
	double dx = 2.84e-03;
	double dt = 0.5*dx*dx/alpha_max;
	double t_init = 0.0;
	unsigned int step = 1;

	particle* particles = eng_init3D_dbg(hdx);
	laplacian lap = {&perform_cspm};
	kernel   kern = {&wendland_quintic3D};

	grid_update_grid_geom(particles);
	grid_sort(particles);

	rkpm_init(geometry->n,NULL);
	cspm_init(geometry->n,NULL);
	nmfs_init(geometry->n,NULL);

	sph_set_kernel(kern);
	cspm_set_kernel(kern);
	rkpm_set_kernel(kern);
	nmfs_set_kernel(kern);

	rkpm_compute_correction_terms(particles);
	cspm_compute_correction_terms(particles);
	nmfs_compute_correction_terms(particles);

	const char *lapl_names = "PSE";
	printf("CORRECTION TERMS COMPUTED SUCCESSFULLY....SIMULATION WITH %s IS GETTING STARTED\n\n",lapl_names);

	double totalT_old;
	double totalT = 0;
	double eps = 5e-2;
	double change = 1;

	clock_t start_heat, end_heat;
	start_heat = clock();//===================================
	while(step<MAX_STEPS && step<20) {

		totalT_old = totalT;
		totalT = eng_heat_perform_cpu3D_dbg(particles,lap,dt,step);
		change = totalT - totalT_old;

		if(step%1 == 0) {
			printf("change=%0.3f step=%d\n",change, step);
			plot_eng(particles,geometry->n,step);
		}
		step++;
	}
	end_heat   = clock();//===================================

	cspm_wipe_out();
	rkpm_wipe_out();
	nmfs_wipe_out();

	// final values
	grid_restore(particles);
	plot_eng(particles,geometry->n,step);

	double time_heat = (end_heat - start_heat)/((double) CLOCKS_PER_SEC);
	printf("CPU took %e seconds\n",time_heat);
};
