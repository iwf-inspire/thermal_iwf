#include "synthetic_heat_test.h"

static bool rk4 = true;

static void plot_results_heat(particle* particles, unsigned int Ntot, unsigned int step) {
	char fname[256];
	sprintf(fname,"./images/heat_step_%04d.txt",step);
	FILE *fp = fopen(fname,"w+");

	for (unsigned int i = 0; i < Ntot; i++) {
		if(!particles[i].bnd)
			fprintf(fp,"%e %e %e %e\n",particles[i].px, particles[i].py, particles[i].pz, particles[i].f);
	}
	fclose(fp);
}

void run_transient_heat_3D(METHOD method, unsigned int N, bool RANDOM) {
	printf("transient state heat transfer simulation with %d particles \n",N*N*N);

	// parameters *************************************************
	unsigned int Ntot;
	double L, alpha, hdx;
	double dt, t_init, t_final;
	unsigned int step;
	// *************************************************************

	// transient settings
	Ntot = N*N*N;
	L = 10.;
	hdx = 1.5;
	alpha = 1.0;
	dt = 0.25*(L/(N-1))*(L/(N-1))/alpha;
	t_init = 0.1;
	t_final = 0.3;
	step = 1;

	unsigned int Nsteps = fabs(t_final - t_init)/dt;
	assert(Nsteps >= 3);
	assert(Nsteps*dt <= t_final);

	// choose your smoothing kernel
	kernel kernel_func  = &wendland_quintic3D;
	particle* particles = transient_init3D(N,hdx,RANDOM);

	// construct the cell list
	grid_update_grid_geom(particles);
	grid_sort(particles);

	// mesh-free approximation
	initialize_meshfree_method(particles,method,kernel_func,Ntot);
	correction_terms_meshfree_method(particles,method);

	double t = t_init;
	while(t <= t_final) {
		// solve the heat conduction equation
		if(rk4) perform_transient_heat_rk4_3D  (particles,method,dt,step,Ntot);
		else    perform_transient_heat_euler_3D(particles,method,dt,step,Ntot);

		step ++;
		printf("t=%f\n",t);
		t += dt;
	}

	transient_print_errors_3D(particles,t-dt,Ntot);
	grid_restore(particles);
	plot_results_heat(particles,Ntot,step);
	methods_wipe_out();
}

void run_steady_heat_3D(METHOD method, unsigned int N, bool RANDOM) {
	printf("steady state heat transfer simulation with %d particles \n",N*N*N);

	// parameters *************************************************
	unsigned int Nbnd, NN, Ntot;
	double L, alpha, hdx;
	double dt, t_init, t_final;
	unsigned int step;
	// ************************************************************

	// steady settings
	Nbnd = 2;
	NN = N + 2*Nbnd;
	Ntot = NN*NN*NN;
	L = 1.0;
	hdx = 1.5;
	alpha = 1.0;
	dt = 0.25*(L/(N-1))*(L/(N-1))/alpha;
	t_init = 0.0;
	t_final = 1.0;
	step = 1;

	unsigned int Nsteps = fabs(t_final - t_init)/dt;
	assert(Nsteps >= 3);
	assert(Nsteps*dt <= t_final);

	// choose your smoothing kernel
	kernel kernel_func  = &wendland_quintic3D;
	particle* particles = steady_init3D(N,Nbnd,hdx,RANDOM);

	// construct the cell list
	grid_update_grid_geom(particles);
	grid_sort(particles);

	// mesh-free approximation
	initialize_meshfree_method(particles,method,kernel_func,Ntot);
	correction_terms_meshfree_method(particles,method);

	double t = t_init;
	while(t <= t_final) {
		// solve the heat conduction equation
		if(rk4) perform_steady_heat_rk4_3D  (particles,method,dt,step,N,Nbnd);
		else    perform_steady_heat_euler_3D(particles,method,dt,step,N,Nbnd);

		step ++;
		printf("t=%f\n",t);
		t += dt;
	}

	steady_print_errors_3D(particles,N,Nbnd);
	grid_restore(particles);
	plot_results_heat(particles,Ntot,step);
	methods_wipe_out();
}
