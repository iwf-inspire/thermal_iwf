#include "synthetic_heat_test.h"

// global definition of simulation parameters
//-------------------------------------------
#define num_bench   2
#define num_methods 4
#define num_kernels 1
#define num_h       1
//-------------------------------------------

static bool rk4 = false;

static void plot_results(particle* particles_cpu, unsigned int Ntot, unsigned int step) {
	char fname[256];
	sprintf(fname,"step_%04d.txt",step);
	FILE *fp = fopen(fname,"w+");

	for (unsigned int i = 0; i < Ntot; i++) {
		if(!particles_cpu[i].bnd)
			fprintf(fp,"%e %e %e %e\n",particles_cpu[i].px, particles_cpu[i].py, particles_cpu[i].pz, particles_cpu[i].f);
	}
	fclose(fp);
}

static void plot_results_anal(particle* particles, unsigned int Ntot, unsigned int step, double t) {
	char fname[256];
	sprintf(fname,"anal_%04d.txt",step);
	FILE *fp = fopen(fname,"w+");

	for (unsigned int i = 0; i < Ntot; i++) {
		double x  = particles[i].px;
		double y  = particles[i].py;
		double z  = particles[i].pz;
		double r2 = x*x + y*y + z*z;
		double sol = solution_transient(r2,t,particles[i].alpha);
		fprintf(fp,"%e %e %e %e\n",particles[i].px, particles[i].py, particles[i].pz, sol);
	}
	fclose(fp);
}

static void print_thermal_errors_cpu_transient(particle* particles_cpu, double t, unsigned int N,  FILE *fp) {
	double L_inf = 0.;
	double L_2 = 0.;
	double EC = 0.;
	unsigned int Ntot = N*N*N;

	for (unsigned int i = 0; i < Ntot; i++) {
		double x = particles_cpu[i].px;
		double y = particles_cpu[i].py;
		double z = particles_cpu[i].pz;

		double r2 = x*x + y*y + z*z;
		double sol = solution_transient(r2,t,particles_cpu[i].alpha);
		double err = particles_cpu[i].f - sol;

		L_2 += err*err;
		L_inf = fmax(L_inf,fabs(err));
		EC += fabs(err);
	}
	L_2 /= Ntot;
	EC  /= Ntot;

	printf("Linf=%e; L2=%e; EC=%e; \n",L_inf,sqrt(L_2),EC);
	if (fp) {
		fprintf(fp,"%e; %e; %e\n",L_inf,sqrt(L_2),EC);
	}
}

static void print_thermal_errors_cpu_steady(particle* particles_cpu, unsigned int N, unsigned int Nbnd, FILE *fp) {
	unsigned int NN = N + 2*Nbnd;
	unsigned int N3 = NN*NN*NN;
	unsigned int Np = N*N*N;
	double L_inf = 0.;
	double L_2 = 0.;
	double EC = 0.;

	for (unsigned int i = 0; i < N3; i++) {

		if (particles_cpu[i].bnd) continue;

		double x = particles_cpu[i].px;
		double y = particles_cpu[i].py;
		double z = particles_cpu[i].pz;

		double err = particles_cpu[i].f - particles_cpu[i].anal;

		L_2   += err*err;
		L_inf  = fmax(L_inf,fabs(err));
		EC    += fabs(err);
	}
	L_2 /= Np;
	EC  /= Np;

	printf("Linf=%e; L2=%e; EC=%e; \n",L_inf,sqrt(L_2),EC);
	if (fp) {
		fprintf(fp,"%e; %e; %e\n",L_inf,sqrt(L_2),EC);
	}
}

static void reset_thermal3D_transient(particle* particles_cpu, double t_init, double alpha, unsigned int Ntot) {
	for (unsigned int i = 0; i < Ntot; i++) {
		particles_cpu[i].LaplF = 0.;
		double x = particles_cpu[i].px;
		double y = particles_cpu[i].py;
		double z = particles_cpu[i].pz;
		particles_cpu[i].f = solution_transient(x*x+y*y+z*z,t_init,alpha);
	}
}

static double run_thermal3D_transient(particle* particles_cpu, laplacian lap, double t, double t_final, double dt, unsigned int step, unsigned int Ntot) {
	while(t <= t_final) {
		if(rk4) transient_heat_perform_rk4_cpu3D(particles_cpu,lap,dt,step,Ntot);
		else    transient_heat_perform_cpu3D(particles_cpu,lap,dt,step,Ntot);
		step ++;
		t += dt;
	}
	return (t-dt);
}

static void reset_thermal3D_steady(particle* particles_cpu, unsigned int Ntot) {
	for (unsigned int i = 0; i < Ntot; i++) {
		particles_cpu[i].LaplF = 0.;
		particles_cpu[i].f     = 0.;
	}
}

static void run_thermal3D_steady(particle* particles_cpu, laplacian lap, double t, double t_final, double dt, unsigned int step, unsigned int N, unsigned int Nbnd) {
	while(t <= t_final) {
		if(rk4) steady_heat_perform_rk4_cpu3D(particles_cpu,lap,dt,step,N,Nbnd);
		else 	steady_heat_perform_cpu3D(particles_cpu,lap,dt,step,N,Nbnd);
		// if(step%40 == 0) steady_print_errors3D(particles_cpu,N,Nbnd);
		step ++;
		t += dt;
	}
}

/*
void full_perform_thermal3D() {
	FILE *fp = fopen("full_run_heat3D_TRNAS.txt","w+");

	// parameters *************************************************
	unsigned int N, Nbnd, NN, Ntot;
	double L, alpha;
	double Tamp, Tbnd;
	double dt, t_init, t_final;
	unsigned int step;
	// *************************************************************

	unsigned int part_nums[num_bench] = {41,51};
	double smoothing_lengths[num_h] = {1.5};
	const char *lapl_names[num_methods] = {"PSE","RKPM","CSPM","NMFS"};
	const char *kernel_names[num_kernels] = {"wendland quintic"};

	laplacian laplacians[num_methods] = {&perform_pse_cpu,&perform_rkpm_cpu,&perform_cspm_cpu,&perform_nmfs_cpu};
	kernel kernels[num_kernels] = {&wendland_quintic3D};

	clock_t start_full, end_full;
	start_full = clock();//===================================
	for (unsigned int i = 0; i < num_bench; i++) {
		unsigned int N = part_nums[i];

		if (steady_bench) {
			// steady
			Nbnd = 2;
			NN = N + 2*Nbnd;
			Ntot = NN*NN*NN;
			L = 1.0;
			alpha = 1.0;
			Tamp = 1.0;
			Tbnd = 0.0;
			dt = 0.5*(L/(N-1))*(L/(N-1))/alpha;
			t_init = 0.;
			t_final = 0.6;
			step = 1;

		} else {
			// transient
			Ntot = N*N*N;
			L = 10.;
			alpha = 1.0;
			dt = 0.25*(L/(N-1))*(L/(N-1))/alpha;
			t_init = 0.1;
			t_final = 0.3;
			step = 1;
		}


		unsigned int Nsteps = fabs(t_final - t_init)/dt;
		assert(Nsteps >= 3);
		assert(Nsteps*dt <= t_final);

		for (unsigned int j = 0; j < num_h; j++) {
			double hdx = smoothing_lengths[j];

			particle* particles_cpu = (steady_bench) ? steady_init3D(N,Nbnd,hdx,false) : transient_init3D(N,hdx,false);

			// build the cell list stuff ==============================
			grid_update_grid_geom(particles_cpu);
			grid_sort(particles_cpu);

			rkpm_init(Ntot,NULL);
			cspm_init(Ntot,NULL);
			nmfs_init(Ntot,NULL);

			for (unsigned int k = 0; k < num_kernels; k++) {
				char preamble[256];
				sprintf(preamble,"N=%d kernel=%s h=%f",N,kernel_names[k],hdx);
				if(fp) fprintf(fp,"%s\n",preamble);
				printf("N=%d, Nsteps=%d, h=%f kernel=%s\n",N,Nsteps,hdx,kernel_names[k]);

				sph_set_kernel(kernels[k]);
				cspm_set_kernel(kernels[k]);
				rkpm_set_kernel(kernels[k]);
				nmfs_set_kernel(kernels[k]);

				rkpm_compute_correction_terms_cpu(particles_cpu);
				cspm_compute_correction_terms_cpu(particles_cpu);
				nmfs_compute_correction_terms_cpu(particles_cpu);

				printf("correction terms computed!\n");

				for (unsigned int l = 0; l < num_methods; l++) {

					if (steady_bench) {
						run_thermal3D_steady(particles_cpu,laplacians[l],t_init,t_final,dt,step,N,Nbnd);
						print_thermal_errors_cpu_steady(particles_cpu,N,Nbnd,fp);
						reset_thermal3D_steady(particles_cpu,Ntot);
					} else {
						double time = run_thermal3D_transient(particles_cpu,laplacians[l],t_init,t_final,dt,step,Ntot);
						print_thermal_errors_cpu_transient(particles_cpu,time,N,fp);
						reset_thermal3D_transient(particles_cpu,t_init,1.0,Ntot);
					}

				}
			}
			rkpm_wipe_out();
			cspm_wipe_out();
			nmfs_wipe_out();
		}
	}
	end_full = clock();//===================================

	double time_full = (end_full - start_full)/((double) CLOCKS_PER_SEC);
	printf("CPU reconstruction took %f seconds\n",time_full);
}

 */

void transient_heat_test(unsigned int N) {
	printf("heat transfer simulation with %d particles in CPU \n",N*N*N);

	// parameters *************************************************
	unsigned int Ntot;
	double L, alpha;
	double dt, t_init, t_final;
	unsigned int step;
	// *************************************************************

	// transient settings
	Ntot = N*N*N;
	L = 10.;
	alpha = 1.0;
	dt = 0.2*(L/(N-1))*(L/(N-1))/alpha;
	t_init = 0.1;
	t_final = 0.3;
	step = 1;

	unsigned int Nsteps = fabs(t_final - t_init)/dt;
	assert(Nsteps >= 3);
	assert(Nsteps*dt <= t_final);

	particle* p   = transient_init3D(N,1.5,false);
	laplacian lap = {&perform_nmfs};
	kernel   kern = {&wendland_quintic3D};

	grid_update_grid_geom(p);
	grid_sort(p);

	rkpm_init(Ntot,NULL);
	cspm_init(Ntot,NULL);
	nmfs_init(Ntot,NULL);

	/*
	 * as long as particles do NOT move, compute the correction terms only once
	 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	 */
	sph_set_kernel(kern);
	cspm_set_kernel(kern);
	rkpm_set_kernel(kern);
	nmfs_set_kernel(kern);

	rkpm_compute_correction_terms(p);
	cspm_compute_correction_terms(p);
	nmfs_compute_correction_terms(p);

	const char *lapl_names = "NMFS";
	printf("SIMULATION WITH %s\n",lapl_names);

	double t = t_init;
	clock_t start_heat, end_heat;
	start_heat = clock();//===================================
	while(t <= t_final) {

		if(rk4) transient_heat_perform_rk4_cpu3D(p,lap,dt,step,Ntot);
		else    transient_heat_perform_cpu3D(p,lap,dt,step,Ntot);

		step ++;
		printf("t=%f\n",t);
		t += dt;
	}
	end_heat   = clock();//===================================

	transient_print_errors3D(p,t-dt,Ntot);

	cspm_wipe_out();
	rkpm_wipe_out();
	nmfs_wipe_out();

	grid_restore(p);

	double time_heat = (end_heat - start_heat)/((double) CLOCKS_PER_SEC);
	printf("CPU took %e seconds\n",time_heat);
}

void steady_heat_test(unsigned int N) {
	printf("heat transfer simulation with %d particles in CPU \n",N*N*N);

	// parameters *************************************************
	unsigned int Nbnd, NN, Ntot;
	double L, alpha;
	double Tamp, Tbnd;
	double dt, t_init, t_final;
	unsigned int step;
	// *************************************************************

	// steady settings
	Nbnd = 2;
	NN = N + 2*Nbnd;
	Ntot = NN*NN*NN;
	L = 1.0;
	alpha = 1.0;
	Tamp = 1.0;
	Tbnd = 0.0;
	dt = 0.5*(L/(N-1))*(L/(N-1))/alpha;
	t_init = 0.0;
	t_final = 1.0;
	step = 1;

	unsigned int Nsteps = fabs(t_final - t_init)/dt;
	assert(Nsteps >= 3);
	assert(Nsteps*dt <= t_final);

	particle* p   = steady_init3D(N,Nbnd,1.5,false);
	laplacian lap = {&perform_nmfs};
	kernel   kern = {&wendland_quintic3D};

	grid_update_grid_geom(p);
	grid_sort(p);

	rkpm_init(Ntot,NULL);
	cspm_init(Ntot,NULL);
	nmfs_init(Ntot,NULL);

	/*
	 * as long as particles do NOT move, compute the correction terms only once
	 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	 */
	sph_set_kernel(kern);
	cspm_set_kernel(kern);
	rkpm_set_kernel(kern);
	nmfs_set_kernel(kern);

	rkpm_compute_correction_terms(p);
	cspm_compute_correction_terms(p);
	nmfs_compute_correction_terms(p);

	const char *lapl_names = "NMFS";
	printf("SIMULATION WITH %s\n",lapl_names);

	double t = t_init;
	clock_t start_heat, end_heat;
	start_heat = clock();//===================================
	while(t <= t_final) {

		if(rk4) steady_heat_perform_rk4_cpu3D(p,lap,dt,step,N,Nbnd);
		else    steady_heat_perform_cpu3D(p,lap,dt,step,N,Nbnd);

		step ++;
		printf("t=%f\n",t);
		t += dt;
	}
	end_heat   = clock();//===================================

	steady_print_errors3D(p,N,Nbnd);

	cspm_wipe_out();
	rkpm_wipe_out();
	nmfs_wipe_out();

	grid_restore(p);

	double time_heat = (end_heat - start_heat)/((double) CLOCKS_PER_SEC);
	printf("CPU took %e seconds\n",time_heat);
}
