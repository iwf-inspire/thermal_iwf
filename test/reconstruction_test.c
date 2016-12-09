/*
 * nbh_bench.c
 *
 *  Created on: Sep 17, 2015
 *      Author: mamzi
 */

#include "reconstruction_test.h"

// global definition of simulation parameters
//-------------------------------------------
#define num_bench   1
#define num_methods 5
#define num_kernels 1
#define num_h       1
//-------------------------------------------

static void plot_results_cpu(particle* particles_cpu, unsigned int Ntot, const char* method) {
	char fname[256];
	sprintf(fname,"%s_%d.txt",method,Ntot);
	FILE *fp = fopen(fname,"w+");

	for (unsigned int i = 0; i < Ntot; i++) {
		fprintf(fp,"%e %e %e %e\n",particles_cpu[i].px, particles_cpu[i].py, particles_cpu[i].pz, particles_cpu[i].LaplF);
	}
	fclose(fp);
}

static void plot_results_anal(particle* particles_cpu, unsigned int Ntot) {
	char fname[256];
	sprintf(fname,"anal.txt");
	FILE *fp = fopen(fname,"w+");

	for (unsigned int i = 0; i < Ntot; i++) {
		fprintf(fp,"%e %e %e %e\n",particles_cpu[i].px, particles_cpu[i].py, particles_cpu[i].pz, particles_cpu[i].anal);
	}
	fclose(fp);
}

static void print_errors_cpu(particle* particles_cpu, unsigned int N,  FILE *fp) {
	double L_inf = 0.;
	double L_2 = 0.;
	double EC = 0.;

	unsigned int Ntot = N*N*N;
	double L = 1.0;
	double dx = L/(N-1);

	//unsigned int iter = 0;
	for (unsigned int i = 0; i < Ntot; i++) {
		double x = particles_cpu[i].px;
		double y = particles_cpu[i].py;
		double z = particles_cpu[i].pz;

//		// error calculation for interior particles only
//		if (x<5.0*dx || x>(L-5.0*dx) || y<5.0*dx || y>(L-5.0*dx) || z<5.0*dx || z>(L-5.0*dx)) continue;

		double err = particles_cpu[i].LaplF - particles_cpu[i].anal;
		L_2  += (err*err);
		EC   += fabs(err);
		L_inf = fmax(L_inf,fabs(err));
//		iter++;
	}
	L_2 = sqrt(L_2/Ntot);
	EC /= Ntot;

	printf("Linf=%e; L2=%e; EC=%e; \n",L_inf,L_2,EC);
	if (fp) {
		fprintf(fp,"%e; %e; %e\n",L_inf,L_2,EC);
	}
}

static void reset_reconstruction3D(particle* particles_cpu, unsigned int Ntot) {
	for (unsigned int i = 0; i < Ntot; i++) {
		particles_cpu[i].LaplF = 0.;
	}
}

static void run_reconstrunction3D(particle* particles_cpu, laplacian lap) {
	lap(particles_cpu);
}

void run_reconstrunction_3D(unsigned int N) {
	printf("3D Laplacian reconstruction with %d particles\n",N*N*N);
	unsigned int Ntot = N*N*N;
	double hdx = 1.5;

	unsigned int part_nums[num_bench] = {10};
	double smoothing_lengths[num_h] = {1.5};
	//const char *lapl_names[num_methods] = {"SPH","PSE","RKPM","CSPM","NMFS"};
	const char *lapl_names[num_methods] = {"RKPM"};
	const char *kernel_names[num_kernels] = {"gaussian"};

	//laplacian laplacians[num_methods] = {&perform_sph, &perform_pse, &perform_rkpm, &perform_cspm, &perform_nmfs};
	laplacian laplacians[num_methods] = {&perform_rkpm};
	kernel kernels[num_kernels] = {&gaussian3D};

	particle* particles_cpu = setup_reconstruction_3D(N,hdx);


	/* -------------------------------------------------------
	 * CPU test
	 * -------------------------------------------------------
	 */
	clock_t start_spatial, end_spatial;
	grid_update_grid_geom(particles_cpu);
	grid_sort(particles_cpu);

	cspm_init(Ntot,NULL);
	rkpm_init(Ntot,NULL);
	nmfs_init(Ntot,NULL);

	sph_set_kernel(kernels[0]);
	cspm_set_kernel(kernels[0]);
	rkpm_set_kernel(kernels[0]);
	nmfs_set_kernel(kernels[0]);

	start_spatial = clock();//===================================
	cspm_compute_correction_terms(particles_cpu);
	rkpm_compute_correction_terms(particles_cpu);
	nmfs_compute_correction_terms(particles_cpu);


	for (unsigned int l = 0; l < num_methods; l++) {
		run_reconstrunction3D(particles_cpu,laplacians[l]);
		//print_errors_cpu(particles_cpu,N);
		reset_reconstruction3D(particles_cpu,Ntot);
	}

	end_spatial = clock();//===================================

	cspm_wipe_out();
	rkpm_wipe_out();
	nmfs_wipe_out();

	//plot_results_cpu(particles_cpu,Ntot);
	//plot_results_anal(particles_cpu,Ntot);
	//print_errors_cpu(particles_cpu,N);

	double time_spatial = (end_spatial - start_spatial)/((double) CLOCKS_PER_SEC);
	printf("CPU spatial took %e seconds\n",time_spatial);
}

void full_perform() {
	FILE *fp = fopen("full_run.txt","w+");

	unsigned int part_nums[num_bench] = {40};
	double smoothing_lengths[num_h] = {1.5};

	const char *lapl_names[num_methods] = {"SPH","PSE","RKPM","CSPM","NMFS"};
	const char *kernel_names[num_kernels] = {"wendland"};

	laplacian laplacians[num_methods] = {&perform_sph, &perform_pse, &perform_rkpm, &perform_cspm, &perform_nmfs};
	kernel kernels[num_kernels] = {&wendland_quintic3D};

	clock_t start_full, end_full;
	start_full = clock();//===================================
	for (unsigned int i = 0; i < num_bench; i++) {
		unsigned int N = part_nums[i];
		unsigned int Ntot = N*N*N;

		for (unsigned int j = 0; j < num_h; j++) {
			double hdx = smoothing_lengths[j];

			particle* particles_cpu = setup_reconstruction_3D(N,hdx);

			// build the cell list stuff ==============================
			grid_update_grid_geom(particles_cpu);
			grid_sort(particles_cpu);

			cspm_init(Ntot,NULL);
			rkpm_init(Ntot,NULL);
			nmfs_init(Ntot,NULL);

			for (unsigned int k = 0; k < num_kernels; k++) {
				char preamble[256];
				sprintf(preamble,"N=%d kernel=%s h=%f",N,kernel_names[k],hdx);
				if(fp) fprintf(fp,"%s\n",preamble);
				printf("N=%d, h=%f kernel=%s\n",N,hdx,kernel_names[k]);

				sph_set_kernel(kernels[k]);
				cspm_set_kernel(kernels[k]);
				rkpm_set_kernel(kernels[k]);
				nmfs_set_kernel(kernels[k]);

				rkpm_compute_correction_terms(particles_cpu);
				cspm_compute_correction_terms(particles_cpu);
				nmfs_compute_correction_terms(particles_cpu);

				printf("correction terms computed!\n");

				for (unsigned int l = 0; l < num_methods; l++) {
					run_reconstrunction3D(particles_cpu,laplacians[l]);
					//plot_results_cpu(particles_cpu,Ntot,lapl_names[l]);
					print_errors_cpu(particles_cpu,N,fp);
					reset_reconstruction3D(particles_cpu,Ntot);
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
