#include "reconstruction_test.h"

static void plot_results_reconstruction(particle* particles, unsigned int Ntot) {
	char fname[256];
	sprintf(fname,"./images/reconstructiopn_results.txt");
	FILE *fp = fopen(fname,"w+");

	for (unsigned int i = 0; i < Ntot; i++) {
		fprintf(fp,"%e %e %e %e\n",particles[i].px, particles[i].py, particles[i].pz, particles[i].LaplF);
	}
	fclose(fp);
}

static void print_errors_reconstruction(particle* particles, unsigned int Ntot) {
	double L_inf = 0.;
	double L_2 = 0.;
	double EC = 0.;

	for (unsigned int i = 0; i < Ntot; i++) {
		if (particles[i].bnd) continue;

		double err = (particles[i].LaplF - particles[i].anal);
		L_inf = fmax(L_inf,fabs(err));
		L_2  += err*err;
		EC   += fabs(err);
	}
	L_2 = sqrt(L_2/Ntot);
	EC /= Ntot;

	printf("Linf=%e; L2=%e; EC=%e; \n",L_inf,L_2,EC);
}

void run_reconstrunction_polynomial_3D(METHOD method, unsigned int N, bool RANDOM) {
	printf("polynomial reconstruction benchmark with %d particles \n",N*N*N);

	unsigned int Ntot = N*N*N;
	double hdx = 1.5;

	// choose your smoothing kernel
	kernel kernel_func  = &gaussian3D;
	particle* particles = setup_reconstruction_polynomial_3D(N,hdx,RANDOM);

	// construct the cell list
	grid_update_grid_geom(particles);
	grid_sort(particles);

	// mesh-free approximation
	initialize_meshfree_method(particles,method,kernel_func,Ntot);
	correction_terms_meshfree_method(particles,method);
	laplacian_meshfree_method(particles,method);
	methods_wipe_out();
	grid_restore(particles);

	// output
	plot_results_reconstruction(particles,Ntot);
	print_errors_reconstruction(particles,Ntot);
}

void run_reconstrunction_trigonometric_3D(METHOD method, unsigned int N, bool RANDOM) {
	printf("trigonometric reconstruction benchmark with %d particles \n",N*N*N);

	unsigned int Ntot = N*N*N;
	double hdx = 1.5;

	// choose your smoothing kernel
	kernel kernel_func  = &cubic_spline3D;
	particle* particles = setup_reconstruction_trigonometric_3D(N,hdx,RANDOM);

	// construct the cell list
	grid_update_grid_geom(particles);
	grid_sort(particles);

	// mesh-free approximation
	initialize_meshfree_method(particles,method,kernel_func,Ntot);
	correction_terms_meshfree_method(particles,method);
	laplacian_meshfree_method(particles,method);
	methods_wipe_out();
	grid_restore(particles);

	// output
	plot_results_reconstruction(particles,Ntot);
	print_errors_reconstruction(particles,Ntot);
}
