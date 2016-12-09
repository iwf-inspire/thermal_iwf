#include "ex3_transient_heat.h"

#define M_PI 3.14159265358979323846

bool trans_print = true;

// RK4 temporary storage arrays
static double *k1,*k2,*k3,*k4;

static bool is_edge1 (unsigned int i, unsigned int j, unsigned int k, unsigned int N) {
	return i == 0 || j == 0 || k == 0 || i == N-1 || j == N-1 || k == N-1;
}

static bool is_edge2 (unsigned int i, unsigned int j, unsigned int k, unsigned int N) {
	return i == 1 || j == 1 || k == 1 || i == N-2 || j == N-2 || k == N-2;
}

double solution_transient(double r2, double t, double alpha) {
	return exp(-r2/(4.0*alpha*t))/((4.0*M_PI*alpha*t)*sqrt(4.*M_PI*alpha*t));
}

void transient_print_errors3D(particle* particles, double t, unsigned int Ntot) {
	double L_inf = 0.;
	double L_2 = 0.;
	double EC = 0.;

	for (unsigned int i = 0; i < Ntot; i++) {

		double x  = particles[i].px;
		double y  = particles[i].py;
		double z  = particles[i].pz;
		double r2 = x*x + y*y + z*z;
		double sol = solution_transient(r2,t,particles[i].alpha);
		double err = particles[i].f - sol;

		L_inf = fmax(L_inf,fabs(err));
		L_2  += err*err;
		EC   += fabs(err);
	}
	L_2 /= Ntot;
	EC  /= Ntot;

	printf("Linf=%e; L2=%e; EC=%e; \n",L_inf,sqrt(L_2),EC);
}

static void plot_results(particle* particles, unsigned int Ntot, unsigned int step) {
	char fname[256];
	sprintf(fname,"step_%04d.txt",step);
	FILE *fp = fopen(fname,"w+");

	for (unsigned int i = 0; i < Ntot; i++) {
		fprintf(fp,"%e %e %e %e %e %e\n",particles[i].px, particles[i].py, particles[i].pz, particles[i].f, particles[i].LaplF, particles[i].anal);
	}
	fclose(fp);
}

static void dump(particle* particles, unsigned int N, unsigned int step) {
	char fname[256];
	sprintf(fname,"step_%04d.txt",step);
	FILE *fp = fopen(fname,"w+");

	for (unsigned int i = 0; i < N; i++) {
		for (unsigned int j = 0; j < N; j++) {
			for (unsigned int k = 0; j < N; j++) {
				fprintf(fp,"%e ",particles[i*N*N+j*N+k].f);
			}
			fprintf(fp,"\n");
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
}

void transient_heat_perform_cpu3D(particle* particles, laplacian lap, double dt, unsigned int step, unsigned int Ntot) {
	// as particles are fixed-in-space throughout the TRANSIENT benchmark,
	// we only calculate the correction terms and cell list stuff once.

	lap(particles);

	double totalT = 0.;
	for (unsigned int i = 0; i < Ntot; i++) {
		if(!particles[i].bnd) {
			// update TEMPERATURE-------------------------------------------------------------------------
			particles[i].f += particles[i].alpha * particles[i].LaplF * dt;
		}

		// total concentration
		if (trans_print) totalT += particles[i].f;
	}

	if (trans_print  && step%1==0) {
		//plot_results(particles,Ntot,step);
		dump(particles,21,step);
		printf("sum=%e \n",totalT);
	}
}

void transient_heat_perform_rk4_cpu3D(particle* particles, laplacian lap, double dt, unsigned int step, unsigned int Ntot) {
	double dt2 = 0.5*dt;

	k1 = (double*) calloc(Ntot, sizeof(double));
	k2 = (double*) calloc(Ntot, sizeof(double));
	k3 = (double*) calloc(Ntot, sizeof(double));
	k4 = (double*) calloc(Ntot, sizeof(double));

	// copy the current values before the increment
	for (unsigned int i = 0; i < Ntot; i++) {
		particles[i].T = particles[i].f;
	}

	// start RK4
	lap(particles);
	for (unsigned int i = 0; i < Ntot; i++) {
		if(!particles[i].bnd) {
			k1[i] = particles[i].alpha * particles[i].LaplF;
			particles[i].f = particles[i].T + dt2*k1[i];
		}
	}

	lap(particles);
	for (unsigned int i = 0; i < Ntot; i++) {
		if(!particles[i].bnd) {
			k2[i] = particles[i].alpha * particles[i].LaplF;
			particles[i].f = particles[i].T + dt2*k2[i];
		}
	}

	lap(particles);
	for (unsigned int i = 0; i < Ntot; i++) {
		if(!particles[i].bnd) {
			k3[i] = particles[i].alpha * particles[i].LaplF;
			particles[i].f = particles[i].T + dt*k3[i];
		}
	}

	lap(particles);
	for (unsigned int i = 0; i < Ntot; i++) {
		if(!particles[i].bnd) {
			k4[i] = particles[i].alpha * particles[i].LaplF;
		}
	}

	// update TEMPERATURE-------------------------------------------------------------------------
	double totalT = 0.;
	for (unsigned int i = 0; i < Ntot; i++) {
		if(!particles[i].bnd) {
			particles[i].f = particles[i].T + (1./6.)*dt*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
		}

		// total concentration
		if (trans_print) totalT += particles[i].f;
	}

	if (trans_print) {
		plot_results(particles,Ntot,step);
		printf("sum=%e \n",totalT);
	}


	free(k1);
	free(k2);
	free(k3);
	free(k4);
}

particle* transient_init3D(unsigned int N, double hdx, bool randomized) {
	singleton_geometry* geometry = get_singleton_geometry();
	particle* particles = (particle*) calloc(N*N*N,sizeof(particle));

	double L   = 10.0;
	double dx  = L/(N-1);
	double rho = 1.0;
	double m   = rho*dx*dx*dx;
	double alpha = 1.0;
	double t_init = 0.1;

	for (unsigned int i = 0; i < N; i++) {
		for (unsigned int j = 0; j < N; j++) {
			for (unsigned int k = 0; k < N; k++) {
				double x = i*dx - L/2.;
				double y = j*dx - L/2.;
				double z = k*dx - L/2.;

				unsigned int ID = i*N*N + j*N + k;

				bool boundary = i == 0 || j == 0 || k == 0 || i == N-1 || j == N-1 || k == N-1;

				double ddx = 0.;
				double ddy = 0.;
				double ddz = 0.;
				if (!boundary && randomized) {
					ddx = rand()/((double) RAND_MAX)*0.5-0.25;
					ddy = rand()/((double) RAND_MAX)*0.5-0.25;
					ddz = rand()/((double) RAND_MAX)*0.5-0.25;
				}

				particles[ID].px = x+ddx*dx;
				particles[ID].py = y+ddy*dx;
				particles[ID].pz = z+ddz*dx;
				particles[ID].rho = rho;
				particles[ID].m = m;
				particles[ID].h = hdx*dx;
				particles[ID].alpha = alpha;

				double x2 = particles[ID].px * particles[ID].px;
				double y2 = particles[ID].py * particles[ID].py;
				double z2 = particles[ID].pz * particles[ID].pz;

				particles[ID].T = 0.0;
				particles[ID].f = solution_transient(x2+y2+z2,t_init,particles[ID].alpha);

				particles[ID].idx = i;
				particles[ID].idy = j;
				particles[ID].idz = k;

				particles[ID].bnd   = false;
				particles[ID].blank = false;

				if (is_edge1(i,j,k,N) || is_edge2(i,j,k,N)) {
					particles[ID].bnd = true;
				}

			}
		}
	}

	geometry->n = N*N*N;

	return particles;
}
