#include "ex4_steady_heat.h"

#define M_PI 3.14159265358979323846

bool steady_print = true;

// RK4 temporary storage arrays
static double *k1,*k2,*k3,*k4;

static bool is_edge1 (unsigned int i, unsigned int j, unsigned int k, unsigned int N) {
	return i == 0 || j == 0 || k == 0 || i == N-1 || j == N-1 || k == N-1;
}

static bool is_edge2 (unsigned int i, unsigned int j, unsigned int k, unsigned int N) {
	return i == 1 || j == 1 || k == 1 || i == N-2 || j == N-2 || k == N-2;
}

static bool is_edge3 (unsigned int i, unsigned int j, unsigned int k, unsigned int N) {
	return i == 2 || j == 2 || k == 2 || i == N-3 || j == N-3 || k == N-3;
}

double solution_steady(double Tamp, double Tbnd, double x, double y, double z) {
	return Tamp*sin(M_PI*x)*sin(M_PI*z)*(sinh(M_PI*y)/sinh(M_PI)) + Tbnd;
}

void steady_print_errors3D(particle* particles, unsigned int N, unsigned int Nbnd) {
	unsigned int NN = N + 2*Nbnd;
	unsigned int N3 = NN*NN*NN;
	unsigned int Np = N*N*N;
	double L_inf = 0.;
	double L_2 = 0.;
	double EC = 0.;

	for (unsigned int i = 0; i < N3; i++) {

		double x = particles[i].px;
		double y = particles[i].py;
		double z = particles[i].pz;

		if (x>=0 && x<=1 && y>=0 && y<=1 && z>=0 && z<=1) {

			double err = particles[i].f - particles[i].anal;

			L_2   += err*err;
			L_inf  = fmax(L_inf,fabs(err));
			EC    += fabs(err);
		}
	}
	L_2 /= Np;
	EC  /= Np;

	printf("Linf=%e; L2=%e; EC=%e; \n",L_inf,sqrt(L_2),EC);
}

static void apply_boundary(particle* particles, unsigned int N, unsigned int Nbnd, double L, double Lbnd) {
	unsigned int NN = N + 2*Nbnd;
	unsigned int N3 = NN*NN*NN;

	for (unsigned int i = 0; i < N3; i++) {
		if (!particles[i].bnd) continue;

		particles[i].f = 0.;

		double xp = particles[i].px;
		double yp = particles[i].py;
		double zp = particles[i].pz;

		if (xp>=0 && xp<=L && zp>=0 && zp<=L && yp>=L)
			particles[i].f = sin(M_PI*xp)*sin(M_PI*zp);
	}
}

static void plot_results(particle* particles, unsigned int N, unsigned int Nbnd, unsigned int step) {
	unsigned int NN = N + 2*Nbnd;
	unsigned int N3 = NN*NN*NN;

	char fname[256];
	sprintf(fname,"step_%04d.txt",step);
	FILE *fp = fopen(fname,"w+");

	for (unsigned int i = 0; i < N3; i++) {
		double xp = particles[i].px;
		double yp = particles[i].py;
		double zp = particles[i].pz;

		fprintf(fp,"%e %e %e %e %e\n",particles[i].px, particles[i].py, particles[i].pz, particles[i].f, particles[i].anal);
	}
	fclose(fp);
}

void plot_test_results(particle* particles, unsigned int N, unsigned int step) {
	unsigned int N3 = N*N*N;

	char fname[256];
	sprintf(fname,"step_%04d.txt",step);
	FILE *fp = fopen(fname,"w+");

	for (unsigned int i = 0; i < N3; i++) {
		double xp = particles[i].px;
		double yp = particles[i].py;
		double zp = particles[i].pz;

		fprintf(fp,"%e %e %e %e\n",particles[i].px, particles[i].py, particles[i].pz, particles[i].f);
	}
	fclose(fp);
}

void steady_heat_perform_cpu3D(particle* particles, laplacian lap, double dt, unsigned int step, unsigned int N, unsigned int Nbnd) {
	unsigned int NN = N + 2*Nbnd;
	unsigned int N3 = NN*NN*NN;

	double L = 1.;
	double dx = L/(N-1);
	double Lbnd = (Nbnd)*dx;

	// particles are fixed-in-space:
	// calculate the correction terms and cell list stuff only once.

	apply_boundary(particles,N,Nbnd,L,Lbnd);
	lap(particles);

	double totalT = 0;
	for (unsigned int i = 0; i < N3; i++) {
		if(!particles[i].bnd) {
			// update TEMPERATURE-------------------------------------------------------------------------
			particles[i].f += dt * particles[i].alpha * particles[i].LaplF;

			// total concentration
			totalT += particles[i].f;
		}
	}

	if (steady_print && step%40==0) {
		//plot_results(particles,N,Nbnd,step);
		printf("step=%d  sum=%e\n",step,totalT);
	}
}

void steady_heat_perform_test_cpu3D(particle* particles, laplacian lap, double dt, unsigned int step, unsigned int N) {
	unsigned int N3 = N*N*N;

	lap(particles);

	double totalT = 0;
	for (unsigned int i = 0; i < N3; i++) {
		if(!particles[i].bnd) {
			// update TEMPERATURE-------------------------------------------------------------------------
			particles[i].f += dt * particles[i].alpha * particles[i].LaplF;
		}
		// total concentration
		totalT += particles[i].f;
	}

	if (steady_print && step%10==0) {
		plot_test_results(particles,N,step);
		printf("step=%d  sum=%e\n",step,totalT);
	}
}

void steady_heat_perform_rk4_cpu3D(particle* particles, laplacian lap, double dt, unsigned int step, unsigned int N, unsigned int Nbnd) {
	double dt2 = 0.5*dt;

	unsigned int NN = N + 2*Nbnd;
	unsigned int N3 = NN*NN*NN;

	double L = 1.;
	double dx = L/(N-1);
	double Lbnd = (Nbnd)*dx;

	k1 = (double*) calloc(N3, sizeof(double));
	k2 = (double*) calloc(N3, sizeof(double));
	k3 = (double*) calloc(N3, sizeof(double));
	k4 = (double*) calloc(N3, sizeof(double));

	apply_boundary(particles,N,Nbnd,L,Lbnd);

	// copy the current values before the increment
	for (unsigned int i = 0; i < N3; i++) {
		if(particles[i].bnd) continue;
		particles[i].T = particles[i].f;
	}

	// start RK4
	lap(particles);
	for (unsigned int i = 0; i < N3; i++) {
		if(particles[i].bnd) continue;
		k1[i] = particles[i].alpha * particles[i].LaplF;
		particles[i].f = particles[i].T + dt2*k1[i];

	}

	lap(particles);
	for (unsigned int i = 0; i < N3; i++) {
		if(particles[i].bnd) continue;
		k2[i] = particles[i].alpha * particles[i].LaplF;
		particles[i].f = particles[i].T + dt2*k2[i];

	}

	lap(particles);
	for (unsigned int i = 0; i < N3; i++) {
		if(particles[i].bnd) continue;
		k3[i] = particles[i].alpha * particles[i].LaplF;
		particles[i].f = particles[i].T + dt*k3[i];

	}

	lap(particles);
	for (unsigned int i = 0; i < N3; i++) {
		if(particles[i].bnd) continue;
		k4[i] = particles[i].alpha * particles[i].LaplF;

	}


	// update TEMPERATURE-------------------------------------------------------------------------
	double totalT = 0.;
	for (unsigned int i = 0; i < N3; i++) {
		if(particles[i].bnd) continue;
		particles[i].f = particles[i].T + (1./6.)*dt*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);

		// total concentration
		if (steady_print) totalT += particles[i].f;

	}

	if (steady_print && step%4==0) {
		plot_results(particles,N,Nbnd,step);
		printf("step=%d  sum=%e\n",step,totalT);
	}

	free(k1);
	free(k2);
	free(k3);
	free(k4);
}

particle* steady_init3D(unsigned int N, unsigned int Nbnd, double hdx, bool randomized) {
	unsigned int NN = N+2*Nbnd;
	unsigned int N3 = NN*NN*NN;

	singleton_geometry *geometry = get_singleton_geometry();
	particle* particles = (particle*) calloc(N3,sizeof(particle));

	double L = 1.;
	double dx = L/(N-1);
	double Lbnd = (Nbnd)*dx;
	double rho = 1.0;
	double m   = rho*dx*dx*dx;
	double Tamp = 1.;
	double Tbnd = 0.;
	double alpha = 1.;

	for (unsigned int i = 0; i < NN; i++) {
		for (unsigned int j = 0; j < NN; j++) {
			for (unsigned int k = 0; k < NN; k++) {
				unsigned int ID = i*NN*NN + j*NN + k;

				double x =(i)*dx - Lbnd;
				double y =(j)*dx - Lbnd;
				double z =(k)*dx - Lbnd;

				double ddx = 0.;
				double ddy = 0.;
				double ddz = 0.;
				if (!is_edge1(i,j,k,NN) && !is_edge2(i,j,k,NN) && is_edge3(i,j,k,NN) && randomized) {
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
				particles[ID].anal = solution_steady(Tamp,Tbnd,particles[ID].px,particles[ID].py,particles[ID].pz);
				particles[ID].f = 0.0;
				particles[ID].T = 0.0;

				particles[ID].idx = i;
				particles[ID].idy = j;
				particles[ID].idz = k;

				particles[ID].bnd   = false;
				particles[ID].blank = false;

				if (is_edge1(i,j,k,NN) || is_edge2(i,j,k,NN)) particles[ID].bnd = true;

			}
		}
	}
	geometry->n = NN*NN*NN;


	return particles;
}

particle* steady_test_init3D(unsigned int N, double hdx, bool randomized) {
	unsigned int N3 = N*N*N;

	singleton_geometry *geometry = get_singleton_geometry();
	particle* particles = (particle*) calloc(N3,sizeof(particle));

	double L = 1.;
	double dx = L/(N-1);
	double rho = 1.0;
	double m   = rho*dx*dx*dx;
	double alpha = 1.;

	for (unsigned int i = 0; i < N; i++) {
		for (unsigned int j = 0; j < N; j++) {
			for (unsigned int k = 0; k < N; k++) {
				unsigned int ID = i*N*N + j*N + k;

				bool rear  = i==0 || i==1;
				bool front = i==N-1 || i==N-2;
				bool edge  = i==0 || j==0 || k==0 || i==N-1 || j==N-1 || k==N-1;

				double x =(i)*dx;
				double y =(j)*dx;
				double z =(k)*dx;

				double ddx = 0.;
				double ddy = 0.;
				double ddz = 0.;
				if (!is_edge1(i,j,k,N) && !is_edge2(i,j,k,N) && randomized) {
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
				particles[ID].T = 0.0;
				particles[ID].f = (rear) ? 1.0 : 0.0;

				particles[ID].idx = i;
				particles[ID].idy = j;
				particles[ID].idz = k;

				particles[ID].bnd   = false;
				particles[ID].blank = false;

			    if(edge) particles[ID].bnd = true;
			}
		}
	}
	geometry->n = N3;

	return particles;
}
