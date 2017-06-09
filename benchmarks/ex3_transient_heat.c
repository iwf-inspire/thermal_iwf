#include "ex3_transient_heat.h"

bool trans_print = true;

// RK4 interim storage arrays
static double *k1,*k2,*k3,*k4;

static bool is_edge1 (unsigned int i, unsigned int j, unsigned int k, unsigned int N) {
	return i == 0 || j == 0 || k == 0 || i == N-1 || j == N-1 || k == N-1;
}

static double solution_transient(double r2, double t, double alpha) {
	return exp(-r2/(4.0*alpha*t))/((4.0*M_PI*alpha*t)*sqrt(4.*M_PI*alpha*t));
}

void transient_print_errors_3D(particle* particles, double t, unsigned int Ntot) {
	double L_inf = 0.;
	double L_2 = 0.;
	double EC = 0.;

	for (unsigned int i = 0; i < Ntot; i++) {
		if(particles[i].bnd) continue;

		double x  = particles[i].px;
		double y  = particles[i].py;
		double z  = particles[i].pz;
		double r2 = x*x + y*y + z*z;
		double sol = solution_transient(r2,t,particles[i].alpha);
		double err = fabs(particles[i].f - sol);

		L_inf = fmax(L_inf,err);
		L_2  += err*err;
		EC   += err;
	}
	L_2  = sqrt(L_2/Ntot);
	EC  /= Ntot;

	printf("Linf=%e; L2=%e; EC=%e; \n",L_inf,L_2,EC);
}

void perform_transient_heat_euler_3D(particle* particles, METHOD method, double dt, unsigned int step, unsigned int Ntot) {

	laplacian_meshfree_method(particles,method);

	double totalT = 0.;
	for (unsigned int i = 0; i < Ntot; i++) {
		if(!particles[i].bnd) {
			// update TEMPERATURE------------------------------------------
			particles[i].f += particles[i].alpha * particles[i].LaplF * dt;
		}
		// total concentration
		totalT += particles[i].f;
	}
	printf("sum=%e \n",totalT);
}

void perform_transient_heat_rk4_3D(particle* particles, METHOD method, double dt, unsigned int step, unsigned int Ntot) {
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
	laplacian_meshfree_method(particles,method);
	for (unsigned int i = 0; i < Ntot; i++) {
		if(!particles[i].bnd) {
			k1[i] = particles[i].alpha * particles[i].LaplF;
			particles[i].f = particles[i].T + dt2*k1[i];
		}
	}

	laplacian_meshfree_method(particles,method);
	for (unsigned int i = 0; i < Ntot; i++) {
		if(!particles[i].bnd) {
			k2[i] = particles[i].alpha * particles[i].LaplF;
			particles[i].f = particles[i].T + dt2*k2[i];
		}
	}

	laplacian_meshfree_method(particles,method);
	for (unsigned int i = 0; i < Ntot; i++) {
		if(!particles[i].bnd) {
			k3[i] = particles[i].alpha * particles[i].LaplF;
			particles[i].f = particles[i].T + dt*k3[i];
		}
	}

	laplacian_meshfree_method(particles,method);
	for (unsigned int i = 0; i < Ntot; i++) {
		if(!particles[i].bnd) {
			k4[i] = particles[i].alpha * particles[i].LaplF;
		}
	}

	// update TEMPERATURE-------------------------------------------------------------------
	double totalT = 0.;
	for (unsigned int i = 0; i < Ntot; i++) {
		if(!particles[i].bnd) {
			particles[i].f = particles[i].T + (1./6.)*dt*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
		}

		// total concentration
		if (trans_print) totalT += particles[i].f;
	}

	printf("sum=%e \n",totalT);

	free(k1);
	free(k2);
	free(k3);
	free(k4);
}

particle* transient_init3D(unsigned int N, double hdx, bool RANDOM) {
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
				unsigned int ID = i*N*N + j*N + k;

				double px = i*dx;
				double py = j*dx;
				double pz = k*dx;

				particles[ID] = make_particle(ID, px, py, pz);

				particles[ID].idx = i;
				particles[ID].idy = j;
				particles[ID].idz = k;

				double x = i*dx - L/2.;
				double y = j*dx - L/2.;
				double z = k*dx - L/2.;

				double ddx = 0.;
				double ddy = 0.;
				double ddz = 0.;
				if (!is_edge1(i,j,k,N) && RANDOM) {
					ddx = rand()/((double) RAND_MAX)*0.5-0.25;
					ddy = rand()/((double) RAND_MAX)*0.5-0.25;
					ddz = rand()/((double) RAND_MAX)*0.5-0.25;
				}

				particles[ID].px  = x + ddx*dx;
				particles[ID].py  = y + ddy*dx;
				particles[ID].pz  = z + ddz*dx;
				particles[ID].rho = rho;
				particles[ID].m = m;
				particles[ID].h = hdx*dx;
				particles[ID].alpha = alpha;

				double x2 = particles[ID].px * particles[ID].px;
				double y2 = particles[ID].py * particles[ID].py;
				double z2 = particles[ID].pz * particles[ID].pz;

				particles[ID].T = 0.0;
				particles[ID].f = solution_transient(x2+y2+z2,t_init,particles[ID].alpha);

				particles[ID].label = NORMAL;
				particles[ID].bnd   = false;

				if (is_edge1(i,j,k,N)) {
					particles[ID].bnd = true;
				}
			}
		}
	}
	geometry->n = N*N*N;

	return particles;
}
