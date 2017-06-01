#include "ex4_steady_heat.h"

#define M_PI 3.14159265358979323846

static bool steady_print = true;

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

static double solution_steady(double Tamp, double Tbnd, double x, double y, double z) {
	return Tamp*sin(M_PI*x)*sin(M_PI*z)*(sinh(M_PI*y)/sinh(M_PI)) + Tbnd;
}

void steady_print_errors_3D(particle* particles, unsigned int N, unsigned int Nbnd) {
	unsigned int NN = N + 2*Nbnd;
	unsigned int N3 = NN*NN*NN;
	unsigned int Np = N*N*N;
	double L_inf = 0.;
	double L_2 = 0.;
	double EC = 0.;

	for (unsigned int i = 0; i < N3; i++) {
		if(particles[i].bnd) continue;

		double err = fabs(particles[i].f - particles[i].anal);

		L_inf = fmax(L_inf,err);
		L_2  += err*err;
		EC   += err;
	}
	L_2  = sqrt(L_2/Np);
	EC  /= Np;

	printf("Linf=%e; L2=%e; EC=%e; \n",L_inf,L_2,EC);
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

void perform_steady_heat_euler_3D(particle* particles, METHOD method, double dt, unsigned int step, unsigned int N, unsigned int Nbnd) {
	unsigned int NN = N + 2*Nbnd;
	unsigned int Ntot = NN*NN*NN;

	double L = 1.;
	double dx = L/(N-1);
	double Lbnd = (Nbnd)*dx;

	apply_boundary(particles,N,Nbnd,L,Lbnd);
	laplacian_meshfree_method(particles,method);

	double totalT = 0;
	for (unsigned int i = 0; i < Ntot; i++) {
		if(!particles[i].bnd) {
			// update TEMPERATURE-------------------------------------------------------------------------
			particles[i].f += dt * particles[i].alpha * particles[i].LaplF;

			// total concentration
			totalT += particles[i].f;
		}
	}

	printf("step=%d  sum=%e\n",step,totalT);
}

void perform_steady_heat_rk4_3D(particle* particles, METHOD method, double dt, unsigned int step, unsigned int N, unsigned int Nbnd) {
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
	laplacian_meshfree_method(particles,method);
	for (unsigned int i = 0; i < N3; i++) {
		if(particles[i].bnd) continue;
		k1[i] = particles[i].alpha * particles[i].LaplF;
		particles[i].f = particles[i].T + dt2*k1[i];

	}

	laplacian_meshfree_method(particles,method);
	for (unsigned int i = 0; i < N3; i++) {
		if(particles[i].bnd) continue;
		k2[i] = particles[i].alpha * particles[i].LaplF;
		particles[i].f = particles[i].T + dt2*k2[i];

	}

	laplacian_meshfree_method(particles,method);
	for (unsigned int i = 0; i < N3; i++) {
		if(particles[i].bnd) continue;
		k3[i] = particles[i].alpha * particles[i].LaplF;
		particles[i].f = particles[i].T + dt*k3[i];

	}

	laplacian_meshfree_method(particles,method);
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

	printf("step=%d  sum=%e\n",step,totalT);

	free(k1);
	free(k2);
	free(k3);
	free(k4);
}

particle* steady_init3D(unsigned int N, unsigned int Nbnd, double hdx, bool RANDOM) {
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

				particles[ID].idx = i;
				particles[ID].idy = j;
				particles[ID].idz = k;

				double x =(i)*dx - Lbnd;
				double y =(j)*dx - Lbnd;
				double z =(k)*dx - Lbnd;

				double ddx = 0.;
				double ddy = 0.;
				double ddz = 0.;
				if (!is_edge1(i,j,k,NN) && !is_edge2(i,j,k,NN) && is_edge3(i,j,k,NN) && RANDOM) {
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

				particles[ID].label = NORMAL;
				particles[ID].bnd   = false;
				particles[ID].blank = false;

				if (is_edge1(i,j,k,NN) || is_edge2(i,j,k,NN)) particles[ID].bnd = true;
			}
		}
	}
	geometry->n = NN*NN*NN;

	return particles;
}
