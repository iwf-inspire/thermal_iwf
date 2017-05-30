#include "ex1_ex2_reconstruction.h"

particle* setup_reconstruction_polynomial_3D(unsigned int N, double hdx, bool RANDOM) {
	singleton_geometry* geometry = get_singleton_geometry();
	particle* particles = (particle*) calloc(N*N*N,sizeof(particle));

	double L = 1.0;
	double dx = L/(N-1);

	for (unsigned int i = 0; i < N; i++) {
		for (unsigned int j = 0; j < N; j++) {
			for (unsigned int k = 0; k < N; k++) {
				unsigned int ID = i*N*N + j*N + k;

				particles[ID].idx = i;
				particles[ID].idy = j;
				particles[ID].idz = k;

				double px = i*dx;
				double py = j*dx;
				double pz = k*dx;

				bool EDGE = i==0 || j==0 || k==0 || i==N-1 || j==N-1 || k==N-1;

				double ddx = 0.;
				double ddy = 0.;
				double ddz = 0.;
				if (!EDGE && RANDOM) {
					ddx = rand()/((double) RAND_MAX)*0.5-0.25;
					ddy = rand()/((double) RAND_MAX)*0.5-0.25;
					ddz = rand()/((double) RAND_MAX)*0.5-0.25;
				}

				particles[ID].px = px + ddx*dx;
				particles[ID].py = py + ddy*dx;
				particles[ID].pz = pz + ddz*dx;

				particles[ID].f    =  (px)*(px) + (py)*(py) + (pz)*(pz) + 2.*(px)*(py) + 2.*(px)*(pz) + 2.*(py)*(pz);
				particles[ID].anal = +6.0;

				particles[ID].h = hdx*dx;
				particles[ID].rho = 1.;
				particles[ID].m = dx*dx*dx;

				particles[ID].bnd   = false;
				particles[ID].blank = false;
			}
		}
	}
	geometry->n = N*N*N;

	return particles;
}

particle* setup_reconstruction_trigonometric_3D(unsigned int N, double hdx, bool RANDOM) {
	singleton_geometry* geometry = get_singleton_geometry();
	particle* particles = (particle*) calloc(N*N*N,sizeof(particle));

	double L = M_PI;
	double dx = L/(N-1);

	for (unsigned int i = 0; i < N; i++) {
		for (unsigned int j = 0; j < N; j++) {
			for (unsigned int k = 0; k < N; k++) {
				unsigned int ID = i*N*N + j*N + k;

				particles[ID].idx = i;
				particles[ID].idy = j;
				particles[ID].idz = k;

				double px = i*dx;
				double py = j*dx;
				double pz = k*dx;

				bool EDGE = i==0 || j==0 || k==0 || i==N-1 || j==N-1 || k==N-1;

				double ddx = 0.;
				double ddy = 0.;
				double ddz = 0.;
				if (!EDGE && RANDOM) {
					ddx = rand()/((double) RAND_MAX)*0.5-0.25;
					ddy = rand()/((double) RAND_MAX)*0.5-0.25;
					ddz = rand()/((double) RAND_MAX)*0.5-0.25;
				}

				particles[ID].px = px + ddx*dx;
				particles[ID].py = py + ddy*dx;
				particles[ID].pz = pz + ddz*dx;

				particles[ID].f    =  sin(px)*sin(py)*sin(pz);
				particles[ID].anal = -3.*sin(px)*sin(py)*sin(pz);

				particles[ID].h = hdx*dx;
				particles[ID].rho = 1.;
				particles[ID].m = dx*dx*dx;

				particles[ID].bnd   = false;
				particles[ID].blank = false;
			}
		}
	}
	geometry->n = N*N*N;

	return particles;
}
