#include "ex1_ex2_reconstruction.h"

particle* setup_reconstruction_random_3D(unsigned int N) {
	singleton_geometry *geometry = get_singleton_geometry();
	particle* particles = (particle*) calloc(N*N*N,sizeof(particle));

	double L = 1.0;
	double dx = L/(N-1);

	for (unsigned int i = 0; i < N; i++) {
		for (unsigned int j = 0; j < N; j++) {
			for (unsigned int k = 0; k < N; k++) {

				double px = rand()/((double) RAND_MAX)*M_PI;
				double py = rand()/((double) RAND_MAX)*M_PI;
				double pz = rand()/((double) RAND_MAX)*M_PI;

				unsigned int ID = i*N*N + j*N + k;

				particle p = make_particle(ID,px,py,pz,N);

				//				p.F    =  sin(p.px)*sin(p.py)*sin(p.pz);
				//				p.anal = -3.*sin(p.px)*sin(p.py)*sin(p.pz);

				p.F    =  (p.px)*(p.px) + (p.py)*(p.py) + (p.pz)*(p.pz) + 2.*(p.px)*(p.py) + 2.*(p.px)*(p.pz) + 2.*(p.py)*(p.pz);
				p.anal = +6.0;


				//				p.F    = (p.px)*(p.px) + (p.py)*(p.py) + (p.pz)*(p.pz);
				//				p.anal = +6.0;


				p.h   = 1.5*dx;
				p.rho = 1.;
				p.m   = dx*dx*dx;

				particles[ID] = p;
			}
		}
	}
	geometry->n = N*N*N;


	return particles;
}

particle* setup_reconstruction_3D(unsigned int N, double hdx) {
	singleton_geometry* geometry = get_singleton_geometry();
	particle* particles = (particle*) calloc(N*N*N,sizeof(particle));

	double L = M_PI;
	double dx = L/(N-1);

	for (unsigned int i = 0; i < N; i++) {
		for (unsigned int j = 0; j < N; j++) {
			for (unsigned int k = 0; k < N; k++) {
				double px = k*dx;
				double py = j*dx;
				double pz = i*dx;

				unsigned int ID = i*N*N + j*N + k;

				particles[ID].idx = k;
				particles[ID].idy = j;
				particles[ID].idz = i;

				particles[ID].px  = px;
				particles[ID].py  = py;
				particles[ID].pz  = pz;

				particles[ID].f    =  sin(px)*sin(py)*sin(pz);
				particles[ID].anal = -3.*sin(px)*sin(py)*sin(pz);
//				particles[ID].f    =  (px)*(px) + (py)*(py) + (pz)*(pz) + 2.*(px)*(py) + 2.*(px)*(pz) + 2.*(py)*(pz);
//				particles[ID].anal = +6.0;

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
