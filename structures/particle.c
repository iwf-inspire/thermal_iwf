#include "particle.h"

particle make_particle(unsigned int idx, double px, double py, double pz, unsigned int N) {
	particle p;
	memset(&p, 0, sizeof(double)*N*N*N);

	p.idx = idx;
	p.px = px;
	p.py = py;
	p.pz = pz;

	return p;
}
