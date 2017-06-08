#include "sph.h"

static kernel _kernel;

void perform_sph(particle* particles) {
	singleton_geometry *grid_geom  = get_singleton_geometry();

	for (unsigned int b = 0; b < grid_geom->num_cell; b++) {

		// finding neighbors
		int gi = 0;
		int gj = 0;
		int gk = 0;
		grid_unhash(b,&gi,&gj,&gk);

		int low_i  = gi-1 < 0 ? 0 : gi-1;
		int low_j  = gj-1 < 0 ? 0 : gj-1;
		int low_k  = gk-1 < 0 ? 0 : gk-1;
		int high_i = gi+2 > grid_geom->nx ? grid_geom->nx : gi+2;
		int high_j = gj+2 > grid_geom->ny ? grid_geom->ny : gj+2;
		int high_k = gk+2 > grid_geom->nz ? grid_geom->nz : gk+2;

		for (unsigned int i = grid_geom->cell_start[b]; i < grid_geom->cell_end[b]; i++) {

			if (particles[i].bnd) continue;

			particle part 	= particles[i];
			double xi 		= part.px;
			double yi 		= part.py;
			double zi 		= part.pz;
			double hi 		= part.h;
			double fi       = part.f;

			double f_x  = 0.; double f_y  = 0.; double f_z  = 0.;
			double f_xx = 0.; double f_yy = 0.; double f_zz = 0.;

			// loop over neighbors
			for (unsigned int ni = low_i; ni < high_i; ni++) {
				for (unsigned int nj = low_j; nj < high_j; nj++) {
					for (unsigned int nk = low_k; nk < high_k; nk++) {

						unsigned int ID = ni*grid_geom->ny*grid_geom->nz + nj*grid_geom->nz + nk;

						for (unsigned int j = grid_geom->cell_start[ID]; j < grid_geom->cell_end[ID]; j++) {

							particle part 	= particles[j];
							double xj 		= part.px;
							double yj 		= part.py;
							double zj 		= part.pz;
							double rhoj     = part.rho;
							double m        = part.m;
							double fj       = part.f;

							kernel_result w = _kernel(xi,yi,zi,xj,yj,zj,hi);

							f_x  += (fj - fi)*w.w_x* m/rhoj;
							f_y  += (fj - fi)*w.w_y* m/rhoj;
							f_z  += (fj - fi)*w.w_z* m/rhoj;
							f_xx += (fj - fi)*w.w_xx*m/rhoj;
							f_yy += (fj - fi)*w.w_yy*m/rhoj;
							f_yy += (fj - fi)*w.w_zz*m/rhoj;
						}
					}
				}
			}
			particles[i].f_x   = f_x;
			particles[i].f_y   = f_y;
			particles[i].f_z   = f_z;
			particles[i].LaplF = f_xx + f_yy + f_zz;
		}
	}
}

void perform_sph_brookshaw(particle* particles) {
	singleton_geometry *grid_geom  = get_singleton_geometry();

	for (unsigned int b = 0; b < grid_geom->num_cell; b++) {
		// finding neighbors
		int gi = 0;
		int gj = 0;
		int gk = 0;
		grid_unhash(b,&gi,&gj,&gk);

		int low_i  = gi-1 < 0 ? 0 : gi-1;
		int low_j  = gj-1 < 0 ? 0 : gj-1;
		int low_k  = gk-1 < 0 ? 0 : gk-1;
		int high_i = gi+2 > grid_geom->nx ? grid_geom->nx : gi+2;
		int high_j = gj+2 > grid_geom->ny ? grid_geom->ny : gj+2;
		int high_k = gk+2 > grid_geom->nz ? grid_geom->nz : gk+2;


		for (unsigned int i = grid_geom->cell_start[b]; i < grid_geom->cell_end[b]; i++) {

			if (particles[i].bnd) continue;

			particle part 	= particles[i];
			double xi 		= part.px;
			double yi 		= part.py;
			double zi 		= part.pz;
			double hi 		= part.h;
			double fi       = part.f;

			double f_x  = 0.; double f_y  = 0.; double f_z  = 0.;
			double f_xx = 0.; double f_yy = 0.; double f_zz = 0.;

			double laplacian  = 0.;

			// loop over neighbors
			for (unsigned int ni = low_i; ni < high_i; ni++) {
				for (unsigned int nj = low_j; nj < high_j; nj++) {
					for (unsigned int nk = low_k; nk < high_k; nk++) {

						unsigned int ID = ni*grid_geom->ny*grid_geom->nz + nj*grid_geom->nz + nk;

						for (unsigned int j = grid_geom->cell_start[ID]; j < grid_geom->cell_end[ID]; j++) {

							particle part 	= particles[j];
							double xj 		= part.px;
							double yj 		= part.py;
							double zj 		= part.pz;
							double rhoj     = part.rho;
							double m        = part.m;
							double fj       = part.f;

							double xij = xi-xj;
							double yij = yi-yj;
							double zij = zi-zj;
							double rij = sqrt(xij*xij + yij*yij + zij*zij);

							if (rij < 1e-8) continue;

							double eijx = xij/rij;
							double eijy = yij/rij;
							double eijz = zij/rij;
							double rij1 = 1./rij;

							kernel_result w = _kernel(xi,yi,zi,xj,yj,zj,hi);

							f_x  += (fi + fj)*w.w_x* m/rhoj;
							f_y  += (fi + fj)*w.w_y* m/rhoj;
							f_z  += (fi + fj)*w.w_z* m/rhoj;
							f_xx += (fj + fi)*w.w_xx*m/rhoj;
							f_yy += (fj + fi)*w.w_yy*m/rhoj;
							f_zz += (fj + fi)*w.w_zz*m/rhoj;

							laplacian += 2.0*(m/rhoj)*(fi-fj)*rij1*(eijx*w.w_x + eijy*w.w_y + eijz*w.w_z);
						}
					}
				}
			}
			particles[i].f_x = f_x;
			particles[i].f_y = f_y;
			particles[i].f_z = f_z;
			particles[i].LaplF = laplacian;
		}
	}
}

void sph_init(kernel kernel_fun) {
	_kernel = kernel_fun;
}
