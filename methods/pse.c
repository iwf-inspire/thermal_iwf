#include "pse.h"

void perform_pse(particle* particles) {
	singleton_geometry *grid_geom  = get_singleton_geometry();

	for (unsigned int b = 0; b < grid_geom->num_cell; b++) {

		//finding neighbors
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
			if (particles[i].bnd || particles[i].label==HEAT_SINK) continue;

			double xi  = particles[i].px;
			double yi  = particles[i].py;
			double zi  = particles[i].pz;
			double fi  = particles[i].f;
			double h   = particles[i].h;
			double h2  = h*h;
			double h3  = h*h2;

			double f_x  = 0.;
			double f_y  = 0.;
			double f_z  = 0.;
			double f_lapl3 = 0.;

			// loop over neighbors
			for (unsigned int ni = low_i; ni < high_i; ni++) {
				for (unsigned int nj = low_j; nj < high_j; nj++) {
					for (unsigned int nk = low_k; nk < high_k; nk++) {

						unsigned int ID = ni*grid_geom->ny*grid_geom->nz + nj*grid_geom->nz + nk;

						for (unsigned int j = grid_geom->cell_start[ID]; j < grid_geom->cell_end[ID]; j++) {

							if (particles[i].label==INLET && particles[j].label!=INLET) continue;
							if (particles[i].label==CONNECTOR_TOP && particles[j].label==INLET) continue;
							if (particles[i].label==CONNECTOR_INLET && particles[j].label==INLET) continue;
							if (particles[i].label==CONNECTOR_TOP && particles[j].label==CONNECTOR_INLET) continue;
							if (particles[i].label==CONNECTOR_INLET && particles[j].label==CONNECTOR_TOP) continue;
							if (particles[i].label==BAR && particles[j].label==CONNECTOR_BAR_CENTER) continue;
							if (particles[i].label==CONNECTOR_BAR_CENTER && particles[j].label==BAR) continue;

							double xj   = particles[j].px;
							double yj   = particles[j].py;
							double zj   = particles[j].pz;
							double fj   = particles[j].f;
							double m    = particles[j].m;
							double rhoj = particles[j].rho;

							double xij = xi-xj;
							double yij = yi-yj;
							double zij = zi-zj;

							double r = sqrt(xij*xij + yij*yij + zij*zij);

							// 2nd order kernel -- 1st spatial derivatives
							double w1x_pse = -2.*xij/(h2*M_PI*sqrt(M_PI))*exp(-r*r/(h2));
							double w1y_pse = -2.*yij/(h2*M_PI*sqrt(M_PI))*exp(-r*r/(h2));
							double w1z_pse = -2.*zij/(h2*M_PI*sqrt(M_PI))*exp(-r*r/(h2));

							// 2nd order kernel -- Laplacian
							double w3_pse  = +4./(h2*M_PI*sqrt(M_PI))*exp(-r*r/(h2));

							f_x += (fj + fi)*w1x_pse*m/rhoj/h3;
							f_y += (fj + fi)*w1y_pse*m/rhoj/h3;
							f_z += (fj + fi)*w1z_pse*m/rhoj/h3;

							f_lapl3 += (fj - fi)*w3_pse*m/rhoj/(h3);
						}
					}
				}
			}
			particles[i].f_x   = f_x;
			particles[i].f_y   = f_y;
			particles[i].f_z   = f_z;
			particles[i].LaplF = f_lapl3;
		}
	}
}
