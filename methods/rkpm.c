#include "rkpm.h"

static double *m0;
static double *m1;
static double *m2;

static double *m0_x;
static double *m1_x;
static double *m2_x;
static double *m0_y;
static double *m1_y;
static double *m2_y;
static double *m0_z;
static double *m1_z;
static double *m2_z;

static double *m0_xx;
static double *m1_xx;
static double *m2_xx;
static double *m0_yy;
static double *m1_yy;
static double *m2_yy;
static double *m0_zz;
static double *m1_zz;
static double *m2_zz;

static kernel _kernel;

void rkpm_wipe_out() {
	free(m0);
	free(m1);
	free(m2);

	free(m0_x);
	free(m1_x);
	free(m2_x);
	free(m0_y);
	free(m1_y);
	free(m2_y);
	free(m0_z);
	free(m1_z);
	free(m2_z);

	free(m0_xx);
	free(m1_xx);
	free(m2_xx);
	free(m0_yy);
	free(m1_yy);
	free(m2_yy);
	free(m0_zz);
	free(m1_zz);
	free(m2_zz);
}

void rkpm_compute_correction_terms(particle* particles) {
	singleton_geometry* grid_geom  =  get_singleton_geometry();

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

			double xi = particles[i].px;
			double yi = particles[i].py;
			double zi = particles[i].pz;
			double hi = particles[i].h;

			double m0i       = 0.f;
			double m1i[3]    = {0.f,0.f,0.f};
			double m2i[9]    = {0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f};

			double m0i_x     = 0.f;
			double m1i_x[3]  = {0.f,0.f,0.f};
			double m2i_x[9]  = {0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f};

			double m0i_y     = 0.f;
			double m1i_y[3]  = {0.f,0.f,0.f};
			double m2i_y[9]  = {0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f};

			double m0i_z     = 0.f;
			double m1i_z[3]  = {0.f,0.f,0.f};
			double m2i_z[9]  = {0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f};

			double m0i_xx    = 0.f;
			double m1i_xx[3] = {0.f,0.f,0.f};
			double m2i_xx[9] = {0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f};

			double m0i_yy    = 0.f;
			double m1i_yy[3] = {0.f,0.f,0.f};
			double m2i_yy[9] = {0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f};

			double m0i_zz    = 0.f;
			double m1i_zz[3] = {0.f,0.f,0.f};
			double m2i_zz[9] = {0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f};

			// loop over neighbors
			for (unsigned int ni = low_i; ni < high_i; ni++) {
				for (unsigned int nj = low_j; nj < high_j; nj++) {
					for (unsigned int nk = low_k; nk < high_k; nk++) {

						unsigned int ID = ni*grid_geom->ny*grid_geom->nz + nj*grid_geom->nz + nk;

						for (unsigned int j = grid_geom->cell_start[ID]; j < grid_geom->cell_end[ID]; j++) {

							double xj   = particles[j].px;
							double yj   = particles[j].py;
							double zj   = particles[j].pz;
							double rhoj = particles[j].rho;
							double m    = particles[j].m;

							double rx   = xi-xj;
							double ry   = yi-yj;
							double rz   = zi-zj;

							double r_xx = 1.;
							double r_xy = 0.;
							double r_xz = 0.;

							double r_yx = 0.;
							double r_yy = 1.;
							double r_yz = 0.;

							double r_zx = 0.;
							double r_zy = 0.;
							double r_zz = 1.;

							kernel_result w = _kernel(xi,yi,zi,xj,yj,zj,hi);

							// moment
							m0i      += w.w*m/rhoj;

							m1i[0]   += rx*w.w*m/rhoj;
							m1i[1]   += ry*w.w*m/rhoj;
							m1i[2]   += rz*w.w*m/rhoj;

							m2i[0]   += rx*rx*w.w*m/rhoj;
							m2i[1]   += rx*ry*w.w*m/rhoj;
							m2i[2]   += rx*rz*w.w*m/rhoj;
							m2i[3]   += ry*rx*w.w*m/rhoj;
							m2i[4]   += ry*ry*w.w*m/rhoj;
							m2i[5]   += ry*rz*w.w*m/rhoj;
							m2i[6]   += rz*rx*w.w*m/rhoj;
							m2i[7]   += rz*ry*w.w*m/rhoj;
							m2i[8]   += rz*rz*w.w*m/rhoj;

							// xder
							m0i_x    += w.w_x*m/rhoj;

							m1i_x[0] += w.w_x*rx*m/rhoj + w.w*r_xx*m/rhoj;
							m1i_x[1] += w.w_x*ry*m/rhoj + w.w*r_xy*m/rhoj;
							m1i_x[2] += w.w_x*rz*m/rhoj + w.w*r_xz*m/rhoj;

							m2i_x[0] += rx*rx*w.w_x*m/rhoj + w.w*(r_xx*rx + rx*r_xx)*m/rhoj;
							m2i_x[1] += rx*ry*w.w_x*m/rhoj + w.w*(r_xx*ry + rx*r_xy)*m/rhoj;
							m2i_x[2] += rx*rz*w.w_x*m/rhoj + w.w*(r_xx*rz + rx*r_xz)*m/rhoj;
							m2i_x[3] += ry*rx*w.w_x*m/rhoj + w.w*(r_xy*rx + ry*r_xx)*m/rhoj;
							m2i_x[4] += ry*ry*w.w_x*m/rhoj + w.w*(r_xy*ry + ry*r_xy)*m/rhoj;
							m2i_x[5] += ry*rz*w.w_x*m/rhoj + w.w*(r_xy*rz + ry*r_xz)*m/rhoj;
							m2i_x[6] += rz*rx*w.w_x*m/rhoj + w.w*(r_xz*rx + rz*r_xx)*m/rhoj;
							m2i_x[7] += rz*ry*w.w_x*m/rhoj + w.w*(r_xz*ry + rz*r_xy)*m/rhoj;
							m2i_x[8] += rz*rz*w.w_x*m/rhoj + w.w*(r_xz*rz + rz*r_xz)*m/rhoj;

							// yder
							m0i_y    += w.w_y*m/rhoj;

							m1i_y[0] += w.w_y*rx*m/rhoj + w.w*r_yx*m/rhoj;
							m1i_y[1] += w.w_y*ry*m/rhoj + w.w*r_yy*m/rhoj;
							m1i_y[2] += w.w_y*rz*m/rhoj + w.w*r_yz*m/rhoj;

							m2i_y[0] += rx*rx*w.w_y*m/rhoj + w.w*(r_yx*rx + rx*r_yx)*m/rhoj;
							m2i_y[1] += rx*ry*w.w_y*m/rhoj + w.w*(r_yx*ry + rx*r_yy)*m/rhoj;
							m2i_y[2] += rx*rz*w.w_y*m/rhoj + w.w*(r_yx*rz + rx*r_yz)*m/rhoj;
							m2i_y[3] += ry*rx*w.w_y*m/rhoj + w.w*(r_yy*rx + ry*r_yx)*m/rhoj;
							m2i_y[4] += ry*ry*w.w_y*m/rhoj + w.w*(r_yy*ry + ry*r_yy)*m/rhoj;
							m2i_y[5] += ry*rz*w.w_y*m/rhoj + w.w*(r_yy*rz + ry*r_yz)*m/rhoj;
							m2i_y[6] += rz*rx*w.w_y*m/rhoj + w.w*(r_yz*rx + rz*r_yx)*m/rhoj;
							m2i_y[7] += rz*ry*w.w_y*m/rhoj + w.w*(r_yz*ry + rz*r_yy)*m/rhoj;
							m2i_y[8] += rz*rz*w.w_y*m/rhoj + w.w*(r_yz*rz + rz*r_yz)*m/rhoj;

							// zder
							m0i_z    += w.w_z*m/rhoj;

							m1i_z[0] += w.w_z*rx*m/rhoj + w.w*r_zx*m/rhoj;
							m1i_z[1] += w.w_z*ry*m/rhoj + w.w*r_zy*m/rhoj;
							m1i_z[2] += w.w_z*rz*m/rhoj + w.w*r_zz*m/rhoj;

							m2i_z[0] += rx*rx*w.w_z*m/rhoj + w.w*(r_zx*rx + rx*r_zx)*m/rhoj;
							m2i_z[1] += rx*ry*w.w_z*m/rhoj + w.w*(r_zx*ry + rx*r_zy)*m/rhoj;
							m2i_z[2] += rx*rz*w.w_z*m/rhoj + w.w*(r_zx*rz + rx*r_zz)*m/rhoj;
							m2i_z[3] += ry*rx*w.w_z*m/rhoj + w.w*(r_zy*rx + ry*r_zx)*m/rhoj;
							m2i_z[4] += ry*ry*w.w_z*m/rhoj + w.w*(r_zy*ry + ry*r_zy)*m/rhoj;
							m2i_z[5] += ry*rz*w.w_z*m/rhoj + w.w*(r_zy*rz + ry*r_zz)*m/rhoj;
							m2i_z[6] += rz*rx*w.w_z*m/rhoj + w.w*(r_zz*rx + rz*r_zx)*m/rhoj;
							m2i_z[7] += rz*ry*w.w_z*m/rhoj + w.w*(r_zz*ry + rz*r_zy)*m/rhoj;
							m2i_z[8] += rz*rz*w.w_z*m/rhoj + w.w*(r_zz*rz + rz*r_zz)*m/rhoj;

							// xxder
							m0i_xx    += w.w_xx*m/rhoj;

							m1i_xx[0] += (2.f*r_xx*w.w_x + rx*w.w_xx)*m/rhoj;
							m1i_xx[1] += (2.f*r_xy*w.w_x + ry*w.w_xx)*m/rhoj;
							m1i_xx[2] += (2.f*r_xz*w.w_x + rz*w.w_xx)*m/rhoj;

							m2i_xx[0] += (rx*rx*w.w_xx + 4.*rx*r_xx*w.w_x + 2.*w.w*r_xx*r_xx)*m/rhoj;
							m2i_xx[1] += (rx*ry*w.w_xx + 4.*rx*r_xy*w.w_x + 2.*w.w*r_xx*r_xy)*m/rhoj;
							m2i_xx[2] += (rx*rz*w.w_xx + 4.*rx*r_xz*w.w_x + 2.*w.w*r_xx*r_xz)*m/rhoj;
							m2i_xx[3] += (ry*rx*w.w_xx + 4.*ry*r_xx*w.w_x + 2.*w.w*r_xy*r_xx)*m/rhoj;
							m2i_xx[4] += (ry*ry*w.w_xx + 4.*ry*r_xy*w.w_x + 2.*w.w*r_xy*r_xy)*m/rhoj;
							m2i_xx[5] += (ry*rz*w.w_xx + 4.*ry*r_xz*w.w_x + 2.*w.w*r_xy*r_xz)*m/rhoj;
							m2i_xx[6] += (rz*rx*w.w_xx + 4.*rz*r_xx*w.w_x + 2.*w.w*r_xz*r_xx)*m/rhoj;
							m2i_xx[7] += (rz*ry*w.w_xx + 4.*rz*r_xy*w.w_x + 2.*w.w*r_xz*r_xy)*m/rhoj;
							m2i_xx[8] += (rz*rz*w.w_xx + 4.*rz*r_xz*w.w_x + 2.*w.w*r_xz*r_xz)*m/rhoj;

							// yyder
							m0i_yy    += w.w_yy*m/rhoj;

							m1i_yy[0] += (2.f*r_yx*w.w_y + rx*w.w_yy)*m/rhoj;
							m1i_yy[1] += (2.f*r_yy*w.w_y + ry*w.w_yy)*m/rhoj;
							m1i_yy[2] += (2.f*r_yz*w.w_y + rz*w.w_yy)*m/rhoj;

							m2i_yy[0] += (rx*rx*w.w_yy + 4.*rx*r_yx*w.w_y + 2.*w.w*r_yx*r_yx)*m/rhoj;
							m2i_yy[1] += (rx*ry*w.w_yy + 4.*rx*r_yy*w.w_y + 2.*w.w*r_yx*r_yy)*m/rhoj;
							m2i_yy[2] += (rx*rz*w.w_yy + 4.*rx*r_yz*w.w_y + 2.*w.w*r_yx*r_yz)*m/rhoj;
							m2i_yy[3] += (ry*rx*w.w_yy + 4.*ry*r_yx*w.w_y + 2.*w.w*r_yy*r_yx)*m/rhoj;
							m2i_yy[4] += (ry*ry*w.w_yy + 4.*ry*r_yy*w.w_y + 2.*w.w*r_yy*r_yy)*m/rhoj;
							m2i_yy[5] += (ry*rz*w.w_yy + 4.*ry*r_yz*w.w_y + 2.*w.w*r_yy*r_yz)*m/rhoj;
							m2i_yy[6] += (rz*rx*w.w_yy + 4.*rz*r_yx*w.w_y + 2.*w.w*r_yz*r_yx)*m/rhoj;
							m2i_yy[7] += (rz*ry*w.w_yy + 4.*rz*r_yy*w.w_y + 2.*w.w*r_yz*r_yy)*m/rhoj;
							m2i_yy[8] += (rz*rz*w.w_yy + 4.*rz*r_yz*w.w_y + 2.*w.w*r_yz*r_yz)*m/rhoj;

							// zzder
							m0i_zz    += w.w_zz*m/rhoj;

							m1i_zz[0] += (2.f*r_zx*w.w_z + rx*w.w_zz)*m/rhoj;
							m1i_zz[1] += (2.f*r_zy*w.w_z + ry*w.w_zz)*m/rhoj;
							m1i_zz[2] += (2.f*r_zz*w.w_z + rz*w.w_zz)*m/rhoj;

							m2i_zz[0] += (rx*rx*w.w_zz + 4.*rx*r_zx*w.w_z + 2.*w.w*r_zx*r_zx)*m/rhoj;
							m2i_zz[1] += (rx*ry*w.w_zz + 4.*rx*r_zy*w.w_z + 2.*w.w*r_zx*r_zy)*m/rhoj;
							m2i_zz[2] += (rx*rz*w.w_zz + 4.*rx*r_zz*w.w_z + 2.*w.w*r_zx*r_zz)*m/rhoj;
							m2i_zz[3] += (ry*rx*w.w_zz + 4.*ry*r_zx*w.w_z + 2.*w.w*r_zy*r_zx)*m/rhoj;
							m2i_zz[4] += (ry*ry*w.w_zz + 4.*ry*r_zy*w.w_z + 2.*w.w*r_zy*r_zy)*m/rhoj;
							m2i_zz[5] += (ry*rz*w.w_zz + 4.*ry*r_zz*w.w_z + 2.*w.w*r_zy*r_zz)*m/rhoj;
							m2i_zz[6] += (rz*rx*w.w_zz + 4.*rz*r_zx*w.w_z + 2.*w.w*r_zz*r_zx)*m/rhoj;
							m2i_zz[7] += (rz*ry*w.w_zz + 4.*rz*r_zy*w.w_z + 2.*w.w*r_zz*r_zy)*m/rhoj;
							m2i_zz[8] += (rz*rz*w.w_zz + 4.*rz*r_zz*w.w_z + 2.*w.w*r_zz*r_zz)*m/rhoj;
						}
					}
				}
			}


			m0[i]    = m0i;
			m0_x[i]  = m0i_x;
			m0_y[i]  = m0i_y;
			m0_z[i]  = m0i_z;
			m0_xx[i] = m0i_xx;
			m0_yy[i] = m0i_yy;
			m0_zz[i] = m0i_zz;

			memcpy(m1    + 3*i,m1i,    sizeof(double)*3);
			memcpy(m1_x  + 3*i,m1i_x,  sizeof(double)*3);
			memcpy(m1_y  + 3*i,m1i_y,  sizeof(double)*3);
			memcpy(m1_z  + 3*i,m1i_z,  sizeof(double)*3);
			memcpy(m1_xx + 3*i,m1i_xx, sizeof(double)*3);
			memcpy(m1_yy + 3*i,m1i_yy, sizeof(double)*3);
			memcpy(m1_zz + 3*i,m1i_zz, sizeof(double)*3);

			memcpy(m2    + 9*i,m2i,    sizeof(double)*9);
			memcpy(m2_x  + 9*i,m2i_x,  sizeof(double)*9);
			memcpy(m2_y  + 9*i,m2i_y,  sizeof(double)*9);
			memcpy(m2_z  + 9*i,m2i_z,  sizeof(double)*9);
			memcpy(m2_xx + 9*i,m2i_xx, sizeof(double)*9);
			memcpy(m2_yy + 9*i,m2i_yy, sizeof(double)*9);
			memcpy(m2_zz + 9*i,m2i_zz, sizeof(double)*9);
		}
	}
}

void perform_rkpm(particle* particles) {
	singleton_geometry  *grid_geom  = get_singleton_geometry();

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

			double xi  = particles[i].px;
			double yi  = particles[i].py;
			double zi  = particles[i].pz;
			double hi  = particles[i].h;
			double fi  = particles[i].f;

			double f_x   = 0.f;
			double f_y   = 0.f;
			double f_z   = 0.f;
			double f_xx  = 0.f;
			double f_yy  = 0.f;
			double f_zz  = 0.f;

			double M[16] = { m0[i    ], m1[3*i+0], m1[3*i+1], m1[3*i+2],
					m1[3*i+0], m2[9*i+0], m2[9*i+1], m2[9*i+2],
					m1[3*i+1], m2[9*i+3], m2[9*i+4], m2[9*i+5],
					m1[3*i+2], m2[9*i+6], m2[9*i+7], m2[9*i+8]};

			double M_x[16] = { m0_x[i    ], m1_x[3*i+0], m1_x[3*i+1], m1_x[3*i+2],
					m1_x[3*i+0], m2_x[9*i+0], m2_x[9*i+1], m2_x[9*i+2],
					m1_x[3*i+1], m2_x[9*i+3], m2_x[9*i+4], m2_x[9*i+5],
					m1_x[3*i+2], m2_x[9*i+6], m2_x[9*i+7], m2_x[9*i+8]};

			double M_y[16] = { m0_y[i    ], m1_y[3*i+0], m1_y[3*i+1], m1_y[3*i+2],
					m1_y[3*i+0], m2_y[9*i+0], m2_y[9*i+1], m2_y[9*i+2],
					m1_y[3*i+1], m2_y[9*i+3], m2_y[9*i+4], m2_y[9*i+5],
					m1_y[3*i+2], m2_y[9*i+6], m2_y[9*i+7], m2_y[9*i+8]};

			double M_z[16] = { m0_z[i    ], m1_z[3*i+0], m1_z[3*i+1], m1_z[3*i+2],
					m1_z[3*i+0], m2_z[9*i+0], m2_z[9*i+1], m2_z[9*i+2],
					m1_z[3*i+1], m2_z[9*i+3], m2_z[9*i+4], m2_z[9*i+5],
					m1_z[3*i+2], m2_z[9*i+6], m2_z[9*i+7], m2_z[9*i+8]};

			double M_xx[16] = { m0_xx[i    ], m1_xx[3*i+0], m1_xx[3*i+1], m1_xx[3*i+2],
					m1_xx[3*i+0], m2_xx[9*i+0], m2_xx[9*i+1], m2_xx[9*i+2],
					m1_xx[3*i+1], m2_xx[9*i+3], m2_xx[9*i+4], m2_xx[9*i+5],
					m1_xx[3*i+2], m2_xx[9*i+6], m2_xx[9*i+7], m2_xx[9*i+8]};

			double M_yy[16] = { m0_yy[i    ], m1_yy[3*i+0], m1_yy[3*i+1], m1_yy[3*i+2],
					m1_yy[3*i+0], m2_yy[9*i+0], m2_yy[9*i+1], m2_yy[9*i+2],
					m1_yy[3*i+1], m2_yy[9*i+3], m2_yy[9*i+4], m2_yy[9*i+5],
					m1_yy[3*i+2], m2_yy[9*i+6], m2_yy[9*i+7], m2_yy[9*i+8]};

			double M_zz[16] = { m0_zz[i    ], m1_zz[3*i+0], m1_zz[3*i+1], m1_zz[3*i+2],
					m1_zz[3*i+0], m2_zz[9*i+0], m2_zz[9*i+1], m2_zz[9*i+2],
					m1_zz[3*i+1], m2_zz[9*i+3], m2_zz[9*i+4], m2_zz[9*i+5],
					m1_zz[3*i+2], m2_zz[9*i+6], m2_zz[9*i+7], m2_zz[9*i+8]};

			double P[4] = {1.0, 0.0, 0.0, 0.0};

			double C[4], C_x[4], C_y[4], C_z[4], C_xx[4], C_yy[4], C_zz[4];

			solve4x4(M, P, C);

			double rhs_1[4] = { -M_x[ 0]*C[0]  -M_x[ 1]*C[1]  -M_x[ 2]*C[2]  -M_x[ 3]*C[3],
					-M_x[ 4]*C[0]  -M_x[ 5]*C[1]  -M_x[ 6]*C[2]  -M_x[ 7]*C[3],
					-M_x[ 8]*C[0]  -M_x[ 9]*C[1]  -M_x[10]*C[2]  -M_x[11]*C[3],
					-M_x[12]*C[0]  -M_x[13]*C[1]  -M_x[14]*C[2]  -M_x[15]*C[3]};

			double rhs_2[4] = { -M_y[0]*C[0]  -M_y[1]*C[1]  -M_y[2]*C[2]  -M_y[3]*C[3],
					-M_y[4]*C[0]  -M_y[5]*C[1]  -M_y[6]*C[2]  -M_y[7]*C[3],
					-M_y[8]*C[0]  -M_y[9]*C[1]  -M_y[10]*C[2] -M_y[11]*C[3],
					-M_y[12]*C[0] -M_y[13]*C[1] -M_y[14]*C[2] -M_y[15]*C[3]};

			double rhs_3[4] = { -M_z[0]*C[0]  -M_z[1]*C[1]  -M_z[2]*C[2]  -M_z[3]*C[3],
					-M_z[4]*C[0]  -M_z[5]*C[1]  -M_z[6]*C[2]  -M_z[7]*C[3],
					-M_z[8]*C[0]  -M_z[9]*C[1]  -M_z[10]*C[2] -M_z[11]*C[3],
					-M_z[12]*C[0] -M_z[13]*C[1] -M_z[14]*C[2] -M_z[15]*C[3]};

			solve4x4(M, rhs_1, C_x);
			solve4x4(M, rhs_2, C_y);
			solve4x4(M, rhs_3, C_z);

			double rhs_a[4] = {-M_xx[0] *C[0] -M_xx[1] *C[1] -M_xx[2] *C[2] -M_xx[3] *C[3] - 2*(M_x[0] *C_x[0] + M_x[1] *C_x[1] + M_x[2] *C_x[2] + M_x[3] *C_x[3]),
					-M_xx[4] *C[0] -M_xx[5] *C[1] -M_xx[6] *C[2] -M_xx[7] *C[3] - 2*(M_x[4] *C_x[0] + M_x[5] *C_x[1] + M_x[6] *C_x[2] + M_x[7] *C_x[3]),
					-M_xx[8] *C[0] -M_xx[9] *C[1] -M_xx[10]*C[2] -M_xx[11]*C[3] - 2*(M_x[8] *C_x[0] + M_x[9] *C_x[1] + M_x[10]*C_x[2] + M_x[11]*C_x[3]),
					-M_xx[12]*C[0] -M_xx[13]*C[1] -M_xx[14]*C[2] -M_xx[15]*C[3] - 2*(M_x[12]*C_x[0] + M_x[13]*C_x[1] + M_x[14]*C_x[2] + M_x[15]*C_x[3])};

			double rhs_b[4] = {-M_yy[0] *C[0] -M_yy[1] *C[1] -M_yy[2] *C[2] -M_yy[3] *C[3] - 2*(M_y[0] *C_y[0] + M_y[1] *C_y[1] + M_y[2] *C_y[2] + M_y[3] *C_y[3]),
					-M_yy[4] *C[0] -M_yy[5] *C[1] -M_yy[6] *C[2] -M_yy[7] *C[3] - 2*(M_y[4] *C_y[0] + M_y[5] *C_y[1] + M_y[6] *C_y[2] + M_y[7] *C_y[3]),
					-M_yy[8] *C[0] -M_yy[9] *C[1] -M_yy[10]*C[2] -M_yy[11]*C[3] - 2*(M_y[8] *C_y[0] + M_y[9] *C_y[1] + M_y[10]*C_y[2] + M_y[11]*C_y[3]),
					-M_yy[12]*C[0] -M_yy[13]*C[1] -M_yy[14]*C[2] -M_yy[15]*C[3] - 2*(M_y[12]*C_y[0] + M_y[13]*C_y[1] + M_y[14]*C_y[2] + M_y[15]*C_y[3])};

			double rhs_c[4] = {-M_zz[0] *C[0] -M_zz[1] *C[1] -M_zz[2] *C[2] -M_zz[3] *C[3] - 2*(M_z[0] *C_z[0] + M_z[1] *C_z[1] + M_z[2] *C_z[2] + M_z[3] *C_z[3]),
					-M_zz[4] *C[0] -M_zz[5] *C[1] -M_zz[6] *C[2] -M_zz[7] *C[3] - 2*(M_z[4] *C_z[0] + M_z[5] *C_z[1] + M_z[6] *C_z[2] + M_z[7] *C_z[3]),
					-M_zz[8] *C[0] -M_zz[9] *C[1] -M_zz[10]*C[2] -M_zz[11]*C[3] - 2*(M_z[8] *C_z[0] + M_z[9] *C_z[1] + M_z[10]*C_z[2] + M_z[11]*C_z[3]),
					-M_zz[12]*C[0] -M_zz[13]*C[1] -M_zz[14]*C[2] -M_zz[15]*C[3] - 2*(M_z[12]*C_z[0] + M_z[13]*C_z[1] + M_z[14]*C_z[2] + M_z[15]*C_z[3])};

			solve4x4(M, rhs_a, C_xx);
			solve4x4(M, rhs_b, C_yy);
			solve4x4(M, rhs_c, C_zz);

			// loop over neighbors
			for (unsigned int ni = low_i; ni < high_i; ni++) {
				for (unsigned int nj = low_j; nj < high_j; nj++) {
					for (unsigned int nk = low_k; nk < high_k; nk++) {

						unsigned int ID = ni*grid_geom->ny*grid_geom->nz + nj*grid_geom->nz + nk;

						for (unsigned int j = grid_geom->cell_start[ID]; j < grid_geom->cell_end[ID]; j++) {

							double xj 		= particles[j].px;
							double yj 		= particles[j].py;
							double zj 		= particles[j].pz;
							double rhoj 	= particles[j].rho;
							double m        = particles[j].m;
							double fj       = particles[j].f;

							double rx   = xi-xj;
							double ry   = yi-yj;
							double rz   = zi-zj;

							double r_xx = 1.;
							double r_xy = 0.;
							double r_xz = 0.;

							double r_yx = 0.;
							double r_yy = 1.;
							double r_yz = 0.;

							double r_zx = 0.;
							double r_zy = 0.;
							double r_zz = 1.;

							kernel_result w = _kernel(xi,yi,zi,xj,yj,zj,hi);

							double Corr    = (rx*C[1] + ry*C[2] + rz*C[3]) + C[0];

							double Corr_x  = (rx*C_x[1] + ry*C_x[2] + rz*C_x[3]) + (r_xx*C[1] + r_xy*C[2] + r_xz*C[3]) + C_x[0];
							double Corr_y  = (rx*C_y[1] + ry*C_y[2] + rz*C_y[3]) + (r_yx*C[1] + r_yy*C[2] + r_yz*C[3]) + C_y[0];
							double Corr_z  = (rx*C_z[1] + ry*C_z[2] + rz*C_z[3]) + (r_zx*C[1] + r_zy*C[2] + r_zz*C[3]) + C_z[0];

							double Corr_xx = (rx*C_xx[1] + ry*C_xx[2] + rz*C_xx[3]) + 2*(r_xx*C_x[1] + r_xy*C_x[2] + r_xz*C_x[3]) + C_xx[0];
							double Corr_yy = (rx*C_yy[1] + ry*C_yy[2] + rz*C_yy[3]) + 2*(r_yx*C_y[1] + r_yy*C_y[2] + r_yz*C_y[3]) + C_yy[0];
							double Corr_zz = (rx*C_zz[1] + ry*C_zz[2] + rz*C_zz[3]) + 2*(r_zx*C_z[1] + r_zy*C_z[2] + r_zz*C_z[3]) + C_zz[0];

							double Nj    = (Corr*w.w)*m/rhoj;
							double Nj_x  = (Corr_x*w.w + Corr*w.w_x)*m/rhoj;
							double Nj_y  = (Corr_y*w.w + Corr*w.w_y)*m/rhoj;
							double Nj_z  = (Corr_z*w.w + Corr*w.w_z)*m/rhoj;
							double Nj_xx = (Corr_xx*w.w + 2*Corr_x*w.w_x + Corr*w.w_xx)*m/rhoj;
							double Nj_yy = (Corr_yy*w.w + 2*Corr_y*w.w_y + Corr*w.w_yy)*m/rhoj;
							double Nj_zz = (Corr_zz*w.w + 2*Corr_z*w.w_z + Corr*w.w_zz)*m/rhoj;

							f_x  += (fj - fi)*Nj_x;
							f_y  += (fj - fi)*Nj_y;
							f_z  += (fj - fi)*Nj_z;
							f_xx += (fj - fi)*Nj_xx;
							f_yy += (fj - fi)*Nj_yy;
							f_zz += (fj - fi)*Nj_zz;
						}
					}
				}
			}
			particles[i].f_x   = f_x;
			particles[i].f_y   = f_y;
			particles[i].f_z   = f_z;
			particles[i].f_xx  = f_xx;
			particles[i].f_yy  = f_yy;
			particles[i].f_zz  = f_zz;
			particles[i].LaplF = f_xx + f_yy + f_zz;
		}
	}
}

void rkpm_init(unsigned int Ntot, kernel kernel_fun) {

	m0    = (double*) calloc(Ntot,   sizeof(double));		//scalar
	m1    = (double*) calloc(Ntot*3, sizeof(double));		//vector
	m2    = (double*) calloc(Ntot*9, sizeof(double));		//tensor

	m0_x  = (double*) calloc(Ntot,   sizeof(double));		//scalar
	m1_x  = (double*) calloc(Ntot*3, sizeof(double));		//vector
	m2_x  = (double*) calloc(Ntot*9, sizeof(double));		//tensor

	m0_y  = (double*) calloc(Ntot,   sizeof(double));		//scalar
	m1_y  = (double*) calloc(Ntot*3, sizeof(double));		//vector
	m2_y  = (double*) calloc(Ntot*9, sizeof(double));		//tensor

	m0_z  = (double*) calloc(Ntot,   sizeof(double));		//scalar
	m1_z  = (double*) calloc(Ntot*3, sizeof(double));		//vector
	m2_z  = (double*) calloc(Ntot*9, sizeof(double));		//tensor

	m0_xx = (double*) calloc(Ntot,   sizeof(double));		//scalar
	m1_xx = (double*) calloc(Ntot*3, sizeof(double));		//vector
	m2_xx = (double*) calloc(Ntot*9, sizeof(double));		//tensor

	m0_yy = (double*) calloc(Ntot,   sizeof(double));		//scalar
	m1_yy = (double*) calloc(Ntot*3, sizeof(double));		//vector
	m2_yy = (double*) calloc(Ntot*9, sizeof(double));		//tensor

	m0_zz = (double*) calloc(Ntot,   sizeof(double));		//scalar
	m1_zz = (double*) calloc(Ntot*3, sizeof(double));		//vector
	m2_zz = (double*) calloc(Ntot*9, sizeof(double));		//tensor

	_kernel = kernel_fun;
}
