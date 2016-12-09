#include "cspm.h"

static double *A;
static double *B;

static kernel _kernel;

void cspm_wipe_out() {
	free(A);
	free(B);
}

void cspm_compute_correction_terms(particle* particles) {
	singleton_geometry* grid_geom = get_singleton_geometry();

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

			if (particles[i].bnd) continue;

			double xi = particles[i].px;
			double yi = particles[i].py;
			double zi = particles[i].pz;
			double hi = particles[i].h;

			double A11 = 0.f; double A12 = 0.f; double A13 = 0.f;
			double A21 = 0.f; double A22 = 0.f; double A23 = 0.f;
			double A31 = 0.f; double A32 = 0.f; double A33 = 0.f;

			double B11 = 0.f; double B12 = 0.f; double B13 = 0.f;
			double B14 = 0.f; double B15 = 0.f; double B16 = 0.f;
			double B21 = 0.f; double B22 = 0.f; double B23 = 0.f;
			double B24 = 0.f; double B25 = 0.f; double B26 = 0.f;
			double B31 = 0.f; double B32 = 0.f; double B33 = 0.f;
			double B34 = 0.f; double B35 = 0.f; double B36 = 0.f;
			double B41 = 0.f; double B42 = 0.f; double B43 = 0.f;
			double B44 = 0.f; double B45 = 0.f; double B46 = 0.f;
			double B51 = 0.f; double B52 = 0.f; double B53 = 0.f;
			double B54 = 0.f; double B55 = 0.f; double B56 = 0.f;
			double B61 = 0.f; double B62 = 0.f; double B63 = 0.f;
			double B64 = 0.f; double B65 = 0.f; double B66 = 0.f;

			//loop over neighbors
			for (unsigned int ni = low_i; ni < high_i; ni++) {
				for (unsigned int nj = low_j; nj < high_j; nj++) {
					for (unsigned int nk = low_k; nk < high_k; nk++) {

						unsigned int ID = ni*grid_geom->ny*grid_geom->nz + nj*grid_geom->nz + nk;

						for (unsigned int j = grid_geom->cell_start[ID]; j < grid_geom->cell_end[ID]; j++) {

							double xj    = particles[j].px;
							double yj    = particles[j].py;
							double zj    = particles[j].pz;
							double rhoj  = particles[j].rho;
							double m     = particles[j].m;

							double rx   = xj-xi;
							double ry   = yj-yi;
							double rz   = zj-zi;

							double rij = sqrt(rx*rx + ry*ry + rz*rz);
							double q = rij/hi;

							kernel_result w = _kernel(xi,yi,zi,xj,yj,zj,hi);

							// A: 1st derivatives (3x3)
							A11 += rx*w.w_x*m/rhoj; A12 += ry*w.w_x*m/rhoj; A13 += rz*w.w_x*m/rhoj;
							A21 += rx*w.w_y*m/rhoj; A22 += ry*w.w_y*m/rhoj; A23 += rz*w.w_y*m/rhoj;
							A31 += rx*w.w_z*m/rhoj; A32 += ry*w.w_z*m/rhoj; A33 += rz*w.w_z*m/rhoj;

							// B: 2nd derivatives (6x6)
							B11 += rx*rx*w.w_xx*m/rhoj;
							B12 += rx*ry*w.w_xx*m/rhoj;
							B13 += rx*rz*w.w_xx*m/rhoj;
							B14 += ry*ry*w.w_xx*m/rhoj;
							B15 += ry*rz*w.w_xx*m/rhoj;
							B16 += rz*rz*w.w_xx*m/rhoj;

							B21 += rx*rx*w.w_xy*m/rhoj;
							B22 += rx*ry*w.w_xy*m/rhoj;
							B23 += rx*rz*w.w_xy*m/rhoj;
							B24 += ry*ry*w.w_xy*m/rhoj;
							B25 += ry*rz*w.w_xy*m/rhoj;
							B26 += rz*rz*w.w_xy*m/rhoj;

							B31 += rx*rx*w.w_xz*m/rhoj;
							B32 += rx*ry*w.w_xz*m/rhoj;
							B33 += rx*rz*w.w_xz*m/rhoj;
							B34 += ry*ry*w.w_xz*m/rhoj;
							B35 += ry*rz*w.w_xz*m/rhoj;
							B36 += rz*rz*w.w_xz*m/rhoj;

							B41 += rx*rx*w.w_yy*m/rhoj;
							B42 += rx*ry*w.w_yy*m/rhoj;
							B43 += rx*rz*w.w_yy*m/rhoj;
							B44 += ry*ry*w.w_yy*m/rhoj;
							B45 += ry*rz*w.w_yy*m/rhoj;
							B46 += rz*rz*w.w_yy*m/rhoj;

							B51 += rx*rx*w.w_yz*m/rhoj;
							B52 += rx*ry*w.w_yz*m/rhoj;
							B53 += rx*rz*w.w_yz*m/rhoj;
							B54 += ry*ry*w.w_yz*m/rhoj;
							B55 += ry*rz*w.w_yz*m/rhoj;
							B56 += rz*rz*w.w_yz*m/rhoj;

							B61 += rx*rx*w.w_zz*m/rhoj;
							B62 += rx*ry*w.w_zz*m/rhoj;
							B63 += rx*rz*w.w_zz*m/rhoj;
							B64 += ry*ry*w.w_zz*m/rhoj;
							B65 += ry*rz*w.w_zz*m/rhoj;
							B66 += rz*rz*w.w_zz*m/rhoj;
						}
					}
				}
			}

			double Ai[9] = {A11,A12,A13,
					A21,A22,A23,
					A31,A32,A33};

			double Bi[36] = {0.5*B11,B12,B13,0.5*B14,B15,0.5*B16,
					0.5*B21,B22,B23,0.5*B24,B25,0.5*B26,
					0.5*B31,B32,B33,0.5*B34,B35,0.5*B36,
					0.5*B41,B42,B43,0.5*B44,B45,0.5*B46,
					0.5*B51,B52,B53,0.5*B54,B55,0.5*B56,
					0.5*B61,B62,B63,0.5*B64,B65,0.5*B66};

			// chunks of correction terms for entire particles
			memcpy(A +  9*i, Ai, sizeof(double)*9);
			memcpy(B + 36*i, Bi, sizeof(double)*36);
		}
	}
}

void perform_cspm(particle* particles) {
	singleton_geometry* grid_geom = get_singleton_geometry();

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

			if (particles[i].bnd) continue;

			double xi  = particles[i].px;
			double yi  = particles[i].py;
			double zi  = particles[i].pz;
			double hi  = particles[i].h;
			double fi  = particles[i].f;

			double nom = 0.f;
			double den = 0.f;

			double F1t  = 0.f;
			double F2t  = 0.f;
			double F3t  = 0.f;

			//loop over neighbors
			for (unsigned int ni = low_i; ni < high_i; ni++) {
				for (unsigned int nj = low_j; nj < high_j; nj++) {
					for (unsigned int nk = low_k; nk < high_k; nk++) {

						unsigned int ID = ni*grid_geom->ny*grid_geom->nz + nj*grid_geom->nz + nk;

						for (unsigned int j = grid_geom->cell_start[ID]; j < grid_geom->cell_end[ID]; j++) {

							//if (particles[i].label==INLET && particles[j].label!=INLET) continue;
							//if (particles[i].label==INLET && particles[j].label==CONNECTOR_TOP) continue;
							//if (particles[i].label==CONNECTOR_TOP && particles[j].label==INLET) continue;
							//if (particles[i].label==CONNECTOR_INLET && particles[j].label==INLET) continue;
							//if (particles[i].label==CONNECTOR_TOP && particles[j].label==CONNECTOR_INLET) continue;
							//if (particles[i].label==CONNECTOR_INLET && particles[j].label==CONNECTOR_TOP) continue;
							//if (particles[i].label==BAR && particles[j].label==CONNECTOR_BAR_CENTER) continue;
							//if (particles[i].label==CONNECTOR_BAR_CENTER && particles[j].label==BAR) continue;

							double xj   = particles[j].px;
							double yj   = particles[j].py;
							double zj   = particles[j].pz;
							double rhoj = particles[j].rho;
							double m    = particles[j].m;
							double fj   = particles[j].f;

							//if(particles[j].bnd) fj = fi;

							kernel_result w = _kernel(xi,yi,zi,xj,yj,zj,hi);

							nom +=  fj*w.w*m/rhoj;
							den +=     w.w*m/rhoj;

							F1t += (fj - fi)*w.w_x*m/rhoj;
							F2t += (fj - fi)*w.w_y*m/rhoj;
							F3t += (fj - fi)*w.w_z*m/rhoj;
						}
					}
				}
			}

			// solving 3x3 system of simultaneous linear equations for 1st derivatives
			double Ai[9];
			memcpy(Ai,A + 9*i,sizeof(double)*9);
			double rhs1[3] = {F1t,F2t,F3t};
			double ders1[3];

			solve3x3(Ai,rhs1,ders1);
			particles[i].f_x = ders1[0];
			particles[i].f_y = ders1[1];
			particles[i].f_z = ders1[2];

			// rhs2 terms for 2nd derivatives
			double Fxx = 0.f; double Bxxx = 0.f; double Bxxy = 0.f; double Bxxz=0.f;
			double Fxy = 0.f; double Bxyx = 0.f; double Bxyy = 0.f; double Bxyz=0.f;
			double Fxz = 0.f; double Bxzx = 0.f; double Bxzy = 0.f; double Bxzz=0.f;
			double Fyy = 0.f; double Byyx = 0.f; double Byyy = 0.f; double Byyz=0.f;
			double Fyz = 0.f; double Byzx = 0.f; double Byzy = 0.f; double Byzz=0.f;
			double Fzz = 0.f; double Bzzx = 0.f; double Bzzy = 0.f; double Bzzz=0.f;

			//loop over neighbors
			for (unsigned int ni = low_i; ni < high_i; ni++) {
				for (unsigned int nj = low_j; nj < high_j; nj++) {
					for (unsigned int nk = low_k; nk < high_k; nk++) {

						unsigned int ID = ni*grid_geom->ny*grid_geom->nz + nj*grid_geom->nz + nk;

						for (unsigned int j = grid_geom->cell_start[ID]; j < grid_geom->cell_end[ID]; j++) {

							//if (particles[i].label==INLET && particles[j].label!=INLET) continue;
							//if (particles[i].label==INLET && particles[j].label==CONNECTOR_TOP) continue;
							//if (particles[i].label==CONNECTOR_TOP && particles[j].label==INLET) continue;
							//if (particles[i].label==CONNECTOR_INLET && particles[j].label==INLET) continue;
							//if (particles[i].label==CONNECTOR_TOP && particles[j].label==CONNECTOR_INLET) continue;
							//if (particles[i].label==CONNECTOR_INLET && particles[j].label==CONNECTOR_TOP) continue;
							//if (particles[i].label==BAR && particles[j].label==CONNECTOR_BAR_CENTER) continue;
							//if (particles[i].label==CONNECTOR_BAR_CENTER && particles[j].label==BAR) continue;

							double xj   = particles[j].px;
							double yj 	= particles[j].py;
							double zj   = particles[j].pz;
							double rhoj = particles[j].rho;
							double m    = particles[j].m;
							double fj   = particles[j].f;

							//if(particles[j].bnd) fj = fi;

							kernel_result w = _kernel(xi,yi,zi,xj,yj,zj,hi);

							Fxx += (fj-fi)*w.w_xx*m/rhoj;
							Fxy += (fj-fi)*w.w_xy*m/rhoj;
							Fxz += (fj-fi)*w.w_xz*m/rhoj;
							Fyy += (fj-fi)*w.w_yy*m/rhoj;
							Fyz += (fj-fi)*w.w_yz*m/rhoj;
							Fzz += (fj-fi)*w.w_zz*m/rhoj;

							Bxxx += (xj-xi)*w.w_xx*m/rhoj;
							Bxyx += (xj-xi)*w.w_xy*m/rhoj;
							Bxzx += (xj-xi)*w.w_xz*m/rhoj;
							Byyx += (xj-xi)*w.w_yy*m/rhoj;
							Byzx += (xj-xi)*w.w_yz*m/rhoj;
							Bzzx += (xj-xi)*w.w_zz*m/rhoj;

							Bxxy += (yj-yi)*w.w_xx*m/rhoj;
							Bxyy += (yj-yi)*w.w_xy*m/rhoj;
							Bxzy += (yj-yi)*w.w_xz*m/rhoj;
							Byyy += (yj-yi)*w.w_yy*m/rhoj;
							Byzy += (yj-yi)*w.w_yz*m/rhoj;
							Bzzy += (yj-yi)*w.w_zz*m/rhoj;

							Bxxz += (zj-zi)*w.w_xx*m/rhoj;
							Bxyz += (zj-zi)*w.w_xy*m/rhoj;
							Bxzz += (zj-zi)*w.w_xz*m/rhoj;
							Byyz += (zj-zi)*w.w_yy*m/rhoj;
							Byzz += (zj-zi)*w.w_yz*m/rhoj;
							Bzzz += (zj-zi)*w.w_zz*m/rhoj;
						}
					}
				}
			}

			// solving 6x6 system of simultaneous linear equations for 2nd derivatives
			double Bi[36];
			memcpy(Bi,B + 36*i,sizeof(double)*36);

			double fdx = particles[i].f_x;
			double fdy = particles[i].f_y;
			double fdz = particles[i].f_z;


			double rhs2[6] = {Fxx - Bxxx*fdx - Bxxy*fdy - Bxxz*fdz,
					Fxy - Bxyx*fdx - Bxyy*fdy - Bxyz*fdz,
					Fxz - Bxzx*fdx - Bxzy*fdy - Bxzz*fdz,
					Fyy - Byyx*fdx - Byyy*fdy - Byyz*fdz,
					Fyz - Byzx*fdx - Byzy*fdy - Byzz*fdz,
					Fzz - Bzzx*fdx - Bzzy*fdy - Bzzz*fdz};

			double ders2[6];
			solve6x6(Bi,rhs2,ders2);

			double f_xx = ders2[0];
			double f_xy = ders2[1];
			double f_xz = ders2[2];
			double f_yy = ders2[3];
			double f_yz = ders2[4];
			double f_zz = ders2[5];

			particles[i].LaplF = f_xx + f_yy + f_zz;
		}
	}
}

void cspm_set_kernel(kernel kernel_fun) {
	_kernel = kernel_fun;
}

void cspm_init(unsigned int Ntot, kernel kernel_fun) {
	//first derivatives
	A = (double*) calloc(Ntot*3*3, sizeof(double));
	//second derivatives
	B = (double*) calloc(Ntot*6*6, sizeof(double));

	_kernel = kernel_fun;
}
