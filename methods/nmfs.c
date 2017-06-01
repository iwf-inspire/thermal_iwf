#include "nmfs.h"

// normalization tensors for 1st and 2nd derivatives
static double *B, *Bhat;

static kernel _kernel;

void nmfs_wipe_out() {
	free(B);
	free(Bhat);
}

void nmfs_compute_correction_terms(particle* particles) {
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

			/*     ref. http://www.sciencedirect.com/science/article/pii/S0898122110009004
			 *     "Error estimation in smoothed particle hydrodynamics and a new scheme for second derivatives"
			 *
			 *	   Eq(40) B_hat (new renormalization tensor)
			 *
			 *     B_hat:[E + A.B.C] = -I
			 *     A.B.C = D
			 *     E + D = F
			 */

			double xi = particles[i].px;
			double yi = particles[i].py;
			double zi = particles[i].pz;
			double hi = particles[i].h;

			// 3rd order tensor in 3D -- 3x3x3 = 27 components
			double A111 = 0.f; double A121 = 0.f; double A131 = 0.f;
			double A211 = 0.f; double A221 = 0.f; double A231 = 0.f;
			double A311 = 0.f; double A321 = 0.f; double A331 = 0.f;
			double A112 = 0.f; double A122 = 0.f; double A132 = 0.f;
			double A212 = 0.f; double A222 = 0.f; double A232 = 0.f;
			double A312 = 0.f; double A322 = 0.f; double A332 = 0.f;
			double A113 = 0.f; double A123 = 0.f; double A133 = 0.f;
			double A213 = 0.f; double A223 = 0.f; double A233 = 0.f;
			double A313 = 0.f; double A323 = 0.f; double A333 = 0.f;

			double b11 = 0.f; double b12 = 0.f; double b13 = 0.f;
			double b21 = 0.f; double b22 = 0.f; double b23 = 0.f;
			double b31 = 0.f; double b32 = 0.f; double b33 = 0.f;

			// 2nd order tensor in 3D -- 3x3 = 9 components
			double B11 = 0.f; double B12 = 0.f; double B13 = 0.f;
			double B21 = 0.f; double B22 = 0.f; double B23 = 0.f;
			double B31 = 0.f; double B32 = 0.f; double B33 = 0.f;

			// 3rd order tensor in 3D -- 3x3x3 = 27 components
			double C111 = 0.f; double C121 = 0.f; double C131 = 0.f;
			double C211 = 0.f; double C221 = 0.f; double C231 = 0.f;
			double C311 = 0.f; double C321 = 0.f; double C331 = 0.f;
			double C112 = 0.f; double C122 = 0.f; double C132 = 0.f;
			double C212 = 0.f; double C222 = 0.f; double C232 = 0.f;
			double C312 = 0.f; double C322 = 0.f; double C332 = 0.f;
			double C113 = 0.f; double C123 = 0.f; double C133 = 0.f;
			double C213 = 0.f; double C223 = 0.f; double C233 = 0.f;
			double C313 = 0.f; double C323 = 0.f; double C333 = 0.f;

			// 4th order tensor in 3D -- 3x3x3x3 = 81 components
			double D1111 = 0.f; double D1211 = 0.f; double D1311 = 0.f;
			double D2111 = 0.f; double D2211 = 0.f; double D2311 = 0.f;
			double D3111 = 0.f; double D3211 = 0.f; double D3311 = 0.f;
			double D1121 = 0.f; double D1221 = 0.f; double D1321 = 0.f;
			double D2121 = 0.f; double D2221 = 0.f; double D2321 = 0.f;
			double D3121 = 0.f; double D3221 = 0.f; double D3321 = 0.f;
			double D1131 = 0.f; double D1231 = 0.f; double D1331 = 0.f;
			double D2131 = 0.f; double D2231 = 0.f; double D2331 = 0.f;
			double D3131 = 0.f; double D3231 = 0.f; double D3331 = 0.f;
			double D1112 = 0.f; double D1212 = 0.f; double D1312 = 0.f;
			double D2112 = 0.f; double D2212 = 0.f; double D2312 = 0.f;
			double D3112 = 0.f; double D3212 = 0.f; double D3312 = 0.f;
			double D1122 = 0.f; double D1222 = 0.f; double D1322 = 0.f;
			double D2122 = 0.f; double D2222 = 0.f; double D2322 = 0.f;
			double D3122 = 0.f; double D3222 = 0.f; double D3322 = 0.f;
			double D1132 = 0.f; double D1232 = 0.f; double D1332 = 0.f;
			double D2132 = 0.f; double D2232 = 0.f; double D2332 = 0.f;
			double D3132 = 0.f; double D3232 = 0.f; double D3332 = 0.f;
			double D1113 = 0.f; double D1213 = 0.f; double D1313 = 0.f;
			double D2113 = 0.f; double D2213 = 0.f; double D2313 = 0.f;
			double D3113 = 0.f; double D3213 = 0.f; double D3313 = 0.f;
			double D1123 = 0.f; double D1223 = 0.f; double D1323 = 0.f;
			double D2123 = 0.f; double D2223 = 0.f; double D2323 = 0.f;
			double D3123 = 0.f; double D3223 = 0.f; double D3323 = 0.f;
			double D1133 = 0.f; double D1233 = 0.f; double D1333 = 0.f;
			double D2133 = 0.f; double D2233 = 0.f; double D2333 = 0.f;
			double D3133 = 0.f; double D3233 = 0.f; double D3333 = 0.f;

			// 4th order tensor in 3D -- 3x3x3x3 = 81 components
			double E1111 = 0.f; double E1211 = 0.f; double E1311 = 0.f;
			double E2111 = 0.f; double E2211 = 0.f; double E2311 = 0.f;
			double E3111 = 0.f; double E3211 = 0.f; double E3311 = 0.f;
			double E1121 = 0.f; double E1221 = 0.f; double E1321 = 0.f;
			double E2121 = 0.f; double E2221 = 0.f; double E2321 = 0.f;
			double E3121 = 0.f; double E3221 = 0.f; double E3321 = 0.f;
			double E1131 = 0.f; double E1231 = 0.f; double E1331 = 0.;
			double E2131 = 0.f; double E2231 = 0.f; double E2331 = 0.f;
			double E3131 = 0.f; double E3231 = 0.f; double E3331 = 0.f;
			double E1112 = 0.f; double E1212 = 0.f; double E1312 = 0.f;
			double E2112 = 0.f; double E2212 = 0.f; double E2312 = 0.f;
			double E3112 = 0.f; double E3212 = 0.f; double E3312 = 0.f;
			double E1122 = 0.f; double E1222 = 0.f; double E1322 = 0.f;
			double E2122 = 0.f; double E2222 = 0.f; double E2322 = 0.f;
			double E3122 = 0.f; double E3222 = 0.f; double E3322 = 0.f;
			double E1132 = 0.f; double E1232 = 0.f; double E1332 = 0.f;
			double E2132 = 0.f; double E2232 = 0.f; double E2332 = 0.f;
			double E3132 = 0.f; double E3232 = 0.f; double E3332 = 0.f;
			double E1113 = 0.f; double E1213 = 0.f; double E1313 = 0.f;
			double E2113 = 0.f; double E2213 = 0.f; double E2313 = 0.f;
			double E3113 = 0.f; double E3213 = 0.f; double E3313 = 0.f;
			double E1123 = 0.f; double E1223 = 0.f; double E1323 = 0.f;
			double E2123 = 0.f; double E2223 = 0.f; double E2323 = 0.f;
			double E3123 = 0.f; double E3223 = 0.f; double E3323 = 0.f;
			double E1133 = 0.f; double E1233 = 0.f; double E1333 = 0.f;
			double E2133 = 0.f; double E2233 = 0.f; double E2333 = 0.f;
			double E3133 = 0.f; double E3233 = 0.f; double E3333 = 0.f;

			// 4th order tensor in 3D -- 3x3x3x3 = 81 components
			double F1111 = 0.f; double F1211 = 0.f; double F1311 = 0.f;
			double F2111 = 0.f; double F2211 = 0.f; double F2311 = 0.f;
			double F3111 = 0.f; double F3211 = 0.f; double F3311 = 0.f;
			double F1121 = 0.f; double F1221 = 0.f; double F1321 = 0.f;
			double F2121 = 0.f; double F2221 = 0.f; double F2321 = 0.f;
			double F3121 = 0.f; double F3221 = 0.f; double F3321 = 0.f;
			double F1131 = 0.f; double F1231 = 0.f; double F1331 = 0.;
			double F2131 = 0.f; double F2231 = 0.f; double F2331 = 0.f;
			double F3131 = 0.f; double F3231 = 0.f; double F3331 = 0.f;
			double F1112 = 0.f; double F1212 = 0.f; double F1312 = 0.f;
			double F2112 = 0.f; double F2212 = 0.f; double F2312 = 0.f;
			double F3112 = 0.f; double F3212 = 0.f; double F3312 = 0.f;
			double F1122 = 0.f; double F1222 = 0.f; double F1322 = 0.f;
			double F2122 = 0.f; double F2222 = 0.f; double F2322 = 0.f;
			double F3122 = 0.f; double F3222 = 0.f; double F3322 = 0.f;
			double F1132 = 0.f; double F1232 = 0.f; double F1332 = 0.f;
			double F2132 = 0.f; double F2232 = 0.f; double F2332 = 0.f;
			double F3132 = 0.f; double F3232 = 0.f; double F3332 = 0.f;
			double F1113 = 0.f; double F1213 = 0.f; double F1313 = 0.f;
			double F2113 = 0.f; double F2213 = 0.f; double F2313 = 0.f;
			double F3113 = 0.f; double F3213 = 0.f; double F3313 = 0.f;
			double F1123 = 0.f; double F1223 = 0.f; double F1323 = 0.f;
			double F2123 = 0.f; double F2223 = 0.f; double F2323 = 0.f;
			double F3123 = 0.f; double F3223 = 0.f; double F3323 = 0.f;
			double F1133 = 0.f; double F1233 = 0.f; double F1333 = 0.f;
			double F2133 = 0.f; double F2233 = 0.f; double F2333 = 0.f;
			double F3133 = 0.f; double F3233 = 0.f; double F3333 = 0.f;


			//loop over neighbors
			for (unsigned int ni = low_i; ni < high_i; ni++) {
				for (unsigned int nj = low_j; nj < high_j; nj++) {
					for (unsigned int nk = low_k; nk < high_k; nk++) {

						unsigned int ID = ni*grid_geom->ny*grid_geom->nz + nj*grid_geom->nz + nk;

						for (unsigned int j = grid_geom->cell_start[ID]; j < grid_geom->cell_end[ID]; j++) {

							double xj 	= particles[j].px;
							double yj 	= particles[j].py;
							double zj   = particles[j].pz;
							double rhoj = particles[j].rho;
							double m 	= particles[j].m;

							if (i == j) continue;

							double rx   = xi-xj;
							double ry   = yi-yj;
							double rz   = zi-zj;

							double rij  = sqrt(rx*rx + ry*ry + rz*rz);
							double ex   = rx/rij;
							double ey   = ry/rij;
							double ez   = rz/rij;

							kernel_result w = _kernel(xi,yi,zi,xj,yj,zj,hi);

							// A: left side of B (3rd order tensor in 3D -> 3x3x3)
							A111 += ex * ex * w.w_x * m/rhoj;
							A121 += ex * ey * w.w_x * m/rhoj;
							A131 += ex * ez * w.w_x * m/rhoj;
							A211 += ey * ex * w.w_x * m/rhoj;
							A221 += ey * ey * w.w_x * m/rhoj;
							A231 += ey * ez * w.w_x * m/rhoj;
							A311 += ez * ex * w.w_x * m/rhoj;
							A321 += ez * ey * w.w_x * m/rhoj;
							A331 += ez * ez * w.w_x * m/rhoj;

							A112 += ex * ex * w.w_y * m/rhoj;
							A122 += ex * ey * w.w_y * m/rhoj;
							A132 += ex * ez * w.w_y * m/rhoj;
							A212 += ey * ex * w.w_y * m/rhoj;
							A222 += ey * ey * w.w_y * m/rhoj;
							A232 += ey * ez * w.w_y * m/rhoj;
							A312 += ez * ex * w.w_y * m/rhoj;
							A322 += ez * ey * w.w_y * m/rhoj;
							A332 += ez * ez * w.w_y * m/rhoj;

							A113 += ex * ex * w.w_z * m/rhoj;
							A123 += ex * ey * w.w_z * m/rhoj;
							A133 += ex * ez * w.w_z * m/rhoj;
							A213 += ey * ex * w.w_z * m/rhoj;
							A223 += ey * ey * w.w_z * m/rhoj;
							A233 += ey * ez * w.w_z * m/rhoj;
							A313 += ez * ex * w.w_z * m/rhoj;
							A323 += ez * ey * w.w_z * m/rhoj;
							A333 += ez * ez * w.w_z * m/rhoj;

							// B: renormalization tensor (2nd order tensor in 3D -> 3x3)
							b11 += rx * w.w_x * m/rhoj;
							b12 += ry * w.w_x * m/rhoj;
							b13 += rz * w.w_x * m/rhoj;
							b21 += rx * w.w_y * m/rhoj;
							b22 += ry * w.w_y * m/rhoj;
							b23 += rz * w.w_y * m/rhoj;
							b31 += rx * w.w_z * m/rhoj;
							b32 += ry * w.w_z * m/rhoj;
							b33 += rz * w.w_z * m/rhoj;

							// C: right side of B (3rd order tensor in 3D -> 3x3x3)
							C111 += rx * rx * w.w_x * m/rhoj;
							C121 += rx * ry * w.w_x * m/rhoj;
							C131 += rx * rz * w.w_x * m/rhoj;
							C211 += ry * rx * w.w_x * m/rhoj;
							C221 += ry * ry * w.w_x * m/rhoj;
							C231 += ry * rz * w.w_x * m/rhoj;
							C311 += rz * rx * w.w_x * m/rhoj;
							C321 += rz * ry * w.w_x * m/rhoj;
							C331 += rz * rz * w.w_x * m/rhoj;

							C112 += rx * rx * w.w_y * m/rhoj;
							C122 += rx * ry * w.w_y * m/rhoj;
							C132 += rx * rz * w.w_y * m/rhoj;
							C212 += ry * rx * w.w_y * m/rhoj;
							C222 += ry * ry * w.w_y * m/rhoj;
							C232 += ry * rz * w.w_y * m/rhoj;
							C312 += rz * rx * w.w_y * m/rhoj;
							C322 += rz * ry * w.w_y * m/rhoj;
							C332 += rz * rz * w.w_y * m/rhoj;

							C113 += rx * rx * w.w_z * m/rhoj;
							C123 += rx * ry * w.w_z * m/rhoj;
							C133 += rx * rz * w.w_z * m/rhoj;
							C213 += ry * rx * w.w_z * m/rhoj;
							C223 += ry * ry * w.w_z * m/rhoj;
							C233 += ry * rz * w.w_z * m/rhoj;
							C313 += rz * rx * w.w_z * m/rhoj;
							C323 += rz * ry * w.w_z * m/rhoj;
							C333 += rz * rz * w.w_z * m/rhoj;

							// E: very left term (4th order tensor in 3D -> 3x3x3x3)
							E1111 += rx * ex * ex * w.w_x * m/rhoj;
							E1211 += rx * ex * ey * w.w_x * m/rhoj;
							E1311 += rx * ex * ez * w.w_x * m/rhoj;
							E2111 += rx * ey * ex * w.w_x * m/rhoj;
							E2211 += rx * ey * ey * w.w_x * m/rhoj;
							E2311 += rx * ey * ez * w.w_x * m/rhoj;
							E3111 += rx * ez * ex * w.w_x * m/rhoj;
							E3211 += rx * ez * ey * w.w_x * m/rhoj;
							E3311 += rx * ez * ez * w.w_x * m/rhoj;
							E1121 += ry * ex * ex * w.w_x * m/rhoj;
							E1221 += ry * ex * ey * w.w_x * m/rhoj;
							E1321 += ry * ex * ez * w.w_x * m/rhoj;
							E2121 += ry * ey * ex * w.w_x * m/rhoj;
							E2221 += ry * ey * ey * w.w_x * m/rhoj;
							E2321 += ry * ey * ez * w.w_x * m/rhoj;
							E3121 += ry * ez * ex * w.w_x * m/rhoj;
							E3221 += ry * ez * ey * w.w_x * m/rhoj;
							E3321 += ry * ez * ez * w.w_x * m/rhoj;
							E1131 += rz * ex * ex * w.w_x * m/rhoj;
							E1231 += rz * ex * ey * w.w_x * m/rhoj;
							E1331 += rz * ex * ez * w.w_x * m/rhoj;
							E2131 += rz * ey * ex * w.w_x * m/rhoj;
							E2231 += rz * ey * ey * w.w_x * m/rhoj;
							E2331 += rz * ey * ez * w.w_x * m/rhoj;
							E3131 += rz * ez * ex * w.w_x * m/rhoj;
							E3231 += rz * ez * ey * w.w_x * m/rhoj;
							E3331 += rz * ez * ez * w.w_x * m/rhoj;

							E1112 += rx * ex * ex * w.w_y * m/rhoj;
							E1212 += rx * ex * ey * w.w_y * m/rhoj;
							E1312 += rx * ex * ez * w.w_y * m/rhoj;
							E2112 += rx * ey * ex * w.w_y * m/rhoj;
							E2212 += rx * ey * ey * w.w_y * m/rhoj;
							E2312 += rx * ey * ez * w.w_y * m/rhoj;
							E3112 += rx * ez * ex * w.w_y * m/rhoj;
							E3212 += rx * ez * ey * w.w_y * m/rhoj;
							E3312 += rx * ez * ez * w.w_y * m/rhoj;
							E1122 += ry * ex * ex * w.w_y * m/rhoj;
							E1222 += ry * ex * ey * w.w_y * m/rhoj;
							E1322 += ry * ex * ez * w.w_y * m/rhoj;
							E2122 += ry * ey * ex * w.w_y * m/rhoj;
							E2222 += ry * ey * ey * w.w_y * m/rhoj;
							E2322 += ry * ey * ez * w.w_y * m/rhoj;
							E3122 += ry * ez * ex * w.w_y * m/rhoj;
							E3222 += ry * ez * ey * w.w_y * m/rhoj;
							E3322 += ry * ez * ez * w.w_y * m/rhoj;
							E1132 += rz * ex * ex * w.w_y * m/rhoj;
							E1232 += rz * ex * ey * w.w_y * m/rhoj;
							E1332 += rz * ex * ez * w.w_y * m/rhoj;
							E2132 += rz * ey * ex * w.w_y * m/rhoj;
							E2232 += rz * ey * ey * w.w_y * m/rhoj;
							E2332 += rz * ey * ez * w.w_y * m/rhoj;
							E3132 += rz * ez * ex * w.w_y * m/rhoj;
							E3232 += rz * ez * ey * w.w_y * m/rhoj;
							E3332 += rz * ez * ez * w.w_y * m/rhoj;

							E1113 += rx * ex * ex * w.w_z * m/rhoj;
							E1213 += rx * ex * ey * w.w_z * m/rhoj;
							E1313 += rx * ex * ez * w.w_z * m/rhoj;
							E2113 += rx * ey * ex * w.w_z * m/rhoj;
							E2213 += rx * ey * ey * w.w_z * m/rhoj;
							E2313 += rx * ey * ez * w.w_z * m/rhoj;
							E3113 += rx * ez * ex * w.w_z * m/rhoj;
							E3213 += rx * ez * ey * w.w_z * m/rhoj;
							E3313 += rx * ez * ez * w.w_z * m/rhoj;
							E1123 += ry * ex * ex * w.w_z * m/rhoj;
							E1223 += ry * ex * ey * w.w_z * m/rhoj;
							E1323 += ry * ex * ez * w.w_z * m/rhoj;
							E2123 += ry * ey * ex * w.w_z * m/rhoj;
							E2223 += ry * ey * ey * w.w_z * m/rhoj;
							E2323 += ry * ey * ez * w.w_z * m/rhoj;
							E3123 += ry * ez * ex * w.w_z * m/rhoj;
							E3223 += ry * ez * ey * w.w_z * m/rhoj;
							E3323 += ry * ez * ez * w.w_z * m/rhoj;
							E1133 += rz * ex * ex * w.w_z * m/rhoj;
							E1233 += rz * ex * ey * w.w_z * m/rhoj;
							E1333 += rz * ex * ez * w.w_z * m/rhoj;
							E2133 += rz * ey * ex * w.w_z * m/rhoj;
							E2233 += rz * ey * ey * w.w_z * m/rhoj;
							E2333 += rz * ey * ez * w.w_z * m/rhoj;
							E3133 += rz * ez * ex * w.w_z * m/rhoj;
							E3233 += rz * ez * ey * w.w_z * m/rhoj;
							E3333 += rz * ez * ez * w.w_z * m/rhoj;
						}
					}
				}
			}

			double binv[9];
			double bMAT[9] = {b11,b12,b13,
					  		  b21,b22,b23,
							  b31,b32,b33};

			// B: re-normalization tensor (i.e. Eq(60) from Fatehi's paper)
			inv3x3(bMAT,binv);

			B11 = -binv[0]; B12 = -binv[1]; B13 = -binv[2];
			B21 = -binv[3]; B22 = -binv[4]; B23 = -binv[5];
			B31 = -binv[6]; B32 = -binv[7]; B33 = -binv[8];

			// D = A.B.C (tensor multiplication ---> 4th order tensor in 3D -----> index: ijkl 3x3x3x3)
			D1111 = A111*B11*C111 + A111*B12*C211 + A111*B13*C311 +
					A112*B21*C111 + A112*B22*C211 + A112*B23*C311 +
					A113*B31*C111 + A113*B32*C211 + A113*B33*C311;
			D1211 = A121*B11*C111 + A121*B12*C211 + A121*B13*C311 +
					A122*B21*C111 + A122*B22*C211 + A122*B23*C311 +
					A123*B31*C111 + A123*B32*C211 + A123*B33*C311;
			D1311 = A131*B11*C111 + A131*B12*C211 + A131*B13*C311 +
					A132*B21*C111 + A132*B22*C211 + A132*B23*C311 +
					A133*B31*C111 + A133*B32*C211 + A133*B33*C311;
			D2111 = A211*B11*C111 + A211*B12*C211 + A211*B13*C311 +
					A212*B21*C111 + A212*B22*C211 + A212*B23*C311 +
					A213*B31*C111 + A213*B32*C211 + A213*B33*C311;
			D2211 = A221*B11*C111 + A221*B12*C211 + A221*B13*C311 +
					A222*B21*C111 + A222*B22*C211 + A222*B23*C311 +
					A223*B31*C111 + A223*B32*C211 + A223*B33*C311;
			D2311 = A231*B11*C111 + A231*B12*C211 + A231*B13*C311 +
					A232*B21*C111 + A232*B22*C211 + A232*B23*C311 +
					A233*B31*C111 + A233*B32*C211 + A233*B33*C311;
			D3111 = A311*B11*C111 + A311*B12*C211 + A311*B13*C311 +
					A312*B21*C111 + A312*B22*C211 + A312*B23*C311 +
					A313*B31*C111 + A313*B32*C211 + A313*B33*C311;
			D3211 = A321*B11*C111 + A321*B12*C211 + A321*B13*C311 +
					A322*B21*C111 + A322*B22*C211 + A322*B23*C311 +
					A323*B31*C111 + A323*B32*C211 + A323*B33*C311;
			D3311 = A331*B11*C111 + A331*B12*C211 + A331*B13*C311 +
					A332*B21*C111 + A332*B22*C211 + A332*B23*C311 +
					A333*B31*C111 + A333*B32*C211 + A333*B33*C311;
			D1121 = A111*B11*C121 + A111*B12*C221 + A111*B13*C321 +
					A112*B21*C121 + A112*B22*C221 + A112*B23*C321 +
					A113*B31*C121 + A113*B32*C221 + A113*B33*C321;
			D1221 = A121*B11*C121 + A121*B12*C221 + A121*B13*C321 +
					A122*B21*C121 + A122*B22*C221 + A122*B23*C321 +
					A123*B31*C121 + A123*B32*C221 + A123*B33*C321;
			D1321 = A131*B11*C121 + A131*B12*C221 + A131*B13*C321 +
					A132*B21*C121 + A132*B22*C221 + A132*B23*C321 +
					A133*B31*C121 + A133*B32*C221 + A133*B33*C321;
			D2121 = A211*B11*C121 + A211*B12*C221 + A211*B13*C321 +
					A212*B21*C121 + A212*B22*C221 + A212*B23*C321 +
					A213*B31*C121 + A213*B32*C221 + A213*B33*C321;
			D2221 = A221*B11*C121 + A221*B12*C221 + A221*B13*C321 +
					A222*B21*C121 + A222*B22*C221 + A222*B23*C321 +
					A223*B31*C121 + A223*B32*C221 + A223*B33*C321;
			D2321 = A231*B11*C121 + A231*B12*C221 + A231*B13*C321 +
					A232*B21*C121 + A232*B22*C221 + A232*B23*C321 +
					A233*B31*C121 + A233*B32*C221 + A233*B33*C321;
			D3121 = A311*B11*C121 + A311*B12*C221 + A311*B13*C321 +
					A312*B21*C121 + A312*B22*C221 + A312*B23*C321 +
					A313*B31*C121 + A313*B32*C221 + A313*B33*C321;
			D3221 = A321*B11*C121 + A321*B12*C221 + A321*B13*C321 +
					A322*B21*C121 + A322*B22*C221 + A322*B23*C321 +
					A323*B31*C121 + A323*B32*C221 + A323*B33*C321;
			D3321 = A331*B11*C121 + A331*B12*C221 + A331*B13*C321 +
					A332*B21*C121 + A332*B22*C221 + A332*B23*C321 +
					A333*B31*C121 + A333*B32*C221 + A333*B33*C321;
			D1131 = A111*B11*C131 + A111*B12*C231 + A111*B13*C331 +
					A112*B21*C131 + A112*B22*C231 + A112*B23*C331 +
					A113*B31*C131 + A113*B32*C231 + A113*B33*C331;
			D1231 = A121*B11*C131 + A121*B12*C231 + A121*B13*C331 +
					A122*B21*C131 + A122*B22*C231 + A122*B23*C331 +
					A123*B31*C131 + A123*B32*C231 + A123*B33*C331;
			D1331 = A131*B11*C131 + A131*B12*C231 + A131*B13*C331 +
					A132*B21*C131 + A132*B22*C231 + A132*B23*C331 +
					A133*B31*C131 + A133*B32*C231 + A133*B33*C331;
			D2131 = A211*B11*C131 + A211*B12*C231 + A211*B13*C331 +
					A212*B21*C131 + A212*B22*C231 + A212*B23*C331 +
					A213*B31*C131 + A213*B32*C231 + A213*B33*C331;
			D2231 = A221*B11*C131 + A221*B12*C231 + A221*B13*C331 +
					A222*B21*C131 + A222*B22*C231 + A222*B23*C331 +
					A223*B31*C131 + A223*B32*C231 + A223*B33*C331;
			D2331 = A231*B11*C131 + A231*B12*C231 + A231*B13*C331 +
					A232*B21*C131 + A232*B22*C231 + A232*B23*C331 +
					A233*B31*C131 + A233*B32*C231 + A233*B33*C331;
			D3131 = A311*B11*C131 + A311*B12*C231 + A311*B13*C331 +
					A312*B21*C131 + A312*B22*C231 + A312*B23*C331 +
					A313*B31*C131 + A313*B32*C231 + A313*B33*C331;
			D3231 = A321*B11*C131 + A321*B12*C231 + A321*B13*C331 +
					A322*B21*C131 + A322*B22*C231 + A322*B23*C331 +
					A323*B31*C131 + A323*B32*C231 + A323*B33*C331;
			D3331 = A331*B11*C131 + A331*B12*C231 + A331*B13*C331 +
					A332*B21*C131 + A332*B22*C231 + A332*B23*C331 +
					A333*B31*C131 + A333*B32*C231 + A333*B33*C331;

			D1112 = A111*B11*C112 + A111*B12*C212 + A111*B13*C312 +
					A112*B21*C112 + A112*B22*C212 + A112*B23*C312 +
					A113*B31*C112 + A113*B32*C212 + A113*B33*C312;
			D1212 = A121*B11*C112 + A121*B12*C212 + A121*B13*C312 +
					A122*B21*C112 + A122*B22*C212 + A122*B23*C312 +
					A123*B31*C112 + A123*B32*C212 + A123*B33*C312;
			D1312 = A131*B11*C112 + A131*B12*C212 + A131*B13*C312 +
					A132*B21*C112 + A132*B22*C212 + A132*B23*C312 +
					A133*B31*C112 + A133*B32*C212 + A133*B33*C312;
			D2112 = A211*B11*C112 + A211*B12*C212 + A211*B13*C312 +
					A212*B21*C112 + A212*B22*C212 + A212*B23*C312 +
					A213*B31*C112 + A213*B32*C212 + A213*B33*C312;
			D2212 = A221*B11*C112 + A221*B12*C212 + A221*B13*C312 +
					A222*B21*C112 + A222*B22*C212 + A222*B23*C312 +
					A223*B31*C112 + A223*B32*C212 + A223*B33*C312;
			D2312 = A231*B11*C112 + A231*B12*C212 + A231*B13*C312 +
					A232*B21*C112 + A232*B22*C212 + A232*B23*C312 +
					A233*B31*C112 + A233*B32*C212 + A233*B33*C312;
			D3112 = A311*B11*C112 + A311*B12*C212 + A311*B13*C312 +
					A312*B21*C112 + A312*B22*C212 + A312*B23*C312 +
					A313*B31*C112 + A313*B32*C212 + A313*B33*C312;
			D3212 = A321*B11*C112 + A321*B12*C212 + A321*B13*C312 +
					A322*B21*C112 + A322*B22*C212 + A322*B23*C312 +
					A323*B31*C112 + A323*B32*C212 + A323*B33*C312;
			D3312 = A331*B11*C112 + A331*B12*C212 + A331*B13*C312 +
					A332*B21*C112 + A332*B22*C212 + A332*B23*C312 +
					A333*B31*C112 + A333*B32*C212 + A333*B33*C312;
			D1122 = A111*B11*C122 + A111*B12*C222 + A111*B13*C322 +
					A112*B21*C122 + A112*B22*C222 + A112*B23*C322 +
					A113*B31*C122 + A113*B32*C222 + A113*B33*C322;
			D1222 = A121*B11*C122 + A121*B12*C222 + A121*B13*C322 +
					A122*B21*C122 + A122*B22*C222 + A122*B23*C322 +
					A123*B31*C122 + A123*B32*C222 + A123*B33*C322;
			D1322 = A131*B11*C122 + A131*B12*C222 + A131*B13*C322 +
					A132*B21*C122 + A132*B22*C222 + A132*B23*C322 +
					A133*B31*C122 + A133*B32*C222 + A133*B33*C322;
			D2122 = A211*B11*C122 + A211*B12*C222 + A211*B13*C322 +
					A212*B21*C122 + A212*B22*C222 + A212*B23*C322 +
					A213*B31*C122 + A213*B32*C222 + A213*B33*C322;
			D2222 = A221*B11*C122 + A221*B12*C222 + A221*B13*C322 +
					A222*B21*C122 + A222*B22*C222 + A222*B23*C322 +
					A223*B31*C122 + A223*B32*C222 + A223*B33*C322;
			D2322 = A231*B11*C122 + A231*B12*C222 + A231*B13*C322 +
					A232*B21*C122 + A232*B22*C222 + A232*B23*C322 +
					A233*B31*C122 + A233*B32*C222 + A233*B33*C322;
			D3122 = A311*B11*C122 + A311*B12*C222 + A311*B13*C322 +
					A312*B21*C122 + A312*B22*C222 + A312*B23*C322 +
					A313*B31*C122 + A313*B32*C222 + A313*B33*C322;
			D3222 = A321*B11*C122 + A321*B12*C222 + A321*B13*C322 +
					A322*B21*C122 + A322*B22*C222 + A322*B23*C322 +
					A323*B31*C122 + A323*B32*C222 + A323*B33*C322;
			D3322 = A331*B11*C122 + A331*B12*C222 + A331*B13*C322 +
					A332*B21*C122 + A332*B22*C222 + A332*B23*C322 +
					A333*B31*C122 + A333*B32*C222 + A333*B33*C322;
			D1132 = A111*B11*C132 + A111*B12*C232 + A111*B13*C332 +
					A112*B21*C132 + A112*B22*C232 + A112*B23*C332 +
					A113*B31*C132 + A113*B32*C232 + A113*B33*C332;
			D1232 = A121*B11*C132 + A121*B12*C232 + A121*B13*C332 +
					A122*B21*C132 + A122*B22*C232 + A122*B23*C332 +
					A123*B31*C132 + A123*B32*C232 + A123*B33*C332;
			D1332 = A131*B11*C132 + A131*B12*C232 + A131*B13*C332 +
					A132*B21*C132 + A132*B22*C232 + A132*B23*C332 +
					A133*B31*C132 + A133*B32*C232 + A133*B33*C332;
			D2132 = A211*B11*C132 + A211*B12*C232 + A211*B13*C332 +
					A212*B21*C132 + A212*B22*C232 + A212*B23*C332 +
					A213*B31*C132 + A213*B32*C232 + A213*B33*C332;
			D2232 = A221*B11*C132 + A221*B12*C232 + A221*B13*C332 +
					A222*B21*C132 + A222*B22*C232 + A222*B23*C332 +
					A223*B31*C132 + A223*B32*C232 + A223*B33*C332;
			D2332 = A231*B11*C132 + A231*B12*C232 + A231*B13*C332 +
					A232*B21*C132 + A232*B22*C232 + A232*B23*C332 +
					A233*B31*C132 + A233*B32*C232 + A233*B33*C332;
			D3132 = A311*B11*C132 + A311*B12*C232 + A311*B13*C332 +
					A312*B21*C132 + A312*B22*C232 + A312*B23*C332 +
					A313*B31*C132 + A313*B32*C232 + A313*B33*C332;
			D3232 = A321*B11*C132 + A321*B12*C232 + A321*B13*C332 +
					A322*B21*C132 + A322*B22*C232 + A322*B23*C332 +
					A323*B31*C132 + A323*B32*C232 + A323*B33*C332;
			D3332 = A331*B11*C132 + A331*B12*C232 + A331*B13*C332 +
					A332*B21*C132 + A332*B22*C232 + A332*B23*C332 +
					A333*B31*C132 + A333*B32*C232 + A333*B33*C332;

			D1113 = A111*B11*C113 + A111*B12*C213 + A111*B13*C313 +
					A112*B21*C113 + A112*B22*C213 + A112*B23*C313 +
					A113*B31*C113 + A113*B32*C213 + A113*B33*C313;
			D1213 = A121*B11*C113 + A121*B12*C213 + A121*B13*C313 +
					A122*B21*C113 + A122*B22*C213 + A122*B23*C313 +
					A123*B31*C113 + A123*B32*C213 + A123*B33*C313;
			D1313 = A131*B11*C113 + A131*B12*C213 + A131*B13*C313 +
					A132*B21*C113 + A132*B22*C213 + A132*B23*C313 +
					A133*B31*C113 + A133*B32*C213 + A133*B33*C313;
			D2113 = A211*B11*C113 + A211*B12*C213 + A211*B13*C313 +
					A212*B21*C113 + A212*B22*C213 + A212*B23*C313 +
					A213*B31*C113 + A213*B32*C213 + A213*B33*C313;
			D2213 = A221*B11*C113 + A221*B12*C213 + A221*B13*C313 +
					A222*B21*C113 + A222*B22*C213 + A222*B23*C313 +
					A223*B31*C113 + A223*B32*C213 + A223*B33*C313;
			D2313 = A231*B11*C113 + A231*B12*C213 + A231*B13*C313 +
					A232*B21*C113 + A232*B22*C213 + A232*B23*C313 +
					A233*B31*C113 + A233*B32*C213 + A233*B33*C313;
			D3113 = A311*B11*C113 + A311*B12*C213 + A311*B13*C313 +
					A312*B21*C113 + A312*B22*C213 + A312*B23*C313 +
					A313*B31*C113 + A313*B32*C213 + A313*B33*C313;
			D3213 = A321*B11*C113 + A321*B12*C213 + A321*B13*C313 +
					A322*B21*C113 + A322*B22*C213 + A322*B23*C313 +
					A323*B31*C113 + A323*B32*C213 + A323*B33*C313;
			D3313 = A331*B11*C113 + A331*B12*C213 + A331*B13*C313 +
					A332*B21*C113 + A332*B22*C213 + A332*B23*C313 +
					A333*B31*C113 + A333*B32*C213 + A333*B33*C313;
			D1123 = A111*B11*C123 + A111*B12*C223 + A111*B13*C323 +
					A112*B21*C123 + A112*B22*C223 + A112*B23*C323 +
					A113*B31*C123 + A113*B32*C223 + A113*B33*C323;
			D1223 = A121*B11*C123 + A121*B12*C223 + A121*B13*C323 +
					A122*B21*C123 + A122*B22*C223 + A122*B23*C323 +
					A123*B31*C123 + A123*B32*C223 + A123*B33*C323;
			D1323 = A131*B11*C123 + A131*B12*C223 + A131*B13*C323 +
					A132*B21*C123 + A132*B22*C223 + A132*B23*C323 +
					A133*B31*C123 + A133*B32*C223 + A133*B33*C323;
			D2123 = A211*B11*C123 + A211*B12*C223 + A211*B13*C323 +
					A212*B21*C123 + A212*B22*C223 + A212*B23*C323 +
					A213*B31*C123 + A213*B32*C223 + A213*B33*C323;
			D2223 = A221*B11*C123 + A221*B12*C223 + A221*B13*C323 +
					A222*B21*C123 + A222*B22*C223 + A222*B23*C323 +
					A223*B31*C123 + A223*B32*C223 + A223*B33*C323;
			D2323 = A231*B11*C123 + A231*B12*C223 + A231*B13*C323 +
					A232*B21*C123 + A232*B22*C223 + A232*B23*C323 +
					A233*B31*C123 + A233*B32*C223 + A233*B33*C323;
			D3123 = A311*B11*C123 + A311*B12*C223 + A311*B13*C323 +
					A312*B21*C123 + A312*B22*C223 + A312*B23*C323 +
					A313*B31*C123 + A313*B32*C223 + A313*B33*C323;
			D3223 = A321*B11*C123 + A321*B12*C223 + A321*B13*C323 +
					A322*B21*C123 + A322*B22*C223 + A322*B23*C323 +
					A323*B31*C123 + A323*B32*C223 + A323*B33*C323;
			D3323 = A331*B11*C123 + A331*B12*C223 + A331*B13*C323 +
					A332*B21*C123 + A332*B22*C223 + A332*B23*C323 +
					A333*B31*C123 + A333*B32*C223 + A333*B33*C323;
			D1133 = A111*B11*C133 + A111*B12*C233 + A111*B13*C333 +
					A112*B21*C133 + A112*B22*C233 + A112*B23*C333 +
					A113*B31*C133 + A113*B32*C233 + A113*B33*C333;
			D1233 = A121*B11*C133 + A121*B12*C233 + A121*B13*C333 +
					A122*B21*C133 + A122*B22*C233 + A122*B23*C333 +
					A123*B31*C133 + A123*B32*C233 + A123*B33*C333;
			D1333 = A131*B11*C133 + A131*B12*C233 + A131*B13*C333 +
					A132*B21*C133 + A132*B22*C233 + A132*B23*C333 +
					A133*B31*C133 + A133*B32*C233 + A133*B33*C333;
			D2133 = A211*B11*C133 + A211*B12*C233 + A211*B13*C333 +
					A212*B21*C133 + A212*B22*C233 + A212*B23*C333 +
					A213*B31*C133 + A213*B32*C233 + A213*B33*C333;
			D2233 = A221*B11*C133 + A221*B12*C233 + A221*B13*C333 +
					A222*B21*C133 + A222*B22*C233 + A222*B23*C333 +
					A223*B31*C133 + A223*B32*C233 + A223*B33*C333;
			D2333 = A231*B11*C133 + A231*B12*C233 + A231*B13*C333 +
					A232*B21*C133 + A232*B22*C233 + A232*B23*C333 +
					A233*B31*C133 + A233*B32*C233 + A233*B33*C333;
			D3133 = A311*B11*C133 + A311*B12*C233 + A311*B13*C333 +
					A312*B21*C133 + A312*B22*C233 + A312*B23*C333 +
					A313*B31*C133 + A313*B32*C233 + A313*B33*C333;
			D3233 = A321*B11*C133 + A321*B12*C233 + A321*B13*C333 +
					A322*B21*C133 + A322*B22*C233 + A322*B23*C333 +
					A323*B31*C133 + A323*B32*C233 + A323*B33*C333;
			D3333 = A331*B11*C133 + A331*B12*C233 + A331*B13*C333 +
					A332*B21*C133 + A332*B22*C233 + A332*B23*C333 +
					A333*B31*C133 + A333*B32*C233 + A333*B33*C333;

			// F = E + D (tensor addition ---> 4th order tensor in 3D)
			F1111=E1111+D1111; F1211=E1211+D1211; F1311=E1311+D1311;
			F2111=E2111+D2111; F2211=E2211+D2211; F2311=E2311+D2311;
			F3111=E3111+D3111; F3211=E3211+D3211; F3311=E3311+D3311;
			F1121=E1121+D1121; F1221=E1221+D1221; F1321=E1321+D1321;
			F2121=E2121+D2121; F2221=E2221+D2221; F2321=E2321+D2321;
			F3121=E3121+D3121; F3221=E3221+D3221; F3321=E3321+D3321;
			F1131=E1131+D1131; F1231=E1231+D1231; F1331=E1331+D1331;
			F2131=E2131+D2131; F2231=E2231+D2231; F2331=E2331+D2331;
			F3131=E3131+D3131; F3231=E3231+D3231; F3331=E3331+D3331;

			F1112=E1112+D1112; F1212=E1212+D1212; F1312=E1312+D1312;
			F2112=E2112+D2112; F2212=E2212+D2212; F2312=E2312+D2312;
			F3112=E3112+D3112; F3212=E3212+D3212; F3312=E3312+D3312;
			F1122=E1122+D1122; F1222=E1222+D1222; F1322=E1322+D1322;
			F2122=E2122+D2122; F2222=E2222+D2222; F2322=E2322+D2322;
			F3122=E3122+D3122; F3222=E3222+D3222; F3322=E3322+D3322;
			F1132=E1132+D1132; F1232=E1232+D1232; F1332=E1332+D1332;
			F2132=E2132+D2132; F2232=E2232+D2232; F2332=E2332+D2332;
			F3132=E3132+D3132; F3232=E3232+D3232; F3332=E3332+D3332;

			F1113=E1113+D1113; F1213=E1213+D1213; F1313=E1313+D1313;
			F2113=E2113+D2113; F2213=E2213+D2213; F2313=E2313+D2313;
			F3113=E3113+D3113; F3213=E3213+D3213; F3313=E3313+D3313;
			F1123=E1123+D1123; F1223=E1223+D1223; F1323=E1323+D1323;
			F2123=E2123+D2123; F2223=E2223+D2223; F2323=E2323+D2323;
			F3123=E3123+D3123; F3223=E3223+D3223; F3323=E3323+D3323;
			F1133=E1133+D1133; F1233=E1233+D1233; F1333=E1333+D1333;
			F2133=E2133+D2133; F2233=E2233+D2233; F2333=E2333+D2333;
			F3133=E3133+D3133; F3233=E3233+D3233; F3333=E3333+D3333;

			/*
			    B_hat calculation by solving a system of 6-linear equations (i.e. Eq(40) from Fatehi paper)

			    in Eq(41) from Fatehi paper: B_hat:F = -I
			    in which: B_hat is 2nd order tensor, F is 4th order tensor and I is 2nd order unity tensor.
			    due to symmetry, the system of equations can be reduced
			    so the rank is 6 (not 9!)
			 */

			double G[36] = {F1111, 2.*F1211, 2.*F1311, F2211, 2.*F2311, F3311,
							F1112, 2.*F1212, 2.*F1312, F2212, 2.*F2312, F3312,
							F1113, 2.*F1213, 2.*F1313, F2213, 2.*F2313, F3313,
							F1122, 2.*F1222, 2.*F1322, F2222, 2.*F2322, F3322,
							F1123, 2.*F1223, 2.*F1323, F2223, 2.*F2323, F3323,
							F1133, 2.*F1233, 2.*F1333, F2233, 2.*F2333, F3333};

			double rhs[6] = {-1.0, 0.0, 0.0, -1.0, 0.0, -1.0};
			double resBhat[6];
			solve6x6(G,rhs,resBhat);

			double Bhati[9] = {resBhat[0],resBhat[1],resBhat[2],
							   resBhat[1],resBhat[3],resBhat[4],
							   resBhat[2],resBhat[4],resBhat[5]};

			// correction terms (in accord. with CSPM definition)
			double Bi[9] = {-b11,-b12,-b13,
							-b21,-b22,-b23,
							-b31,-b32,-b33};


			// chunks of correction terms for entire particles
			memcpy(B    + 9*i, Bi, sizeof(double)*9);
			memcpy(Bhat + 9*i, Bhati, sizeof(double)*9);
		}
	}
}

void perform_nmfs(particle* particles) {
	singleton_geometry  *grid_geom  = get_singleton_geometry();

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

			double xi   = particles[i].px;
			double yi   = particles[i].py;
			double zi   = particles[i].pz;
			double hi   = particles[i].h;
			double fi   = particles[i].f;

			double nom = 0.f;
			double den = 0.f;
			double F1 = 0.f;
			double F2 = 0.f;
			double F3 = 0.f;

			//loop over neighbors
			for (unsigned int ni = low_i; ni < high_i; ni++) {
				for (unsigned int nj = low_j; nj < high_j; nj++) {
					for (unsigned int nk = low_k; nk < high_k; nk++) {

						unsigned int ID = ni*grid_geom->ny*grid_geom->nz + nj*grid_geom->nz + nk;

						for (unsigned int j = grid_geom->cell_start[ID]; j < grid_geom->cell_end[ID]; j++) {

							double xj   = particles[j].px;
							double yj   = particles[j].py;
							double zj   = particles[j].pz;
							double fj   = particles[j].f;
							double m    = particles[j].m;
							double rhoj = particles[j].rho;

							if (particles[i].label==INLET && particles[j].label!=INLET) fj=fi;
							if (particles[i].label==CONNECTOR_TOP && particles[j].label==INLET) fj=fi;
							if (particles[i].label==CONNECTOR_TOP && particles[j].label==CONNECTOR_INLET) fj=fi;
							if (particles[i].label==CONNECTOR_INLET && particles[j].label==INLET) fj=fi;
							if (particles[i].label==CONNECTOR_INLET && particles[j].label==CONNECTOR_TOP) fj=fi;
							if (particles[i].label==BAR && particles[j].label==CONNECTOR_BAR_CENTER) fj=fi;
							if (particles[i].label==CONNECTOR_BAR_CENTER && particles[j].label==BAR) fj=fi;
							if (particles[j].bnd) fj = fi;

							kernel_result w = _kernel(xi,yi,zi,xj,yj,zj,hi);

							nom +=  fj*w.w*m/rhoj;
							den +=     w.w*m/rhoj;

							F1 += (fj - fi) * w.w_x * m/rhoj;
							F2 += (fj - fi) * w.w_y * m/rhoj;
							F3 += (fj - fi) * w.w_z * m/rhoj;
						}
					}
				}
			}

			// solving 3x3 system of simultaneous linear equations for 1st derivatives
			double Bi[9];
			memcpy(Bi,B + 9*i,sizeof(double)*9);
			double rhs1[3] = {F1, F2, F3};
			double ders1[3];

			solve3x3(Bi,rhs1,ders1);
			double der_x = ders1[0];
			double der_y = ders1[1];
			double der_z = ders1[2];
			particles[i].f_x = der_x;
			particles[i].f_y = der_y;
			particles[i].f_z = der_z;

			// rhs terms for 2nd derivatives & of Laplacian
			double fac11 = 0.f; double fac12 = 0.f; double fac13 = 0.f;
			double fac21 = 0.f; double fac22 = 0.f; double fac23 = 0.f;
			double fac31 = 0.f; double fac32 = 0.f; double fac33 = 0.f;

			//loop over neighbors
			for (unsigned int ni = low_i; ni < high_i; ni++) {
				for (unsigned int nj = low_j; nj < high_j; nj++) {
					for (unsigned int nk = low_k; nk < high_k; nk++) {

						unsigned int ID = ni*grid_geom->ny*grid_geom->nz + nj*grid_geom->nz + nk;

						for (unsigned int j = grid_geom->cell_start[ID]; j < grid_geom->cell_end[ID]; j++) {

							double xj 	= particles[j].px;
							double yj 	= particles[j].py;
							double zj   = particles[j].pz;
							double fj   = particles[j].f;
							double rhoj = particles[j].rho;
							double	m   = particles[j].m;

							if (particles[i].label==INLET && particles[j].label!=INLET) fj=fi;
							if (particles[i].label==CONNECTOR_TOP && particles[j].label==INLET) fj=fi;
							if (particles[i].label==CONNECTOR_TOP && particles[j].label==CONNECTOR_INLET) fj=fi;
							if (particles[i].label==CONNECTOR_INLET && particles[j].label==INLET) fj=fi;
							if (particles[i].label==CONNECTOR_INLET && particles[j].label==CONNECTOR_TOP) fj=fi;
							if (particles[i].label==BAR && particles[j].label==CONNECTOR_BAR_CENTER) fj=fi;
							if (particles[i].label==CONNECTOR_BAR_CENTER && particles[j].label==BAR) fj=fi;
							if (particles[j].bnd) fj = fi;

							if (i == j) continue;

							double rx   = xi-xj;
							double ry   = yi-yj;
							double rz   = zi-zj;

							double rij  = sqrt(rx*rx + ry*ry + rz*rz);
							double rij1 = 1./rij;
							double ex   = rx*rij1;
							double ey   = ry*rij1;
							double ez   = rz*rij1;

							kernel_result w = _kernel(xi,yi,zi,xj,yj,zj,hi);

							// rhs terms (i.e. Eq(41) from Fatehi's paper)
							fac11 += 2.* (m/rhoj) * ex * w.w_x * ((fi - fj)*rij1 - (ex*der_x + ey*der_y + ez*der_z));
							fac12 += 2.* (m/rhoj) * ex * w.w_y * ((fi - fj)*rij1 - (ex*der_x + ey*der_y + ez*der_z));
							fac13 += 2.* (m/rhoj) * ex * w.w_z * ((fi - fj)*rij1 - (ex*der_x + ey*der_y + ez*der_z));

							fac21 += 2.* (m/rhoj) * ey * w.w_x * ((fi - fj)*rij1 - (ex*der_x + ey*der_y + ez*der_z));
							fac22 += 2.* (m/rhoj) * ey * w.w_y * ((fi - fj)*rij1 - (ex*der_x + ey*der_y + ez*der_z));
							fac23 += 2.* (m/rhoj) * ey * w.w_z * ((fi - fj)*rij1 - (ex*der_x + ey*der_y + ez*der_z));

							fac31 += 2.* (m/rhoj) * ez * w.w_x * ((fi - fj)*rij1 - (ex*der_x + ey*der_y + ez*der_z));
							fac32 += 2.* (m/rhoj) * ez * w.w_y * ((fi - fj)*rij1 - (ex*der_x + ey*der_y + ez*der_z));
							fac33 += 2.* (m/rhoj) * ez * w.w_z * ((fi - fj)*rij1 - (ex*der_x + ey*der_y + ez*der_z));
						}
					}
				}
			}

			// unpack the chunk of Bhat
			double B_hat[9];
			memcpy(B_hat,Bhat + 9*i,sizeof(double)*9);

			// Laplacian (i.e. Eq(41) from Fatehi paper) -> double contraction
			particles[i].LaplF = B_hat[0]*fac11 + B_hat[1]*fac12 + B_hat[2]*fac13 +
								 B_hat[3]*fac21 + B_hat[4]*fac22 + B_hat[5]*fac23 +
								 B_hat[6]*fac31 + B_hat[7]*fac32 + B_hat[8]*fac33;
		}
	}
}

void nmfs_init(unsigned int Ntot, kernel kernel_fun) {
	//first derivatives
	B    = (double*) calloc(Ntot*3*3, sizeof(double));
	//second derivatives
	Bhat = (double*) calloc(Ntot*3*3, sizeof(double));

	_kernel = kernel_fun;
}
