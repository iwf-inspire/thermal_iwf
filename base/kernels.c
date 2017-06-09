//Copyright ETH Zurich, IWF

//This file is part of iwf_thermal.

//iwf_thermal is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.

//iwf_thermal is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.

//You should have received a copy of the GNU General Public License
//along with iwf_thermal.  If not, see <http://www.gnu.org/licenses/>.

#include "kernels.h"

kernel_result gaussian3D(double xi, double yi, double zi, double xj, double yj, double zj, double h) {
	double h1   = 1./h;
	double xij  = xi-xj;
	double yij  = yi-yj;
	double zij  = zi-zj;
	double rij  = sqrt(xij*xij + yij*yij + zij*zij);
	double rij3 = rij*rij*rij;

	double fac = 1./((M_PI*h*h)*sqrt(M_PI*h*h));
	double q   = rij*h1;

	kernel_result w;

	w.w    = 0.;
	w.w_x  = 0.;
	w.w_y  = 0.;
	w.w_z  = 0.;
	w.w_xy = 0.;
	w.w_xz = 0.;
	w.w_yz = 0.;
	w.w_xx = 0.;
	w.w_yy = 0.;
	w.w_zz = 0.;

	if (q < 3.) {

		double expq = exp(-q*q);
		double w_q  = -2.*q*expq;
		double w_qq = (4.*q*q-2)*expq;

		w.w = fac*expq;

		if (rij > 1e-8) {
			double q_x  = xij*h1/rij;
			double q_y  = yij*h1/rij;
			double q_z  = zij*h1/rij;
			double q_xy = -xij*yij*h1/rij3;
			double q_xz = -xij*zij*h1/rij3;
			double q_yz = -yij*zij*h1/rij3;
	        double q_xx = (yij*yij + zij*zij)*h1/rij3;
	        double q_yy = (xij*xij + zij*zij)*h1/rij3;
	        double q_zz = (xij*xij + yij*yij)*h1/rij3;

	        w.w_x  = (w_q*q_x*fac);
	        w.w_y  = (w_q*q_y*fac);
	        w.w_z  = (w_q*q_z*fac);
	        w.w_xy = (w_qq*q_x*q_y + w_q*q_xy)*fac;
	        w.w_xz = (w_qq*q_x*q_z + w_q*q_xz)*fac;
	        w.w_yz = (w_qq*q_y*q_z + w_q*q_yz)*fac;
	        w.w_xx = (w_qq*q_x*q_x + w_q*q_xx)*fac;
	        w.w_yy = (w_qq*q_y*q_y + w_q*q_yy)*fac;
	        w.w_zz = (w_qq*q_z*q_z + w_q*q_zz)*fac;
		} else {
			w.w_xx = (w_qq*h1*h1)*fac;
			w.w_yy = (w_qq*h1*h1)*fac;
			w.w_zz = (w_qq*h1*h1)*fac;
		}
	}

	return w;
}

kernel_result cubic_spline3D(double xi, double yi, double zi, double xj, double yj, double zj, double h) {
	double h1   = 1./h;
	double xij  = xi-xj;
	double yij  = yi-yj;
	double zij  = zi-zj;
	double rij  = sqrt(xij*xij + yij*yij + zij*zij);
	double rij3 = rij*rij*rij;

	double fac = (M_1_PI)*h1*h1*h1;
	double q   = rij*h1;

	kernel_result w;

	w.w    = 0.;
	w.w_x  = 0.;
	w.w_y  = 0.;
	w.w_z  = 0.;
	w.w_xy = 0.;
	w.w_xz = 0.;
	w.w_yz = 0.;
	w.w_xx = 0.;
	w.w_yy = 0.;
	w.w_zz = 0.;

	if (q > 2.) {
		return w;

	} else if (q >= 1. && rij > 1e-8) {
	    double q_x  = xij*h1/rij;
	    double q_y  = yij*h1/rij;
	    double q_z  = zij*h1/rij;
	    double q_xy = -xij*yij*h1/rij3;
	    double q_xz = -xij*zij*h1/rij3;
	    double q_yz = -yij*zij*h1/rij3;
	    double q_xx = (yij*yij + zij*zij)*h1/rij3;
	    double q_yy = (xij*xij + zij*zij)*h1/rij3;
	    double q_zz = (xij*xij + yij*yij)*h1/rij3;

		double w_q  = -0.75*(2.-q)*(2.-q);
		double w_qq = 3.-1.5*q;

		w.w =  fac*(0.25*(2.-q)*(2.-q)*(2.-q));

	    w.w_x  = (w_q*q_x*fac);
	    w.w_y  = (w_q*q_y*fac);
	    w.w_z  = (w_q*q_z*fac);
	    w.w_xy = (w_qq*q_x*q_y + w_q*q_xy)*fac;
	    w.w_xz = (w_qq*q_x*q_z + w_q*q_xz)*fac;
	    w.w_yz = (w_qq*q_y*q_z + w_q*q_yz)*fac;
	    w.w_xx = (w_qq*q_x*q_x + w_q*q_xx)*fac;
	    w.w_yy = (w_qq*q_y*q_y + w_q*q_yy)*fac;
	    w.w_zz = (w_qq*q_z*q_z + w_q*q_zz)*fac;

		return w;

	} else {
		double w_q  = -3.0*q*(1.-0.75*q);
		double w_qq = 4.5*q-3.;

	    w.w = fac*(1 - 1.5*q*q*(1-0.5*q));

		if (rij > 1e-8) {
	        double q_x  = xij*h1/rij;
	        double q_y  = yij*h1/rij;
	        double q_z  = zij*h1/rij;
	        double q_xy = -xij*yij*h1/rij3;
	        double q_xz = -xij*zij*h1/rij3;
	        double q_yz = -yij*zij*h1/rij3;
	        double q_xx = (yij*yij + zij*zij)*h1/rij3;
	        double q_yy = (xij*xij + zij*zij)*h1/rij3;
	        double q_zz = (xij*xij + yij*yij)*h1/rij3;

	        w.w_x  = (w_q*q_x*fac);
	        w.w_y  = (w_q*q_y*fac);
	        w.w_z  = (w_q*q_z*fac);
	        w.w_xy = (w_qq*q_x*q_y + w_q*q_xy)*fac;
	        w.w_xz = (w_qq*q_x*q_z + w_q*q_xz)*fac;
	        w.w_yz = (w_qq*q_y*q_z + w_q*q_yz)*fac;
	        w.w_xx = (w_qq*q_x*q_x + w_q*q_xx)*fac;
	        w.w_yy = (w_qq*q_y*q_y + w_q*q_yy)*fac;
	        w.w_zz = (w_qq*q_z*q_z + w_q*q_zz)*fac;

		} else {
			w.w_xx = fac*w_qq*h1*h1;
			w.w_yy = fac*w_qq*h1*h1;
			w.w_zz = fac*w_qq*h1*h1;
		}

		return w;
	}
}

kernel_result wendland_quintic3D(double xi, double yi, double zi, double xj, double yj, double zj, double h) {
	double h1   = 1./h;
	double xij  = xi-xj;
	double yij  = yi-yj;
	double zij  = zi-zj;
	double rij  = sqrt(xij*xij + yij*yij + zij*zij);
	double rij3 = rij*rij*rij;

	double fac = (21./(16.*M_PI))*h1*h1*h1;
	double q   = rij*h1;

	kernel_result w;

	w.w    = 0.;
	w.w_x  = 0.;
	w.w_y  = 0.;
	w.w_z  = 0.;
	w.w_xy = 0.;
	w.w_xz = 0.;
	w.w_yz = 0.;
	w.w_xx = 0.;
	w.w_yy = 0.;
	w.w_zz = 0.;

	if (q < 2.) {
		double tmp = (1.-0.5*q);
		double tmp2 = (q-2.);
		double w_q  = 0.625*tmp2*tmp2*tmp2*q;
		double w_qq = 2.5*tmp2*tmp2*(q-0.5);

		w.w = fac*tmp*tmp*tmp*tmp*(2.*q+1.0);

		if (rij > 1e-8) {
		    double q_x  = xij*h1/rij;
		    double q_y  = yij*h1/rij;
		    double q_z  = zij*h1/rij;
		    double q_xy = -xij*yij*h1/rij3;
		    double q_xz = -xij*zij*h1/rij3;
		    double q_yz = -yij*zij*h1/rij3;
		    double q_xx = (yij*yij + zij*zij)*h1/rij3;
		    double q_yy = (xij*xij + zij*zij)*h1/rij3;
		    double q_zz = (xij*xij + yij*yij)*h1/rij3;

		    w.w_x  = (w_q*q_x*fac);
		    w.w_y  = (w_q*q_y*fac);
		    w.w_z  = (w_q*q_z*fac);
		    w.w_xy = (w_qq*q_x*q_y + w_q*q_xy)*fac;
		    w.w_xz = (w_qq*q_x*q_z + w_q*q_xz)*fac;
		    w.w_yz = (w_qq*q_y*q_z + w_q*q_yz)*fac;
		    w.w_xx = (w_qq*q_x*q_x + w_q*q_xx)*fac;
		    w.w_yy = (w_qq*q_y*q_y + w_q*q_yy)*fac;
		    w.w_zz = (w_qq*q_z*q_z + w_q*q_zz)*fac;
		} else {
			w.w_xx = (w_qq*h1*h1)*fac;
			w.w_yy = (w_qq*h1*h1)*fac;
			w.w_zz = (w_qq*h1*h1)*fac;
		}
		return w;
	}
	return w;
}
