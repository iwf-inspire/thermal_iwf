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

/*
 * kernels_cpu.h
 *
 *  Created on: Sep 16, 2015
 *      Author: mamzi
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#ifndef KERNELS_H_
#define KERNELS_H_

typedef struct kernel_result {
	double w;
	double w_x,w_y,w_z;
	double w_xy,w_xz,w_yz,w_xx,w_yy,w_zz;
} kernel_result;

typedef kernel_result (*kernel)(double, double, double, double, double, double, double);

kernel_result gaussian3D             (double xi, double yi, double zi, double xj, double yj, double zj, double h);
kernel_result cubic_spline3D         (double xi, double yi, double zi, double xj, double yj, double zj, double h);
kernel_result wendland_quintic3D     (double xi, double yi, double zi, double xj, double yj, double zj, double h);

#endif /* KERNELS_H_ */
