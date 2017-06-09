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

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef UTILS_H_
#define UTILS_H_

void inv3x3(double A[9], double Ainv[9]);

void solve3x3(double M[3*3], double P[3], double C[3]);
void solve4x4(double M[4*4], double P[4], double C[4]);
void solve6x6(double M[6*6], double P[6], double C[6]);

#endif /* UTILS_H_ */
