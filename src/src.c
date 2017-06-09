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

#include "../test/eng_heat_model.h"
#include "../test/reconstruction_test.h"
#include "../test/synthetic_heat_test.h"

int main(void) {
	// choose your mesh-free scheme among: {SPH, PSE, RKPM, CSPM, NMFS}
	METHOD meshfree_method = NMFS;
	// define the number of particles in one direction
	unsigned int Nx = 11;

	// ******************************
	// example 1: polynomial function
	// ******************************
	run_reconstrunction_polynomial_3D(meshfree_method,Nx,false);

	// *********************************
	// example 2: trigonometric function
	// *********************************
	// run_reconstrunction_trigonometric_3D(meshfree_method,Nx,false);

	// ***********************************
	// example 3: synthetic transient heat
	// ***********************************
	// run_transient_heat_3D(meshfree_method,Nx,false);

	// ********************************
	// example 4: synthetic steady heat
	// ********************************
	// run_steady_heat_3D(meshfree_method,Nx,false);

	// ***********************************
	// example 5: engineering model
	// ***********************************
	// run_eng_model_3D(meshfree_method);

	printf("\n");
	printf("-----------------finished successfully----------------");
	return EXIT_SUCCESS;
}
