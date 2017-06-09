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
