#include "../base/utils.h"
#include "../test/eng_heat_model.h"
#include "../test/reconstruction_test.h"
#include "../test/synthetic_heat_test.h"

int main() {
	// ******************************
	// example 1: polynomial function
	// ******************************
	// run_reconstrunction_3D(20);
	// full_perform();

	// *********************************
	// example 2: trigonometric function
	// *********************************
	//run_thermal_cpu3D(21);

	// ***********************************
	// example 3: synthetic transient heat
	// ***********************************
	transient_heat_test(21);

	// ********************************
	// example 4: synthetic steady heat
	// ********************************
	// steady_heat_test(21);

	// ***********************************
	// example 5: engineering model
	// ***********************************
	// run_eng_model3D();
	// full_perform_thermal3D();

	printf("\n");
	printf("-----------------finished successfully----------------");
	return EXIT_SUCCESS;
}
