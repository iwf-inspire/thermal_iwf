#include <stdlib.h>

#ifndef SINGLETON_PHYSICS_H_
#define SINGLETON_PHYSICS_H_

typedef struct singleton_physics {
	// density (kg/m3)
	double RHO_STL;
	double RHO_ALU;
	double RHO_PUR;
	// heat conductivity (W/m.k)
	double KA_STL;
	double KA_ALU;
	double KA_PUR;
	// specific heat capacity (J/kg.k)
	double CP_STL;
	double CP_ALU;
	double CP_PUR;

} singleton_physics;

singleton_physics * get_singleton_physics();

#endif /* SINGLETON_PHYSICS_H_ */
