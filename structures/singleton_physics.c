#include "../structures/singleton_physics.h"

singleton_physics * get_singleton_physics() {
    static singleton_physics *instance = NULL;
    if(instance == NULL) {
    	instance = malloc(sizeof(*instance));

    	instance->RHO_STL = 7850.0;
    	instance->RHO_ALU = 2770.0;
    	instance->RHO_PUR = 2170.0;

    	instance->KA_STL = 175.0;
    	instance->KA_ALU = 150.0;
    	instance->KA_PUR = 0.100;

    	instance->CP_STL = 434.0;
    	instance->CP_ALU = 875.0;
    	instance->CP_PUR = 1010.0;
    }

    return instance;
}
