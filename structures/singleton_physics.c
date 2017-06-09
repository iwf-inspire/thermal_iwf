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
