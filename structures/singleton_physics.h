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
