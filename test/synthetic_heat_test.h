/*
 * thermal_test.h
 *
 *  Created on: Jan 12, 2016
 *      Author: mamzi
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "../benchmarks/ex1_ex2_reconstruction.h"
#include "../benchmarks/ex3_transient_heat.h"
#include "../benchmarks/ex4_steady_heat.h"

#ifndef THERMAL_TEST_H_
#define THERMAL_TEST_H_

void transient_heat_test(unsigned int N);
void steady_heat_test(unsigned int N);

void full_perform_thermal3D();

#endif /* THERMAL_TEST_H_ */
