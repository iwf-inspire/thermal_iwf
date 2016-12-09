/*
 * reconstruction.h
 *
 *  Created on: Jan 12, 2016
 *      Author: mamzi
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "../benchmarks/ex1_ex2_reconstruction.h"
#include "../methods/cspm.h"
#include "../methods/laplacian.h"
#include "../methods/nmfs.h"
#include "../methods/pse.h"
#include "../methods/rkpm.h"
#include "../methods/sph.h"

#ifndef RECONSTRUCTION_H_
#define RECONSTRUCTION_H_

void run_reconstrunction_3D(unsigned int N);
void full_perform();

#endif /* RECONSTRUCTION_H_ */
