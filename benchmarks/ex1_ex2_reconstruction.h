#include "../structures/singleton_geometry.h"
#include "../structures/particle.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <time.h>

#ifndef BENCHMARKS_H_
#define BENCHMARKS_H_

particle* setup_reconstruction_polynomial_3D(unsigned int N, double hdx, bool RANDOM);
particle* setup_reconstruction_trigonometric_3D(unsigned int N, double hdx, bool RANDOM);

#endif /* BENCHMARKS_H_ */
