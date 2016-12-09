#include "../structures/singleton_geometry.h"
#include "../structures/particle.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifndef BENCHMARKS_H_
#define BENCHMARKS_H_

particle* setup_random_benchmark3D(unsigned int N);
particle* setup_reconstruction_3D(unsigned int N, double hdx);

#endif /* BENCHMARKS_H_ */
