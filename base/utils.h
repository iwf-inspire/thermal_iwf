#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef UTILS_H_
#define UTILS_H_

void inv3x3(double A[9], double Ainv[9]);

void solve3x3(double M[3*3], double P[3], double C[3]);
void solve4x4(double M[4*4], double P[4], double C[4]);
void solve6x6(double M[6*6], double P[6], double C[6]);

#endif /* UTILS_H_ */
