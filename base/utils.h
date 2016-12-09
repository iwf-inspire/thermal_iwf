#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_linalg.h>

#ifndef UTILS_H_
#define UTILS_H_

void inv3x3(double A[9], double Ainv[9]);

double det4x4(double M[4*4]);
double det5x5(double M[5*5]);
double det6x6(double M[6*6]);
double det50x50_limp(double M[50][50]);

void solve3x3(double M[3*3], double P[3], double C[3]);
void solve4x4(double M[4*4], double P[4], double C[4]);
void solve6x6(double M[6*6], double P[6], double C[6]);

void solve4x4_gsl(double A[4*4], double B[4], double res[4]);
void solve6x6_gsl(double A[6*6], double B[6], double res[6]);
void solve9x9_gsl(double A[9*9], double B[9], double res[9]);

int solve_quadratic(double coeff[3], double sol[2]);
int mod(int a, int b);

#endif /* UTILS_H_ */
