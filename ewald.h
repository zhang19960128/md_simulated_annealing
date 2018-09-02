#ifndef ewald.h
#define ewald.h
#include <math.h>
#include <Eigen/Dense>
#include <stdlib.h>
#include <complex>
double EwaldSum(atom* input, int size, double* lattice, double sigma, double rcut, double gcut);
#endif


