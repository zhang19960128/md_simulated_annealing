#ifndef sa_h
#define sa_h
#include <iostream>
namespace saconst{
	double sa_temp=0.0001;
	double sa_ratio=0.95;
	double sa_max=3;
	double sa_eweight=5.0;
	double sa_fweight=1.0;
	double sa_sweight=0.0;
	double sa_nt=3;
	double sa_ns=3;
}
namespace control{
	double* bvvmatrix;
	double* bvvmatrixmap;
	double* bvvrange;
	double* charge;
	int* type;
}
#endif
