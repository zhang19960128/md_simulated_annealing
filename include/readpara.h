#ifndef readpara_h
#define readpara_h
/*read control.PT*/
#include <string>
#include <vector>
#include <iostream>
namespace control{
extern double** bvvmatrix;
extern int** bvvmatrixmap;
extern double* lb;
extern double* ub;
extern double* charge;
extern int* chargemap;
extern int* type;
extern int pair_num;
extern int paracount_bvv;
extern int paracount_charge;
extern double* xop;
extern std::vector<std::string> ionfile;
}
namespace ewaldsum{
extern double cutoff;
extern double k_cutoff;
extern double alpha;
}
namespace species{
	extern std::vector<std::string> spe;
	extern std::vector<int> num;
}
void readPT(std::string PTfile);
void readvmmap(std::string mapfile);
void readbound(std::string boundfile);
#endif
