#ifndef readpara_h
#define readpara_h
/*read control.PT*/
#include <string>
#include <vector>
#include "atom.h"
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
extern std::vector<box*> database;
extern std::vector<int> minienergytick;/*store the minimum energy of this Ion files*/
extern std::vector<int> ionsize;/*store the structure numbers of different files*/
extern std::vector<std::string> deopt;
extern std::vector<std::string> dfopt;
extern std::vector<std::string> dsopt;
}
namespace ewaldsum{
extern double cutoff;
extern double k_cutoff;
extern double alpha;
}
namespace species{
	extern std::vector<std::string> spe;
	extern std::vector<int> nametag;
    extern std::vector<int> site;/*0 for asite 1 for bsite o for O site*/ 
	extern int** num;
}
void readPT(std::string PTfile);
void readvmmap(std::string mapfile);
void readbound(std::string boundfile);
#endif
