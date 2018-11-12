#ifndef readpara_h
#define readpara_h
/*read control.PT*/
#include <string>
namespace control{
extern double* bvvmatrix;
extern double* bvvmatrixmap;
extern double* bvvrange;
extern double* charge;
extern int* type;
}
void readPT(std::string PTfile);
void readvmmap(std::string mapfile);
void readbound(std::string boundfile);
#endif
