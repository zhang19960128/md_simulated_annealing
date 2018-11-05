#ifdef readpara_h
#define readpara_h
/*read control.PT*/
#include <string>
namespace control{
	double* bvvmatrix;
	double* bvvmatrixmap;
	double* bvvrange;
	double* charge;
	int* type;
}
void readPT(std::string PTfile);
void readvmmap(std::string mapfile);
void readbound(std::string boundfile);
#endif
