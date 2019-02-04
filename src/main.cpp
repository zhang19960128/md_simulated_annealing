#include "atom.h"
#include <iostream>
#include <fstream>
#include <string>
#include "readion.h"
#include <ctime>
#include "interface.h"
#include "readpara.h"
#include "penalty.h"
#include "sa.h"
#include <time.h>
int main(){
	 readPT("control.PT");
	 int size_box;
	 /*
   SimulatedAnnealing(&PenaltyFunc,
			 control::database[0],
			 control::xop,
			 control::paracount_bvv+control::paracount_charge,
			 saconst::sa_nt,
			 saconst::sa_ns,
			 saconst::sa_max,
			 saconst::sa_temp,
			 saconst::sa_ratio,
			 control::vm,
			 control::ub,
			 control::lb,
			 control::c);
	 */
	  int i=0;
		clock_t start=clock();
		double penaltyp;
		for(size_t k=0;k<2;k++){
			std::cout<<k<<std::endl;
		penaltyp = PenaltyFunc(control::xop,control::database[i],control::ionsize[i],control::minienergytick[i]);//Zhenbang
		}
		clock_t end=clock();
		std::cout<<"the minimum energy tick is"<<control::minienergytick[0]<<" the penalty is"<<penaltyp<<std::endl;
		std::cout<<"the time used is:"<<(double)(end-start)/CLOCKS_PER_SEC<<std::endl;
		return 0;
}
