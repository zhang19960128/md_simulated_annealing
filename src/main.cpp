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
		double penaltyp = PenaltyFunc(control::xop,control::database[i],control::ionsize[i],control::minienergytick[i]);//Zhenbang
		return 0;
}
