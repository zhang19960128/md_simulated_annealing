#include "atom.h"
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include "interface.h"
int main(){
   double** inputbv=getparameter("bv","in.BTO");
	 double** inputbvv=getparameter("bvv","in.BTO");
	 double** inputljcoul=getparameter("12lj/cut/coul/long","in.BTO");
	 /*input the configuration for lammps simulation*/
	 atom* testconfig;//the test atom configuration for bond-valence model in lammps
	 testconfig=configuration("mixdata.BTO");
	 double period[3]={8.08,8.08,8.08};
	 box test(testconfig,4,40,period,inputbv,inputbvv,inputljcoul);

	 /*this testconfig is useless after initialization*/
	 test.updatelistbv();
	 test.updatelistbvv();
         test.updatelistlj();
         std::clock_t c_start=std::clock();
	 for(size_t i=0;i<1;i++){
		 	test.freezeforce();
   		test.computebv();
   		test.computebvv();
	    test.computelj();
//        test.computelong();
	 }
   return 0;
}
