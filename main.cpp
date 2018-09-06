#include "atom.h"
#include <iostream>
#include <fstream>
#include <string>
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
     test.freezeforce();
     test.computebv();
     return 0;
}
