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
	 double period[3]={40.4,40.4,40.4};
	 box test(testconfig,4,40,period,inputbv,inputbvv);
	 test.freezeforce();
	 test.
	 return 0;
}
