#include "atom.h"
#include <iostream>
#include <fstream>
#include <string>
#include "readion.h"
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
				test.freezeforce();
				test.computelong(1e-5);
		 	test.computelj();
				test.computebv();
				test.computebvv();
				test.computestress();

		//		test.printlj();
				std::clock_t c_end=std::clock();
				std::cout<<"the total time used is: "<<(c_end-c_start)/CLOCKS_PER_SEC<<std::endl;
        readion("IonFor.dat",40);
				return 0;
}
