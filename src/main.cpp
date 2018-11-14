#include "atom.h"
#include <iostream>
#include <fstream>
#include <string>
#include "readion.h"
#include <ctime>
#include "interface.h"
#include "readpara.h"
#include "penalty.h"
int main(){
	/*
   double** inputbv=getparameter("bv","in.BTO");
	 double** inputbvv=getparameter("bvv","in.BTO");
	 double** inputljcoul=getparameter("12lj/cut/coul/long","in.BTO");
	 atom* testconfig;//the test atom configuration for bond-valence model in lammps
	 testconfig=configuration("mixdata.BTO");
	 double period[3]={8.08,8.08,8.08};
	 box test(testconfig,4,40,period,inputbv,inputbvv,inputljcoul);
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
	 std::clock_t c_end=std::clock();
	 std::cout<<"the total time used is: "<<(c_end-c_start)/CLOCKS_PER_SEC<<std::endl;
   */
	// readion("IonFor.dat",40,3,8);
	 readPT("control.PT");
	 readvmmap("param.map");
	 readbound("param.vmbound");
	// std::cout<<"I am here"<<std::endl;
	 /*
	for(size_t i=0;i<control::pair_num;i++){
		for(size_t j=0;j<12;j++){
			std::cout<<control::bvvmatrixmap[i][j]<<"\t";
		}
		std::cout<<std::endl;
	}
	for(size_t i=0;i<control::paracount_bvv;i++){
		std::cout<<control::lb[i]<<"\t"<<control::ub[i]<<std::endl;
	}
	int i,j;
	indexbvvmap(control::bvvmatrixmap,6,i,j);
	indexchargemap(control::chargemap,8,i);
	std::cout<<i<<" "<<std::endl;
	 return 0;
	 */
}
