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
	 test.updatelistbv();
	 /*
	 for(size_t i=0;i<40;i++){
	 	std::cout<<"atom "<<i<<" ";
		for(std::list<int>::iterator a=testconfig[i].neibv.begin();a!=testconfig[i].neibv.end();a++){
			std::cout<<*a<<" ";
		}
		std::cout<<std::endl;
	 }
	 */
	 return 0;
}
