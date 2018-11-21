#include "atom.h"
#include <string>
#include <iostream>
#include <fstream>
#include "sa.h"
#include <cmath>
#include <iomanip>
/*
 * size is the size ot input box, tick is the minimum DFT energy, deoptfile is the difference energy file, dfoptfile is the difference Force file, dsoptfile is the difference stress file
 * 
 */
void writeoutput(box* input,int size,int tick,int saiter,std::string deoptfile,std::string dfoptfile,std::string dsoptfile){
	std::fstream dfopt,deopt,dsopt;
	dfopt.open(dfoptfile,std::fstream::out);
	deopt.open(deoptfile,std::fstream::out);
	dsopt.open(dsoptfile,std::fstream::out);
	deopt<<"#"<<saconst::sa_temp<<std::endl;
	deopt<<saiter+1<<std::endl;
	deopt<<std::setw(18)<<"ABS(E(DFT)-E(MD))\t"<<std::setw(18)<<"E(DFT)\t"<<std::setw(18)<<"E(MD)"<<std::endl;
	dfopt<<std::setw(18)<<"ABS(F(DFT)-F(MD))\t"<<std::setw(18)<<"F(DFT)\t"<<std::setw(18)<<"F(MD)"<<std::endl;
	double dftetemp;
	double mdetemp;
	for(size_t i=0;i<size;i++){
		mdetemp=input[i].mdenergy-input[tick].mdenergy;
		dftetemp=input[i].dftenergy-input[tick].dftenergy;
		deopt<<std::setw(18)<<"\t"<<std::fabs(dftetemp-mdetemp)<<std::setw(18)<<"\t"<<dftetemp<<std::setw(18)<<"\t"<<mdetemp<<std::endl;
		for(size_t j=0;j<saconst::sa_atom_num;j++){
			for(size_t k=0;k<3;k++){
			dfopt<<std::setw(10)<<"\t"<<fabs(input[i].allatom[j].force[k]-input[i].allatom[j].dftforce[k]);
			}
			}
		for(size_t j=0;j<saconst::sa_atom_num;j++){
			for(size_t k=0;k<3;k++){
			dfopt<<std::setw(10)<<"\t"<<input[i].allatom[j].dftforce[k];
			}
			}
		for(size_t j=0;j<saconst::sa_atom_num;j++){
			for(size_t k=0;k<3;k++){
			dfopt<<std::setw(10)<<"\t"<<input[i].allatom[j].force[k];
			}
			}
		std::cout<<std::endl;
	}
}
