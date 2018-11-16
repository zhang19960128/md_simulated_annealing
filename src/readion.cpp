#include "readion.h"
#include "atom.h"
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include "readpara.h"
box* readion(std::string inputfile,int number,int& boxnumber,double cutoff){
	std::fstream fs;
	std::cout<<"---------------------------------------STARTING READING ION COORDINATES----------------------------"<<std::endl;
	fs.open(inputfile,std::fstream::in);
	std::string line;
	std::istringstream stream1;
	getline(fs,line);
	stream1.str(line);
	int flag=0;
	stream1>>flag;
	boxnumber=flag;
	stream1.clear();
    std::cout<<"There are how many atoms "<<flag<<std::endl;
	box* ionall=new box[flag];
	atom* atomconfig;
	double* period=new double[3];
	double** stress_dft;
	double dftenergy;
	double weight;
	std::cout<<"you are READING "<<inputfile<<" ionic files"<<std::endl;
	/************************set the type tick**************************/
	/*the type tick go from 0,1,2,3,4,5,6,7,.....*/
	int* type_tick=new int [number];
	for(size_t i=0;i<number;i++){
		getline(fs,line);
		for(size_t j=0;j<species::spe.size();j++){
			if(line.find(species::spe[j])!=std::string::npos){
				type_tick[i]=j;
				break;
			}
		}
	}
	std::cout<<"the input sequence for atom is: "<<std::endl;
	for(size_t i=0;i<number;i++){
		if(i%5==0){
			std::cout<<std::endl;
		}
		std::cout<<type_tick[i]+1<<"\t";
	}
	std::cout<<std::endl;
	/*******************************************************************/
	for(size_t tick=0;tick<flag;tick++){
		getline(fs,line);
		atomconfig=new atom [number];
		for(size_t j=0;j<number;j++){
			getline(fs,line);
			stream1.str(line);
		/*reading atomic positions*/
			for(size_t k=0;k<3;k++){
				stream1>>atomconfig[j].position[k];
			}
			stream1.clear();
			/*end reading that*/
		}
		getline(fs,line);
		/*reading periodical boudary condition*/
		for(size_t k=0;k<3;k++){
				getline(fs,line);
				stream1.str(line);
				for(size_t m=0;m<=k;m++){
						stream1>>period[k];
				}
				stream1.clear();
			}
		/*end reading that*/
		getline(fs,line);
		getline(fs,line);
		getline(fs,line);
		/*readling forces*/
		for(size_t j=0;j<number;j++){
			getline(fs,line);
			stream1.str(line);
			for(size_t k=0;k<3;k++){
				stream1>>atomconfig[j].dftforce[k];
			}
			stream1.clear();
		}
			/*end reading that*/
		getline(fs,line);
		/*reading energy*/
		getline(fs,line);
		stream1.str(line);
		stream1>>dftenergy;
		stream1.clear();
		/*end reading energy*/
		getline(fs,line);
		/*reading stress tensor*/
		stress_dft=new double* [3];
		for(size_t i=0;i<3;i++){
			stress_dft[i]=new double [3];
			getline(fs,line);
			stream1.str(line);
			for(size_t j=0;j<3;j++){
				stream1>>stress_dft[i][j];
			}
			stream1.clear();
		}
		/*end reading stress tensor*/
		getline(fs,line);
		/*reading weight*/
		getline(fs,line);
		stream1.str(line);
		stream1>>weight;
		for(size_t i=0;i<number;i++){
			for(size_t j=0;j<3;j++){
				atomconfig[i].position[j]=atomconfig[i].position[j]*period[j];
			}
			atomconfig[i].type=type_tick[i];
		}
		ionall[tick].init(atomconfig,number,species::spe.size(),cutoff,period,dftenergy,stress_dft,weight);
		delete [] atomconfig;
	}
	delete [] type_tick;
	std::cout<<"---------------------------------------------------END--------------------------------------------"<<std::endl;
	return ionall;
}
