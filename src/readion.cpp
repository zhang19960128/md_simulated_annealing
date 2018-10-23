#include "readion.h"
#include "atom.h"
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
box* readion(std::string inputfile,int number){
	std::fstream fs;
	fs.open(inputfile,std::fstream::in);
	std::string line;
	std::istringstream stream1;
	getline(fs,line);
	stream1.str(line);
	int flag=0;
	stream1>>flag;
	stream1.clear();
	box* ionall=new box[flag];
	atom* atomconfig;
	double* period=new double[3];
	double** stress_dft;
	double dftenergy;
	double weight;
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
		}
		ionall[tick].init(atomconfig,number,period,dftenergy,stress_dft,weight);
	}
	return ionall;
}
