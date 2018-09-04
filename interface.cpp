#include <fstream>
#include <string>
#include <sstream>
double** getparameter(std::string pair_style,std::fstream fs){
	std::string line;
	istringstream stream1;
	std::string temp;
	do{
		getline(fs,line);
		stream1.str(line);
		stream1 >> temp;
		if(temp=="pair_coeff"){
			std::cout<<"Ok I found it"<<std::endl;
		}
	}while(line);
	double** p;
	return p;
}
