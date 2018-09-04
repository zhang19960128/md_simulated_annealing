#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
double** getparameter(std::string pair_style,std::fstream& fs){
	std::string line;
  std::istringstream stream1;
	std::string temp;
	do{
		getline(fs,line);
		std::cout<<line<<std::endl;
//		stream1.str(line);
//		stream1 >> temp;
//		std::cout<<temp<<std::endl;
//		if(temp=="bv"){
//			std::cout<<"Ok I found it"<<std::endl;
//		}
	}while(!fs.eof());
	double** p;
	return p;
}
