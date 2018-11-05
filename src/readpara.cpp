#include "readpara.h"
#include <string>
#include <fstream>
#include <sstream>
#include "sa.h"
#include <vector>
/*read parameters in control.PT,and param.map, this will mainly 
 determined the total parameters.*/
namespace control{
	extern double* bvvmatrix;
	extern double* bvvmatrixmap;
	extern double* bvvrange;
	extern double* charge;
	extern int* type;
}
std::vector<std::string> split(std::string s, std::string delimiter) {
		    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
				std::string token;
				std::vector<std::string> res;
				while ((pos_end = s.find (delimiter, pos_start)) != std::string::npos) {
					token = s.substr(pos_start, pos_end - pos_start);
					pos_start = pos_end + delim_len;
					res.push_back(token);
				}
				res.push_back(s.substr(pos_start));
			  return res;
}
double findvalue(std::vector<std::string>& input,std::string key){
	   std::istringstream temp_stream;
		 std::string flag=input[1];
		 double keyvalue;
	   if(input[0].find(key)!=std::string::npos){
				flag=flag.substr(0,flag.find(","));
				temp_stream.str(flag);
				temp_stream>>keyvalue;
				temp_stream.clear();
		    return keyvalue;
			}
		 else{
		 	return -1;
		 }
}
void readPT(std::string PTfile){
	std::fstream fs;
	fs.open(PTfile,std::fstream::in);
	/*reading sa parameters*/
	std::string temp;
	std::istringstream temp_stream;
	getline(fs,temp);
	std::vector<std::string> input;
	while(!fs.eof()){
		std::cout<<temp<<std::endl;
		getline(fs,temp);
		input=split(temp,"=");
		/*reading parameters for simulated annealing*/
		if(input.size()>=1){
			saconst::sa_temp=findvalue(input,"sa_temp");
			saconst::sa_ratio=findvalue(input,"sa_ratio");
			saconst::sa_max=findvalue(input,"sa_eweigh");
			saconst::sa_fweight=findvalue(input,"sa_fweight");
			saconst::sa_sweight=findvalue(input,"sa_sweight");
			saconst::sa_nt=findvalue(input,"sa_nt");
			saconst::sa_ns=findvalue(input,"sa_ns");
		}
	}
	/*end reading sa parameters*/
};
void readvmmap(std::string mapfile){
};
void readbound(std::string boundfile){
}
