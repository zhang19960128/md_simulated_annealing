#include "readpara.h"
#include <string>
#include <fstream>
#include <sstream>
#include "sa.h"
#include <vector>
/*read parameters in control.PT,and param.map, this will mainly 
 determined the total parameters.*/
namespace control{
 double* bvvmatrix;
 double* bvvmatrixmap;
 double* bvvrange;
 double* charge;
 int* type;
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
void findvalue(std::vector<std::string>& input,std::string key,double& keyvalue){
	   std::istringstream temp_stream;
		 std::string flag=input[1];
	   if(input[0].find(key)!=std::string::npos){
			 if(flag.find(",")!=std::string::npos){
				flag=flag.substr(0,flag.find(","));
			 }
				temp_stream.str(flag);
				temp_stream>>keyvalue;
				temp_stream.clear();
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
		getline(fs,temp);
		std::cout<<temp<<std::endl;
		if(!temp.empty()){
	    	input=split(temp,"=");
				if(input.size()>1){
	    	findvalue(input,"sa_temp",saconst::sa_temp);
	    	findvalue(input,"sa_ratio",saconst::sa_ratio);
	    	findvalue(input,"sa_eweigh",saconst::sa_max);
	    	findvalue(input,"sa_fweight",saconst::sa_fweight);
	    	findvalue(input,"sa_sweight",saconst::sa_sweight);
	    	findvalue(input,"sa_nt",saconst::sa_nt);
	    	findvalue(input,"sa_ns",saconst::sa_ns);
				}
		}
		else{
			continue;
		}
		/*reading parameters for simulated annealing*/
	}
	std::cout<<"sa_temp: "<<saconst::sa_temp<<std::endl;
	/*end reading sa parameters*/
};
void readvmmap(std::string mapfile){
};
void readbound(std::string boundfile){
}
