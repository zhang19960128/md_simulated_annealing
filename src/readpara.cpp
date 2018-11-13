#include "readpara.h"
#include <string>
#include <fstream>
#include <sstream>
#include "sa.h"
#include <vector>
#include <iostream>
#include <new>
/*read parameters in control.PT,and param.map, this will mainly 
 determined the total parameters.*/
namespace control{
 double** bvvmatrix;
 int** bvvmatrixmap;
 double* lb;
 double* ub;
 double* charge;
 int* chargemap;
 int* type;
 int pair_num;
 int paracount;
}
namespace ewaldsum{
	double cutoff;
	double k_cutoff;
	double alpha;
}
namespace species{
	std::vector<std::string> spe;
	std::vector<int> num;
}
std::vector<std::string> split(std::string temp, std::string delimiter) {
		    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
				std::string token;
				std::string s;
				if(temp.find("#")!=std::string::npos){
					s=temp.substr(pos_start,temp.find("#")-pos_start);
				}
				else{
					s=temp;
				}
				std::vector<std::string> res;
				while ((pos_end = s.find (delimiter, pos_start)) != std::string::npos) {
					token = s.substr(pos_start, pos_end - pos_start);
					pos_start = pos_end + delim_len;
					res.push_back(token);
				}
				res.push_back(s.substr(pos_start));
			  return res;
}
std::vector<std::string> splitspace(std::string& s){
	std::vector<std::string> result;
	std::istringstream iss(s);
	for(std::string s; iss >> s; )
		    result.push_back(s);
	return result;
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
	std::vector<std::string> input_spe;
	size_t elements=0;
	size_t pair;
	int tick;
	while(!fs.eof()){
		getline(fs,temp);
		//std::cout<<temp<<std::endl;
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
				findvalue(input,"cutoff",ewaldsum::cutoff);
				findvalue(input,"k-cutoff",ewaldsum::k_cutoff);
				findvalue(input,"alpha",ewaldsum::alpha);
				}
				input_spe=splitspace(temp);
				if(input_spe.size()==2 && input.size()!=2){
					species::spe.push_back(input_spe[0]);
					std::cout<<"the first is: "<<input_spe[0]<<" the second is "<<input_spe[1]<<std::endl;
        	temp_stream.str(input_spe[1]);
					temp_stream>>tick;
					temp_stream.clear();
					species::num.push_back(tick);
				}
				if(temp.find("bvmodel")!=std::string::npos){
					getline(fs,temp);
					while(temp.substr(0,1).find('#')!=std::string::npos){
						getline(fs,temp);
					}
					pair=species::num.size()*(species::num.size()+1)/2;
					control::pair_num=pair;
					control::bvvmatrix=new double* [pair];
					for(size_t i=0;i<pair;i++){
						control::bvvmatrix[i]=new double[10];
						temp_stream.str(temp);
						std::cout<<temp<<std::endl;
						temp_stream>>tick;
						temp_stream>>tick;
						for(size_t j=0;j<10;j++)
							temp_stream>>control::bvvmatrix[i][j];
						temp_stream.clear();
						getline(fs,temp);
					}
				}
		}
		else{
			continue;
		}
		/*reading parameters for simulated annealing*/
	}
	/*end reading sa parameters*/
};
void readvmmap(std::string mapfile){
	std::fstream fs;
	fs.open(mapfile,std::fstream::in);
	size_t pair=species::num.size();
	pair=(pair+1)*pair/2;
	control::bvvmatrixmap=new int* [pair];
	std::istringstream temp_stream;
	std::string temp;
	for(size_t i=0;i<pair;i++){
		control::bvvmatrixmap[i]=new int[10];
		getline(fs,temp);
		temp_stream.str(temp);
		for(size_t j=0;j<10;j++)
			temp_stream>>control::bvvmatrixmap[i][j];
		temp_stream.clear();
	}
	getline(fs,temp);
	control::chargemap=new int [species::num.size()];
	temp_stream.str(temp);
	for(size_t j=0;j<species::num.size();j++){
		temp_stream>>control::chargemap[j];
	}
};
void readbound(std::string boundfile){
	std::fstream fs;
	fs.open(boundfile,std::fstream::in);
	int sum=0;
	std::istringstream temp_stream;
	std::string temp;
	for(size_t i=0;i<control::pair_num;i++)
		for(size_t j=0;j<10;j++){
			sum=sum+control::bvvmatrixmap[i][j];
		}
	for(size_t i=0;i<species::num.size();i++)
		sum=sum+control::chargemap[i];
	control::paracount=sum;
	control::lb=new double [sum];
	control::ub=new double [sum];
	double rang;
	std::cout<<"the total variable need to change is "<<sum<<std::endl;
	for(size_t i=0;i<sum;i++){
		getline(fs,temp);
		temp_stream.str(temp);
		temp_stream>>rang;
		control::lb[i]=rang;
		temp_stream>>rang;
		control::ub[i]=rang;
		temp_stream.clear();
	}
}
