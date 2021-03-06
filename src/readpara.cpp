#include "readpara.h"
#include <string>
#include <fstream>
#include <sstream>
#include "sa.h"
#include <vector>
#include <iostream>
#include <iomanip>
#include <new>
#include <stdlib.h> 
#include "atom.h"
#include "readion.h"
/*read parameters in control.PT,and param.map, this will mainly 
 determined the total parameters.*/
int referenceStruct(box* system, int systemSize){
    int index = -1;
    double reference = 1e25;
    int i = 0;
    for (i=0; i<systemSize; i++){
        if (system[i].dftenergy < reference){
           index = i; 
					 reference=system[i].dftenergy;
        }
    }
    return index;
}
int oldname(std::vector<std::string>& input,std::string& temp){
	int size=input.size();
	for(int i=0;i<size;i++){
		if(input[i]==temp){
			return i;/*if it's a old name, return the reference position*/
		}
	}
	return -1;/*if it's not a old name return -1*/
}
namespace control{
 double** bvvmatrix;/*bvv matrix parameters*/
 int** bvvmatrixmap;/*bvv matrix map change parameters*/
 double* lb;/*lower bound of change variable*/
 double* ub;/*upper bound of change variable*/
 double* vm;/*ub - lb*/
 double* c;/*all values are 2.0*/
 double* charge;/*it's the same length as input species*/
 int* chargemap;/*only maps the asite, bsite, osite charge, so it's length is 3*/
 int* type;/*type matrix*/
 int pair_num;/*how many pair of parameteres in BVV*/
 int paracount_bvv;/*count how many parameteres that need to change in bvv*/
 int paracount_charge;/*count how many parameters that need to change in charge*/
 std::vector<int> para_site_charge_change;/*store the signal whether this site charge change*/
 std::vector<double> para_site_charge;/*store the charge of this site*/
 std::vector<std::string> site_name;/*store the name of this site*/
 std::vector<double> chemical_formula;/*store the chemical formula of this site*/
 double* xop;/*map the all simulation optimized parameter to one line array*/
 /*******************************for fast map************************/
 std::vector<std::vector<int> > mapXpTickToBvvTick;/*fast map the index in xp to BvvMatrix map*/
 std::vector<std::vector<int> > mapXpTickToChargeTick;/*fast map the index in xp to para_site_charge_change*/
 int lastchargetick;/*store the last tick of charge change due to charge-neutral*/
 /***********************************************************************/
 int neutral;/*bool value to show whether you want to force charge neutral*/
 std::vector<std::string> ionfile;/*how many ionfiles are there*/
 std::vector<box*> database;/*the box* data to optimized*/
 std::vector<int> minienergytick;/*store the minimum energy of this Ion files*/
 std::vector<int> ionsize;/*store the structure numbers of different files*/
 std::vector<std::string> deopt;/*output difference energy*/
 std::vector<std::string> dfopt;/*output force difference*/
 std::vector<std::string> dsopt;/*output stress tensor difference*/
}
namespace ewaldsum{
	double cutoff;
	double k_cutoff;
	double alpha;
}
namespace species{
	std::vector<std::string> spe;
	std::vector<int> nametag; /*Go from 0 to 3*/
    std::vector<int> site;/*0 for asite 1 for bsite 2 for O site*/ 
	int** num;/*
						 *the first dimention specifies how many Ion data files are there.
						 * the second dienntion specifies the sequence of the data files.
						 **this is specified in the readion file function*****/
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
std::string decomment(std::string temp){
	std::string s;
	size_t pos_start=0;
	if(temp.find("#")!=std::string::npos){
		s=temp.substr(pos_start,temp.find("#")-pos_start);
	}
	else{
		s=temp;
	}
	return s;
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
				temp_stream.str(flag);
				temp_stream>>keyvalue;
				temp_stream.clear();
			}
}
void findvalue(std::vector<std::string>& input,std::string key,int& keyvalue){
	   std::istringstream temp_stream;
		 std::string flag=input[1];
	   if(input[0].find(key)!=std::string::npos){
				temp_stream.str(flag);
				temp_stream>>keyvalue;
				temp_stream.clear();
			}
}


void readPT(std::string PTfile){
	std::fstream fs;
	fs.open(PTfile,std::fstream::in);
	std::istringstream temp_stream;
	std::vector<std::string> input;
	std::vector<std::string> input_spe;
	int m;
	double charge_temp;
	std::string temp;
	getline(fs,temp);
	std::cout<<"-------------------------------------------STARTING READING CONTROL FILE IN THE INPUT :----------------------"<<std::endl;
	if(temp.find("&sa")!=std::string::npos){
		getline(fs,temp);
		do{
		    temp=decomment(temp);
            if(temp.find_first_not_of("\t\n ")!=std::string::npos){
	            input=split(temp,"=");
	            findvalue(input,"sa_temp",saconst::sa_temp);
	            findvalue(input,"sa_ratio",saconst::sa_ratio);
	            findvalue(input,"sa_eweight",saconst::sa_eweight);
	            findvalue(input,"sa_fweight",saconst::sa_fweight);
	            findvalue(input,"sa_sweight",saconst::sa_sweight);
	            findvalue(input,"sa_nt",saconst::sa_nt);
	            findvalue(input,"sa_ns",saconst::sa_ns);
                findvalue(input,"sa_atom_num",saconst::sa_atom_num);
		}
		getline(fs,temp);
		}while(temp.find("/")==std::string::npos);
	}
	std::cout<<"the Simulated Annealing Parameters are: "<<std::endl;
	std::cout<<"sa_temp: "<<saconst::sa_temp<<std::endl;
	std::cout<<"sa_ratio: "<<saconst::sa_ratio<<std::endl;
	std::cout<<"sa_eweight: "<<saconst::sa_eweight<<std::endl;
	std::cout<<"sa_fweight: "<<saconst::sa_fweight<<std::endl;
	std::cout<<"sa_nt: "<<saconst::sa_nt<<std::endl;
	std::cout<<"sa_ns: "<<saconst::sa_ns<<std::endl;
	std::cout<<"sa_atom_num:"<<saconst::sa_atom_num<<std::endl;
    getline(fs,temp);
	if(temp.find("&ewald")!=std::string::npos){
		getline(fs,temp);
		do{
		temp=decomment(temp);
    if(temp.find_first_not_of("\t\n ")!=std::string::npos){
	  input=split(temp,"=");
		findvalue(input,"r-cutoff",ewaldsum::cutoff);
		findvalue(input,"k-cutoff",ewaldsum::k_cutoff);
		findvalue(input,"alpha",ewaldsum::alpha);
		}
		getline(fs,temp);
		}while(temp.find("/")==std::string::npos);
	}
	std::cout<<"cutoff for ewald is "<<ewaldsum::cutoff<<std::endl;
	std::cout<<"k-cutoff for ewald is "<<ewaldsum::k_cutoff<<std::endl;
	std::cout<<"alpha for ewald is "<<ewaldsum::alpha<<std::endl;
	getline(fs,temp);
	int fraction;
	if(temp.find("&chemical_formula")!=std::string::npos){
	  getline(fs,temp);
		while(temp.find("/")==std::string::npos){
			temp=decomment(temp);
			temp_stream.str(temp);
			temp_stream>>temp;
			temp_stream>>fraction;
			temp_stream.clear();
			control::site_name.push_back(temp);
			control::chemical_formula.push_back(fraction);
			getline(fs,temp);
		};
	}
	for(size_t i=0;i<control::site_name.size();i++){
		std::cout<<"site is: "<<control::site_name[i]<<std::endl;
	}
	std::cout<<"I am here"<<std::endl;
	getline(fs,temp);
	std::vector<double> starting_charge;
	int siteinfo;
	if(temp.find("&species")!=std::string::npos){
		getline(fs,temp);
		do{
		temp=decomment(temp);
    if(temp.find_first_not_of("\t\n ")!=std::string::npos){
			temp_stream.str(temp);
			temp_stream>>temp;
			temp_stream>>m;
			temp_stream>>charge_temp;
			species::spe.push_back(temp);
			species::nametag.push_back(m);
			starting_charge.push_back(charge_temp);
      temp_stream>>temp;
			for(size_t i=0;i<control::site_name.size();i++){
				if(temp==control::site_name[i]){
					species::site.push_back(i);
				}
			}
			temp_stream.clear();
		}
			getline(fs,temp);
		}while(temp.find("/")==std::string::npos);
		std::cout<<"There are many sites: "<<std::endl;
		for(size_t i=0;i<4;i++){
			std::cout<<species::site[i]<<std::endl;
		}
		control::charge=new double [species::spe.size()];
		std::cout<<"the size of the system is: "<<species::spe.size()<<std::endl;
		for(size_t i=0;i<species::spe.size();i++){
			control::charge[i]=starting_charge[i];
		}
	}
	std::cout<<"the species are: "<<std::endl;
	for(size_t i=0;i<species::spe.size();i++){
		std::cout<<species::spe[i]<<"\t"<<species::nametag[i]<<std::endl;
	}
    /*the site information is: */
    std::cout<<"the site information is:---------------------------"<<std::endl;
    for(size_t i=0;i<species::site.size();i++){
        std::cout<<species::site[i]<<"\t";
    }
    std::cout<<std::endl;
	for(size_t i=0;i<species::spe.size();i++)
	{
		std::cout<<control::charge[i]<<" ";
	}
	std::cout<<std::endl;
	getline(fs,temp);
	if(temp.find("&charge")!=std::string::npos){
		do{
		  getline(fs,temp);
			if(temp.find("charge_neutral")!=std::string::npos){
			input=split(temp,"=");
			findvalue(input,"charge_neutral",control::neutral);
			}
		}while(temp.find("/")==std::string::npos);
	}
	if(control::neutral==1){
		std::cout<<"=========================================================FORCE CHARGE NEUTRAL==========================================="<<std::endl;
	}
	else{
		std::cout<<"=================================WARNING!!!!CHARGE NOT NEUTRAL, MAKE SURE THIS IS WAHT YOU WANT========================="<<std::endl;
	};
	getline(fs,temp);
	size_t pair=0;
	size_t count=0;
	size_t tick=0;
	if(temp.find("&bvvmodel")!=std::string::npos){
					pair=species::nametag.size()*(species::nametag.size()+1)/2;
					control::pair_num=pair;
					control::bvvmatrix=new double* [pair];
					getline(fs,temp);
					do{
						temp=decomment(temp);
						if(temp.find_first_not_of("\t\n ")!=std::string::npos){
							control::bvvmatrix[count]=new double[12];
							temp_stream.str(temp);
							temp_stream>>tick;
							temp_stream>>tick;
							for(size_t j=0;j<12;j++)
								temp_stream>>control::bvvmatrix[count][j];
							count++;
						}
						temp_stream.clear();
						getline(fs,temp);
					}while(temp.find("/")==std::string::npos);
	}
	std::cout<<"the starting BDV model parameters are:"<<std::endl;
	for(size_t i=0;i<control::pair_num;i++){
		for(size_t j=0;j<12;j++){
			std::cout<<std::setprecision(7)<<std::setw(7)<<control::bvvmatrix[i][j]<<"\t";
		}
		std::cout<<std::endl;
	}
	getline(fs,temp);
	if(temp.find("&map")!=std::string::npos){
	readvmmap(fs);
	}
	getline(fs,temp);
	while(temp.find("/")==std::string::npos){
		getline(fs,temp);
	};
	getline(fs,temp);
	if(temp.find("&lowbound")!=std::string::npos){
		getline(fs,temp);
		tick=0;
		control::lb=new double [control::paracount_bvv+control::paracount_charge];
		control::c=new double [control::paracount_bvv+control::paracount_charge];
		control::vm=new double [control::paracount_bvv+control::paracount_charge];
		for(size_t i=0;i<control::paracount_bvv+control::paracount_charge;i++){
			control::c[i]=2.0;
		}
		std::cout<<"lower  bound is: ";
		do{
			temp=decomment(temp);
			std::cout<<temp<<"\t";
			if(tick<control::paracount_bvv+control::paracount_charge){
			control::lb[tick]=std::stof(temp);
			}
			tick=tick+1;
			getline(fs,temp);
		}while(temp.find("/")==std::string::npos);
	}
	std::cout<<std::endl;
	getline(fs,temp);
	if(temp.find("&highbound")!=std::string::npos){
		getline(fs,temp);
		tick=0;
		control::ub=new double [control::paracount_bvv+control::paracount_charge];
		std::cout<<"higher bound is: ";
		do{
			temp=decomment(temp);
			std::cout<<temp<<"\t";
			if(tick<control::paracount_bvv+control::paracount_charge){
			control::ub[tick]=std::stof(temp);
			}
			tick=tick+1;
			getline(fs,temp);
		}while(temp.find("/")==std::string::npos);
	}
	std::cout<<std::endl;
	getline(fs,temp);
	std::string systemname;
	if(temp.find("&datafile")!=std::string::npos){
		getline(fs,temp);
		do{
		temp=decomment(temp);
    if(temp.find_first_not_of("\t\n ")!=std::string::npos){
		temp_stream.str(temp);
		temp_stream>>systemname;
		temp_stream>>systemname;
		control::ionfile.push_back(systemname);
		temp_stream>>systemname;
		control::deopt.push_back(systemname);
		temp_stream>>systemname;
		control::dfopt.push_back(systemname);
		temp_stream>>systemname;
		control::dsopt.push_back(systemname);
		temp_stream.clear();
		}
		getline(fs,temp);
		}while(temp.find("/")==std::string::npos);
	}
	for(size_t i=0;i<control::paracount_bvv+control::paracount_charge;i++){
		control::vm[i]=control::ub[i]-control::lb[i];
		}
	std::cout<<"----------------------------------------------------------------------------------------------------------------"<<std::endl;
	std::cout<<"The total Ion file you are reading is: "<<std::endl;
	for(size_t i=0;i<control::ionfile.size();i++){
		std::cout<<control::ionfile[i]<<"\t";
	}
	std::cout<<std::endl;
	std::cout<<"the corresponding file E_difference, F_difference and S_stress are: "<<std::endl;
	for(size_t i=0;i<control::ionfile.size();i++){
		std::cout<<control::deopt[i]<<"\t"<<control::dfopt[i]<<"\t"<<control::dsopt[i]<<std::endl;
	}
	species::num=new int* [control::ionfile.size()];
	std::cout<<std::endl;
	std::cout<<"  AT THIS STAGE, YOU NEED TO READ ION FILES"<<std::endl;
	std::cout<<"-------------------------------------------------------------END------------------------------------------------"<<std::endl;
    std::cout<<"starting readling ion files----------------------------"<<std::endl;
    int temp_size;
    box* temp_box;
    for(size_t i=0;i<control::ionfile.size();i++){
        temp_box=readion(control::ionfile[i],saconst::sa_atom_num,temp_size,ewaldsum::cutoff);
        control::database.push_back(temp_box);
        control::ionsize.push_back(temp_size);
				std::cout<<"the length of database is: "<<temp_size<<" ";
        temp_size=referenceStruct(temp_box,temp_size);
				std::cout<<"the minimum tick of database is: "<<temp_size<<std::endl;
        control::minienergytick.push_back(temp_size);
    }
		std::cout<<"End Preparing process------"<<std::endl;
}
int map_xptick_chargetick(int xptick){
	int sum=0;
	int site=control::para_site_charge_change.size();
	int i=0;
	for(i=0;i<site;i++){
		sum=sum+control::para_site_charge_change[i];
		if(sum-1==xptick-control::paracount_bvv){
			break;
		}
	}
	return i;
}
void readvmmap(std::fstream &fs){
	std::cout<<"---------------------------------------------START READING MAP MATRIX ------------------------------------------"<<std::endl;
	size_t pair=control::pair_num;
	control::bvvmatrixmap=new int* [pair];
	std::istringstream temp_stream;
	std::string temp;
	for(size_t i=0;i<pair;i++){
		control::bvvmatrixmap[i]=new int[12];
		getline(fs,temp);
		temp_stream.str(temp);
		for(size_t j=0;j<12;j++)
			temp_stream>>control::bvvmatrixmap[i][j];
		temp_stream.clear();
	}
	getline(fs,temp);
    /*we fixed that only asite/bsite/osite charge change, so it's only three elements, BEAR IN MIND
    * if you ever find better vaiation principle, you can change it!!!!!!!!!!!!!!!!!!!!!!!!
    */
	control::chargemap=new int [species::spe.size()];
	temp_stream.str(temp);
	for(size_t j=0;j<species::spe.size();j++){
		temp_stream>>control::chargemap[j];
	}
    int sum=0;
    for(size_t i=0;i<species::spe.size();i++){
        sum=sum+control::chargemap[i];
//        std::cout<<sum<<std::endl;
    }
	std::cout<<"the map matrix is: "<<std::endl;
	for(size_t i=0;i<pair;i++){
		for(size_t j=0;j<12;j++){
			std::cout<<control::bvvmatrixmap[i][j]<<"\t";
		}
		std::cout<<std::endl;
	}
	for(size_t i=0;i<3;i++){
		control::para_site_charge_change.push_back(0);
		control::para_site_charge.push_back(0.0);
	}
	for(size_t j=0;j<species::spe.size();j++){
		control::para_site_charge_change[species::site[j]]+=control::chargemap[j];
		control::para_site_charge_change[species::site[j]]=control::para_site_charge_change[species::site[j]]>0?1:0;
		control::para_site_charge[species::site[j]]=control::charge[j];
	}
	for(size_t j=0;j<3;j++){
		std::cout<<"the site: "<<j<<" :"<<"{{charge change(1) or not(0):}} "<<control::para_site_charge_change[j]<<" the starting charge is: "<<control::para_site_charge[j]<<std::endl;
	}
	int temp_sum=0;
	for(size_t i=0;i<3;i++){
		temp_sum=temp_sum+control::para_site_charge_change[i];
	}
	if(temp_sum==3){
		control::paracount_charge=control::neutral==1 ? temp_sum-1:temp_sum;
	}
	else if(temp_sum==2){
		control::paracount_charge=control::neutral==1 ? temp_sum-1:temp_sum;
	}
	else if(temp_sum==1){
		std::cout<<"ONLY ONE SITE CHARGE CHANGE"<<std::endl;
		exit(EXIT_FAILURE);
	}
	else if(temp_sum==0){
		control::paracount_charge=0;
	}
    /*We only consider all the site charge change or all the site charge not change*/
  /*
	if(!(sum == species::spe.size()-1 || sum ==0)){
       std::cout<<species::spe.size()<<"  "<<"not charge change neutrol"<<std::endl;
       std::cout<<"!!!!!!!!!!!!!1If you specify one atom on Asite/Bsite/Osite charge change, please also specify other elements that on the same site change charge too!, this is very important!!!!"<<std::endl;
       exit(EXIT_FAILURE);
    }
		*/
	std::cout<<"------------------------------------------------------------END--------------------------------------------------"<<std::endl;
	/*******starting to map XpTickToBvvTick*******/
	/*store the first map function within control::paracount_bvv*/
	std::vector<int> tempxpbvv(3,0);
	int count=0;
	for(size_t m=0;m<control::pair_num;m++)
		for(size_t n=0;n<12;n++){
			if(control::bvvmatrixmap[m][n]==0)
				continue;
			else{
				tempxpbvv[0]=count;
				tempxpbvv[1]=m;
				tempxpbvv[2]=n;
				control::mapXpTickToBvvTick.push_back(tempxpbvv);
				count++;
			}
		}
	control::paracount_bvv=count;
	/*store the second map function within control::paracount_charge*/
	std::vector<int> tempxpcharge(2,0);
	count=0;
	if(control::paracount_charge==0){
	}
	else{
		if(control::neutral==1){
			for(size_t i=control::paracount_bvv;i<control::paracount_charge+control::paracount_bvv;i++){
				tempxpcharge[0]=i;
				tempxpcharge[1]=map_xptick_chargetick(i);
				control::mapXpTickToChargeTick.push_back(tempxpcharge);
			}
			for(size_t start=species::site.size()-1;start>=0;start--){
				if(control::para_site_charge_change[start]==1){
				control::lastchargetick=start;
				break;
				}
			}
		}
		else{
			for(size_t i=control::paracount_bvv;i<control::paracount_charge+control::paracount_bvv;i++){
				tempxpcharge[0]=i;
				tempxpcharge[1]=map_xptick_chargetick(i);
				control::mapXpTickToChargeTick.push_back(tempxpcharge);
			}
		}
	}
	/*malloc the storage space for control::paracount_charge and control::paracount_bvv*/
  control::xop=new double [control::paracount_charge+control::paracount_bvv];
	for(size_t i=0;i<control::paracount_bvv;i++){
		control::xop[i]=control::bvvmatrix[control::mapXpTickToBvvTick[i][1]][control::mapXpTickToBvvTick[i][2]];
	}
	for(size_t j=control::paracount_bvv;j<control::paracount_bvv+control::paracount_charge;j++){
		control::xop[j]=control::para_site_charge[control::mapXpTickToChargeTick[j-control::paracount_bvv][1]];
	}
	for(size_t i=0;i<control::paracount_charge+control::paracount_bvv;i++){
		std::cout<<control::xop[i]<<std::endl;
	}
	/*debug*/
	std::cout<<"MAP xp to what we want........"<<std::endl;
	for(size_t i=0;i<control::mapXpTickToBvvTick.size();i++){
		std::cout<<std::endl;
		for(size_t j=0;j<control::mapXpTickToBvvTick[i].size();j++){
			std::cout<<control::mapXpTickToBvvTick[i][j]<<"\t";
		}
	}
	for(size_t i=0;i<control::mapXpTickToChargeTick.size();i++){
		std::cout<<std::endl;
		for(size_t j=0;j<control::mapXpTickToChargeTick[i].size();j++){
			std::cout<<control::mapXpTickToChargeTick[i][j]<<"\t"	;
		}
	}
	std::cout<<std::endl<<"Exit the function without any erros"<<std::endl;
	/**************** end storing map function **/
};
