#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <stdlib.h>
double** getparameter(std::string pair_style,std::string file){
      std::fstream fs;
      fs.open(file,std::fstream::in);
      std::string line;
      std::istringstream stream1;
      std::string temp;
      int count=0;
      do{
	    getline(fs,line);
	    stream1.str(line);
	    while(stream1>>temp){
               if(temp=="pair_coeff"){
                  stream1>>temp;
                  stream1>>temp;
                  stream1>>temp;
                  if(temp==pair_style){
                     count++;
                  }
               }
         }
         stream1.clear();
	}while(!fs.eof());
        /*test how many atoms are in the input*/
        size_t type=1;
        for(type=1;type<30;type++){
            if(type*(type+1)/2==count){
               break;
            }
        }
        fs.close();
        if(type==30){
           std::cout<<"error, input pair loss or exceed for pair potential: "<<pair_style<<std::endl;
           exit(EXIT_FAILURE);
        }
        else{
          std::cout<<"In this ### "<<pair_style<<" ###: pair potential, we found " <<type<<" types of atoms"<<std::endl;
        }
      double** para=new double* [count];
      fs.open(file,std::fstream::in);
      int tick1,tick2;
      int onedim;
      int* base=new int[type];
      for(size_t j=0;j<4;j++){
         base[j]=type-j;
      }
      do{
	    getline(fs,line);
	    stream1.str(line);
	    while(stream1>>temp){
               if(temp=="pair_coeff"){
                  stream1>>tick1;
                  stream1>>tick2;
                  if(tick1>tick2){
                     std::cout<<"when input the pair parameter, please make the second type tick bigger than first"<<std::endl;
                     exit(EXIT_FAILURE);
                  }
                  stream1>>temp;
                  onedim=tick1*base[tick1-1]+(tick2-tick1);
                  if(pair_style=="12lj/cut/coul/long"){
                     para[onedim]=new double[2];
                     stream1>>para[onedim][0];
                     stream1>>para[onedim][1];
                  }
                  else if(pair_style=="bv"){
                     para[onedim]=new double[5];
                     for(size_t i=0;i<5;i++){
                        stream1>>para[onedim][i];
                     }
                  }
                  else if(pair_style=="bvv"){
                     para[onedim]=new double[5];
                     for(size_t i=0;i<5;i++){
                        stream1>>para[onedim][i];
                     }
                  }
                  else{
                     std::cout<<"unknow pair style"<<std::endl;
                     exit(EXIT_FAILURE);
                  }
               }
         }
         stream1.clear();
	}while(!fs.eof());
        fs.close();
        size_t k=0;
        for(size_t i=0;i<5;i++)
           for(size_t j=i;j<5;j++){
              std::cout<<para[k][0]<<" "<<para[k][1]<<" "<<para[k][2]<<" "<<para[k][3]<<" "<<para[k][4]<<" "<<std::endl;
              k++;
           }
	return para;
}
