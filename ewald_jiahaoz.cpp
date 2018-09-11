#include "atom.h"
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <list>
#include <iomanip>
#define pi 3.14159265359
/*Perform the Ewald summation. The lattice paramters should be given*/
/*Lattice should be a flattened array*/
void box::computelong(){
    double ShortRange = 0;
    double LongRange = 0;
    double selfe = 0;
    double epsil = 0;
		double coul_prefactor=180.9512801302711;
		/*e/(epsilon0*A)=180.951  ev when use U=k*q1*q2/r*/
    double volume=p[0]*p[1]*p[2];
    double delx,dely,delz,rsq,r;
    double root2pi=pow(2*pi,0.5);
    double* xall=new double[size];
    double* yall=new double [size];
    double* zall=new double [size];
    double* allcharge=new double [size];
    for(size_t i=0;i<size;i++){
       xall[i]=allatom[i].position[0];
       yall[i]=allatom[i].position[1];
       zall[i]=allatom[i].position[2];
       allcharge[i]=allatom[i].charge;
       for(std::list<int>::iterator j=allatom[i].neilj.begin();j!=allatom[i].neilj.end();j++){
         delx=allatom[i].position[0]-virtatom[*j].position[0];
         dely=allatom[i].position[1]-virtatom[*j].position[1];
         delz=allatom[i].position[2]-virtatom[*j].position[2];
         rsq=delx*delx+dely*dely+delz*delz;
         r=sqrt(rsq);
         ShortRange+=coul_prefactor/4/pi*allatom[i].charge*virtatom[*j].charge/r*erfc(r/sqrt(2)/sigma);
      }
    }
    ShortRange=ShortRange/2.0;
		//std::cout<<"the short range coulumb potential is: "<<ShortRange<<std::endl;
    for(size_t i=0;i<size;i++){
      selfe=selfe-1/root2pi/sigma/4/pi*coul_prefactor*(allatom[i].charge)*allatom[i].charge;
    }
		//std::cout<<"the self energy is: "<<selfe<<std::endl;
    /*calculate the longrange energy*/
    /*exp(i*k_v*r)
     *=exp(i*k_x*r_x)*exp(i*k_y*r_y)*exp(i*k_z*r_z)
     *=cos(k_x*r_x)*cos(k_y*r_y)*cos(k_z*r_z)-sin(k_x*r_x)*sin(k_y*r_y)*cos(k_z*r_z)-sin(k_x*r_x)*cos(k_y*r_y)*cos(k_z*r_z)-cos(k_x*r_x)*sin(k_y*r_y)*sin(k_z*r_z)
     *+I(********************************************************) similar term just expand it
     *
     * */
    /*define array for cos(k_x*r_x),cos(k_y*r_y),cos(k_z*r_z),sin(k_x*r_x),sin(k_y*r_y),sin(k_z*r_z)*/
    double** cskxrx=new double* [gmax+1];
    double** cskyry=new double* [gmax+1];
    double** cskzrz=new double* [gmax+1];
    double** snkxrx=new double* [gmax+1];
    double** snkyry=new double* [gmax+1];
    double** snkzrz=new double* [gmax+1];
    for(size_t i=0;i<=gmax;i++){
      cskxrx[i]=new double [size];
      cskyry[i]=new double [size];
      cskzrz[i]=new double [size];
      snkxrx[i]=new double [size];
      snkyry[i]=new double [size];
      snkzrz[i]=new double [size];
    }
    for(int i=0;i<=gmax;i++)
       for(size_t j=0;j<size;j++){
          cskxrx[i][j]=cos(2*pi/p[0]*i*xall[j]);
          cskyry[i][j]=cos(2*pi/p[1]*i*yall[j]);
          cskzrz[i][j]=cos(2*pi/p[2]*i*zall[j]);
          snkxrx[i][j]=sin(2*pi/p[0]*i*xall[j]);
          snkyry[i][j]=sin(2*pi/p[1]*i*yall[j]);
          snkzrz[i][j]=sin(2*pi/p[2]*i*zall[j]);
			//		std::cout<<cskxrx[i][j]<<" "<<cskyry[i][j]<<" "<<cskzrz[i][j]<<std::endl;
				  //std::cout<<xall[j]<<" "<<yall[j]<<" "<<zall[j]<<std::endl;
       }
    double skre,skim,skmodsq,ksq;
    for(int h=-1*gmax;h<=gmax;h++)
       for(int k=-1*gmax;k<=gmax;k++)
          for(int l=-1*gmax;l<=gmax;l++){
             skre=0.0;
             skim=0.0;
	     if(h==0&&k==0&&l==0){
		  continue;
	       }
             for(size_t i=0;i<size;i++){
	        /*finish computing real part*/
                skre=skre+(cskxrx[abs(h)][i]*cskyry[abs(k)][i]*cskzrz[abs(l)][i])*allcharge[i];
                skre=skre-allcharge[i]*cskxrx[abs(h)][i]*snkyry[abs(k)][i]*snkzrz[abs(l)][i]*(k>0 ? 1:-1)*(l>0 ? 1:-1);
	        skre=skre-allcharge[i]*snkxrx[abs(h)][i]*cskyry[abs(k)][i]*snkzrz[abs(l)][i]*(h>0 ? 1:-1)*(l>0 ? 1:-1);
	        skre=skre-allcharge[i]*snkxrx[abs(h)][i]*snkyry[abs(k)][i]*cskzrz[abs(l)][i]*(h>0 ? 1:-1)*(k>0 ? 1:-1);
	        /*start to compute the imaginary part*/
	        skim=skim+cskxrx[abs(h)][i]*cskyry[abs(k)][i]*snkzrz[abs(l)][i]*(l>0 ? 1:-1);
	        skim=skim+snkxrx[abs(h)][i]*cskyry[abs(k)][i]*cskzrz[abs(l)][i]*(h>0 ? 1:-1);
	        skim=skim+cskxrx[abs(h)][i]*snkyry[abs(k)][i]*cskzrz[abs(l)][i]*(k>0 ? 1:-1);
	       skim=skim-snkxrx[abs(h)][i]*snkyry[abs(k)][i]*snkzrz[abs(l)][i]*(h>0 ? 1:-1)*(k>0 ? 1:-1)*(l>0 ? 1:-1);
             }
		// std::cout<<"the real part is: "<<skre<<" the imaginary part is: "<<skim<<std::endl;
		 ksq=h*2*pi/p[0]*h*2*pi/p[0]+k*2*pi/p[1]*k*2*pi/p[1]+l*2*pi/p[2]*l*2*pi/p[2];
		 skmodsq=skre*skre+skim*skim;
		 LongRange=LongRange+1/volume/2*coul_prefactor*exp(-1*sigma*sigma*ksq/2)/ksq*skmodsq;
          }
		//std::cout<<"the long range energy is: "<<LongRange<<std::endl;
		std::cout<<"the toatl energy is: "<<std::setprecision(10)<<std::setw(10)<<LongRange+selfe+ShortRange<<std::endl;
		/*this is the final step*/
    for(size_t i=0;i<gmax;i++){
      delete [] cskxrx[i];
      delete [] cskyry[i];
      delete [] cskzrz[i];
      delete [] snkxrx[i];
      delete [] snkyry[i];
      delete [] snkzrz[i];
    }
    delete [] cskxrx;
    delete [] cskyry;
    delete [] cskzrz;
    delete [] snkxrx;
    delete [] snkyry;
    delete [] snkzrz;
    delete [] xall;
    delete [] yall;
    delete [] zall;
    delete [] allcharge;
    /*end memory allocation*/
}
