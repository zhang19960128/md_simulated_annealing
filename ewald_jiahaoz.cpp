#include "atom.h"
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <list>
#include <iomanip>
#define pi 3.14159265359
/*Perform the Ewald summation. The lattice paramters should be given*/
#define EWALD_F   1.12837917//actually this is 2/sqrt(3)
#define EWALD_P   0.3275911
#define A1        0.254829592
#define A2       -0.284496736
#define A3        1.421413741
#define A4       -1.453152027
#define A5        1.061405429
/*Lattice should be a flattened array*/
void box::computelong(){
    double ShortRange = 0;
    double LongRange = 0;
    double selfe = 0;
    double epsil = 0.0055263885;
		double coul_prefactor=180.9512801302711;
		double ewald_alpha=1/sqrt(2)/sigma;
		/*e/(epsilon0*A)=180.951  ev when use U=k*q1*q2/r*/
    double volume=p[0]*p[1]*p[2];
    double delx,dely,delz,rsq,r;
    double root2pi=sqrt(2*pi);
		double root2=sqrt(2);
		double rootpi=sqrt(pi);
    double* xall=new double[size];
    double* yall=new double [size];
    double* zall=new double [size];
    double* allcharge=new double [size];
		double* fx=new double [size];
		double* fy=new double [size];
		double* fz=new double [size];
		double chargei,chargej,temp,temp2,r3,erfc_interpolate,expm2,grij,t;//erfc_exact; use to debug when compare different erfc function.
    for(size_t i=0;i<size;i++){
       xall[i]=allatom[i].position[0];
       yall[i]=allatom[i].position[1];
       zall[i]=allatom[i].position[2];
			 fx[i]=0.00;
			 fy[i]=0.00;
			 fz[i]=0.00;
       allcharge[i]=allatom[i].charge;
			 chargei=allcharge[i];
       for(std::list<int>::iterator j=allatom[i].neilj.begin();j!=allatom[i].neilj.end();j++){
         chargej=virtatom[*j].charge;
				 delx=allatom[i].position[0]-virtatom[*j].position[0];
         dely=allatom[i].position[1]-virtatom[*j].position[1];
         delz=allatom[i].position[2]-virtatom[*j].position[2];
         rsq=delx*delx+dely*dely+delz*delz;
         r=sqrt(rsq);
				 r3=r*rsq;
				// temp=erfc(r/root2/sigma);
        // ShortRange+=1/epsil/4/pi*chargei*chargej/r*temp;
				 //temp2=1.0/8.0/epsil/pi*chargei*chargej/r3*erfc(r/root2/sigma); this is the standard erfc function but now I want to use the interpolation method to accelerate.
				 grij=ewald_alpha*r;
				 expm2=exp(-grij*grij);//=exp(-1*rsq/2/sigma^2)
				 t= 1.0 / (1.0 + EWALD_P*grij);
				 erfc_interpolate = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
				 //std::cout<<"the difference is: "<<erfc_interpolate-erfc(r/root2/sigma)<<std::endl;
				 ShortRange+=1/epsil/4/pi*chargei*chargej/r*erfc_interpolate;
				 temp2=1.0/8.0/epsil/pi*chargei*chargej/r3*erfc_interpolate;
				 fx[i]=delx*temp2+fx[i];
				 fy[i]=dely*temp2+fy[i];
				 fz[i]=delz*temp2+fz[i];
				 //temp=1.0/4/pi/epsil*1.0/2.0*chargei*chargej/rsq*(1.0/rootpi*2.0*expm2)*1.0/root2/sigma;
				 temp=1.0/4/pi/epsil*1.0/2.0*chargei*chargej/r3*1.0/rootpi*2.0*(EWALD_F*expm2*grij+erfc_interpolate);
				 fx[i]=fx[i]+temp*delx;
				 fy[i]=fy[i]+temp*dely;
				 fz[i]=fz[i]+temp*delz;
			 }
			 std::cout<<fx[i]<<" "<<fy[i]<<" "<<fz[i]<<std::endl;
    }
    ShortRange=ShortRange/2.0;
		//std::cout<<"the short range coulumb potential is: "<<std::setprecision(10)<<ShortRange<<std::endl;
    for(size_t i=0;i<size;i++){
      selfe=selfe-1/root2pi/sigma/4/pi/epsil*(allatom[i].charge)*allatom[i].charge;
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
    double skre,skim,skmodsq,ksq,snkr,cskr,kr;
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
	        skim=skim+allcharge[i]*cskxrx[abs(h)][i]*cskyry[abs(k)][i]*snkzrz[abs(l)][i]*(l>0 ? 1:-1);
	        skim=skim+allcharge[i]*snkxrx[abs(h)][i]*cskyry[abs(k)][i]*cskzrz[abs(l)][i]*(h>0 ? 1:-1);
	        skim=skim+allcharge[i]*cskxrx[abs(h)][i]*snkyry[abs(k)][i]*cskzrz[abs(l)][i]*(k>0 ? 1:-1);
	        skim=skim-allcharge[i]*snkxrx[abs(h)][i]*snkyry[abs(k)][i]*snkzrz[abs(l)][i]*(h>0 ? 1:-1)*(k>0 ? 1:-1)*(l>0 ? 1:-1);
						 }
          ksq=h*2*pi/p[0]*h*2*pi/p[0]+k*2*pi/p[1]*k*2*pi/p[1]+l*2*pi/p[2]*l*2*pi/p[2];
		      skmodsq=skre*skre+skim*skim;
		      temp=exp(-1*sigma*sigma*ksq/2)/ksq;
		      LongRange=LongRange+1/volume/2/epsil*temp*skmodsq;
						 for(size_t i=0;i<size;i++){
		//std::cout<<"the real part is: "<<skre<<" the imaginary part is: "<<skim<<std::endl;
		 			chargei=allcharge[i];
		 			kr=h*2*pi/p[0]*xall[i]+k*2*pi/p[1]*yall[i]+l*2*pi/p[2]*zall[i];
					snkr=sin(kr);
					//snkr=(snkxrx[abs(h)][i]*cskyry[abs(k)][i])*cskzrz[abs(l)][i]*(h>0?1:-1);
					//snkr=(snkxrx[abs(h)][i]*cskyry[abs(k)][i]*(h>0?1:-1)+cskxrx[abs(h)][i]*snkyry[abs(k)][i]*(k>0?1:-1))*cskzrz[abs(l)][i];
					//snkr=snkr
					/*
					 * the full formula require this, but due to symmetry, the net sum of snkzrz is basically zeros, so this is useless
					snkr=snkr+(cskxrx[abs(h)][i]*cskyry[abs(k)][i]-snkxrx[abs(h)][i]*snkyry[abs(k)][i]*(h>0?1:-1)*(k>0?1:-1))*snkzrz[abs(l)][i];
		 			*/
					cskr=cos(kr);
				  //	cskr=cskxrx[abs(h)][i]*cskyry[abs(k)][i]*cskzrz[abs(l)][i];
					//  cskr=(cskxrx[abs(h)][i]*cskyry[abs(k)][i]-snkxrx[abs(h)][i]*snkyry[abs(k)][i]*(h>0?1:-1)*(k>0?1:-1))*cskzrz[abs(l)][i];
					//	cskr=cskr-(snkxrx[abs(h)][i]*cskyry[abs(k)][i]*(h>0?1:-1)+cskxrx[abs(h)][i]*snkyry[abs(k)][i]*(k>0?1:-1))*snkzrz[abs(l)][i];
					fx[i]=fx[i]+1/volume/2/epsil*temp*chargei*2*pi/p[0]*h*2*(snkr*skre-cskr*skim);
					    fy[i]=fy[i]+1/volume/2/epsil*temp*chargei*2*pi/p[1]*k*2*(snkr*skre-cskr*skim);
					    fz[i]=fz[i]+1/volume/2/epsil*temp*chargei*2*pi/p[2]*l*2*(snkr*skre-cskr*skim);
						 }
						 }
		std::cout<<"the long range energy is: "<<std::setprecision(10)<<std::setw(10)<<LongRange<<std::endl;
		std::cout<<"the self energy is: "<<std::setprecision(10)<<std::setw(10)<<selfe<<std::endl;
		std::cout<<"the Short Range energy is: "<<std::setprecision(10)<<std::setw(10)<<ShortRange<<std::endl;
		std::cout<<"the toatl energy is: "<<std::setprecision(10)<<std::setw(10)<<LongRange+selfe+ShortRange<<std::endl;
		for(size_t i=0;i<size;i++){
			std::cout<<fx[i]<<" "<<fy[i]<<" "<<fz[i]<<std::endl;
		}
		for(size_t i=0;i<size;i++){
			allatom[i].force[0]+=fx[i];
			allatom[i].force[1]+=fy[i];
			allatom[i].force[2]+=fz[i];
		}
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
		delete [] fx;
		delete [] fy;
		delete [] fz;
    /*end memory allocation*/
}
