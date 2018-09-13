#include "atom.h"
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <complex>
#include <iomanip>
#include <list>
#define pi 3.14159265359
#define sqrt2 1.41421356237
#define sqrtpi 1.77245385091

/*Perform the Ewald summation. The lattice paramters should be given*/
/*Lattice should be a flattened array*/
void box::computelong(){
    double ShortRange = 0;
    double LongRange = 0;
    double selfe = 0;
    double *lattice = p;
    double epsil = 0.0055263885;
    double gvecs[3]={0};
		atom *input = allatom;
    double volume = lattice[0]*lattice[1]*lattice[2];
		std::cout<<volume<<std::endl;
		gvecs[0] = 2*pi/lattice[0]; gvecs[1] = 2*pi/lattice[1]; gvecs[2] = 2*pi/lattice[2];
		/*By doing this, we actually got the transpose of the g lattice, which
    is however preferred due to following matrices multiplication*/
    /*Calculate the short-range energy*/
    double delx, dely, delz,rsq,r;
    for (size_t i=0; i<size;i++){
        for (std::list<int>::iterator j=allatom[i].neilj.begin(); j!=allatom[i].neilj.end(); j++){
            delx=allatom[i].position[0]-virtatom[*j].position[0];
						dely=allatom[i].position[1]-virtatom[*j].position[1];
			      delz=allatom[i].position[2]-virtatom[*j].position[2];
			      rsq=delx*delx+dely*dely+delz*delz;
			      r=sqrt(rsq);
            ShortRange += 1/pi/epsil/4*allatom[i].charge*virtatom[*j].charge/r*erfc(r/sqrt2/sigma);
            allatom[i].force[0] += -1/pi/epsil/4*allatom[i].charge*virtatom[*j].charge/(r*r*r)*delx*erfc(r/sqrt(2)/sigma) - 1/pi/epsil/4*allatom[i].charge*virtatom[*j].charge/r/r*exp(-r*r/2/sigma/sigma)/sqrt(pi*2)/sigma*delx;
            allatom[i].force[1] += -1/pi/epsil/4*allatom[i].charge*virtatom[*j].charge/(r*r*r)*dely*erfc(r/sqrt(2)/sigma) - 1/pi/epsil/4*allatom[i].charge*virtatom[*j].charge/r/r*exp(-r*r/2/sigma/sigma)/sqrt(pi*2)/sigma*dely;	
            allatom[i].force[2] += -1/pi/epsil/4*allatom[i].charge*virtatom[*j].charge/(r*r*r)*delz*erfc(r/sqrt(2)/sigma) - 1/pi/epsil/4*allatom[i].charge*virtatom[*j].charge/r/r*exp(-r*r/2/sigma/sigma)/sqrt(pi*2)/sigma*delz;	
        }
    }
    ShortRange = ShortRange/2;



    /*Calculate the self energy*/
    for (int i=0; i<size; i++){
        selfe += 1/(sqrt2*sqrtpi)/sigma/4/pi/epsil*(input+i)->charge*(input+i)->charge;
    }
		double ri[3] = {0};
    double rij[3] = {0}; 
    for (int gx=-gmax; gx<=gmax; gx++){
        for (int gy=-gmax; gy<=gmax; gy++){
            for (int gz=-gmax; gz<=gmax; gz++){
                if (gx==0 && gy==0 && gz == 0){
                    continue;
                }
                else{
                    std::complex<double> sk(0,0);
										double dgx =gx, dgy=gy, dgz=gz;
										double gvec[3] = {gvecs[0]*dgx, gvecs[1]*dgy, gvecs[2]*dgz};
                    for (int i=0; i<size; i++){
                        std::complex<double> unit_i(0,1);
                        for (int j=0; j<3; j++){
                            ri[j] = (input+i)->position[j];
                        }
                        sk += (input+i)->charge*exp(unit_i*(gvec[0]*ri[0] + gvec[1]*ri[1] + gvec[2]*ri[2]));
                    }
										//std::cout<<"the sk is"<<sk<<std::endl;

                    double k = pow(gvec[0]*gvec[0]+gvec[1]*gvec[1]+gvec[2]*gvec[2],0.5);	
										LongRange += 1/epsil/2*exp(-sigma*sigma*k*k/2)*norm(sk)/volume/k/k;
										//std::cout<<norm(sk)*norm(sk)<<"  "<<norm(sk)<<"  "<<sk<<std::endl;
										//std::cout<<volume<<"  "<<sigma<<"  "<<norm(sk)*norm(sk)<<"  "<<k*k<<std::endl;
                    /*Finish the calculation of energy, continue to calculate the force*/
                    for (size_t i=0; i<size; i++){
                        for (size_t j=0; j<size; j++){ 
                            for (size_t kk=0; kk<3; kk++){
                                rij[kk] = (input+i)->position[kk] - (input+j)->position[kk];
                            } 
                            for (size_t kk=0; kk<3; kk++){
                                (input+i)->force[kk] +=
                                1/volume/epsil*(input+i)->charge*(input+j)->charge*sin(gvec[0]*rij[0] + gvec[1]*rij[1] + gvec[2]*rij[2])*exp(-sigma*sigma*k*k/2)*gvec[kk]/k/k;
                            }
                        }
                    }
                }
            }
        }
    }

    epsilonenergy =  ShortRange - selfe + LongRange;   
		std::cout<<std::setprecision(10)<<"the shortrange is: "<<ShortRange<<" the long range is: "<<LongRange<<" the self energy is: "<<selfe<<std::endl;
		std::cout<<"the total energy is: "<<std::setprecision(10)<<epsilonenergy<<std::endl;
		/*for (size_t i=0; i<size; i++){
			std::cout<<epsilonenergy<<"  "<<(input+i)->force[0]<<"  "<<(input+i)->force[1]<<"  "<<(input+i)->force[2]<<std::endl;
		}
		*/
}
