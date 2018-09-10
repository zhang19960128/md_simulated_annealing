#include "atom.h"
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <list>
#define pi 3.14159265359
/*Perform the Ewald summation. The lattice paramters should be given*/
/*Lattice should be a flattened array*/
void box::computelong(){
    double ShortRange = 0;
    double LongRange = 0;
    double selfe = 0;
    double epsil = 0;
    double volume=p[0]*p[1]*p[2];
    double delx,dely,delz,rsq,r;
    double root2pi=1/pow(2*pi,0.5);
    for(size_t i=0;i<size;i++){
      for(std::list<int>::iterator j=allatom[i].neilj.begin();j!=allatom[i].neilj.end();j++){
         delx=allatom[i].position[0]-virtatom[*j].position[0];
         dely=allatom[i].position[1]-virtatom[*j].position[1];
         delz=allatom[i].position[3]-virtatom[*j].position[2];
         rsq=delx*delx+dely*dely+delz*delz;
         r=sqrt(rsq);
         epsil=coul[allatom[i].type][allatom[*j].type];
         ShortRange+=1/pi/epsil/4*allatom[i].charge*virtatom[*j].charge/r;
      }
    }
    ShortRange=ShortRange/2.0;
    for(size_t i=0;i<size;i++){
      selfe=selfe+1/root2pi/sigma/4/pi/epsil*(allatom[i].charge)*allatom[i].charge;
    }
    /*calculate the longrange energy*/
    /*exp(i*k_v*r)
     *=exp(i*k_x*r_x)*exp(i*k_y*r_y)*exp(i*k_z*r_z)
     *=cos(k_x*r_x)*cos(k_y*r_y)*cos(k_z*r_z)-sin(k_x*r_x)*sin(k_y*r_y)*cos(k_z*r_z)-sin(k_x*r_x)*cos(k_y*r_y)*cos(k_z*r_z)-cos(k_x*r_x)*sin(k_y*r_y)*sin(k_z*r_z)
     *+I(********************************************************) similar term just expand it
     *
     * */

}
