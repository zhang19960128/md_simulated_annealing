#include "atom.h"
#include <stdio.h>
#include "image.h"
#include <new>
#include <list>
#include <iostream>
#include <math.h>
void box::updatelistlj(){
    double tmp;
    double paircut;
    for (size_t i=0; i<size;i++){
        allatom[i].neilj.clear();
        for (size_t j=0; j<virtsize; j++){
            tmp = distance(allatom[i].position,virtatom[j].position);
            if (tmp<ljrcut && tmp>0.00001){
                allatom[i].neilj.push_back(j);
            }
        }
    }
}

void box::computelj(){
    double e_lj = 0.0;
    double delx, dely, delz,rsq,r;
    for (size_t i=0; i<size;i++){
        for (std::list<int>::iterator j=allatom[i].neilj.begin(); j!=allatom[i].neilj.end(); j++){
  			delx=allatom[i].position[0]-virtatom[*j].position[0];
			dely=allatom[i].position[1]-virtatom[*j].position[1];
			delz=allatom[i].position[2]-virtatom[*j].position[2];
			rsq=delx*delx+dely*dely+delz*delz;
			r=sqrt(rsq);
            ljenergy += pow(bij[allatom[i].type][allatom[*j].type]/r,12);
            allatom[i].force[0] += 12*pow(bij[allatom[i].type][allatom[*j].type],12)/pow(r,14)*delx;
            allatom[i].force[1] += 12*pow(bij[allatom[i].type][allatom[*j].type],12)/pow(r,14)*dely;
            allatom[i].force[2] += 12*pow(bij[allatom[i].type][allatom[*j].type],12)/pow(r,14)*delz;
        }
    }
    ljenergy = ljenergy/2;
}

