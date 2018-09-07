#include "atom.h"
#include <stdio.h>
#include "image.h"
#include <new>
#include <list>
#include <iomanip>
#include <iostream>
#include <math.h>
void box::updatelistbvv(){
	double temp;
	double paircut;
	for(size_t i=0;i<size;i++){
		allatom[i].neibv.clear();
		for(size_t j=0;j<virtsize;j++){
			temp=distance(allatom[i].position,virtatom[j].position);
            paircut=bvvrcut[allatom[i].type][virtatom[j].type];
            if(temp<paircut && temp>0.0000001){
				allatom[i].neibv.push_back(j);
			}
		}
	}
}
void box::updatebvv(double** pairbvv_input){
	size_t temp=0;
	for(size_t i=0;i<type;i++)
		for(size_t j=i;j<type;j++){
			r0[i][j]=pairbvv_input[temp][0];
			cij[i][j]=pairbvv_input[temp][1];
			svvij[i][j]=pairbvv_input[temp][2];
			v0[i][j]=pairbvv_input[temp][3];
			bvvrcut[i][j]=pairbvv_input[temp][4];
			temp++;
		}
}
void box::computebvv(){
   /* bond valence vector parameters
    * typeone typetwo r0 Nij S V00 rcut
    *
    *
    *
    *
    *
    *
    */
double delx,dely,delz,rsq,recip,r,s,ss;
	bvvenergy=0.00;
	double* fp=new double[size];
	for(size_t i=0;i<size;i++){
		allatom[i].s0=0;
		allatom[i].s0x=0.0;
		allatom[i].s0y=0.0;
		allatom[i].s0z=0.0;
		for(std::list<int>::iterator j=allatom[i].neibv.begin();j!=allatom[i].neibv.end();j++){
			delx=allatom[i].position[0]-virtatom[*j].position[0];
			dely=allatom[i].position[1]-virtatom[*j].position[1];
			delz=allatom[i].position[2]-virtatom[*j].position[2];
			rsq=delx*delx+dely*dely+delz*delz;
			r=sqrt(rsq);
			recip=1.0/r;
			allatom[i].s0+=pow(r0[allatom[i].type][virtatom[*j].type]/r,cij[allatom[i].type][virtatom[*j].type]);
			allatom[i].s0x+=allatom[i].s0*delx/r;
			allatom[i].s0y+=allatom[i].s0*dely/r;
			allatom[i].s0z+=allatom[i].s0*delz/r;
		}
		s=allatom[i].s0x*allatom[i].s0x+allatom[i].s0y+allatom[i].s0y+allatom[i].s0z*allatom[i].s0z-vv0[allatom[i].type][allatom[i].type];
		ss=s*s;
		bvvenergy=bvvenergy+svvij[allatom[i].type][allatom[i].type]*ss;
	}
	std::cout<<"bond valence energy is: "<<bvvenergy<<std::endl;
}
