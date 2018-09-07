#include "atom.h"
#include <stdio.h>
#include "image.h"
#include <new>
#include <list>
#include <iomanip>
#include <iostream>
#include <math.h>
void box::updatelistbv(){
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
			sij[i][j]=pairbvv_input[temp][2];
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
}
