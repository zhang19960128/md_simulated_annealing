#include "atom.h"
#include <stdio.h>
#include "image.h"
#include <new>
box(atom* inputallatom,
				int t,
				int s,
				double* period,
				double** pairbvco,
				double** cutoff
				){
	p=new double[3];
	for(size_t i=0;i<3;i++){
		p[i]=period[i];
	}
	allatom=new atom[s];
	type=t;
	size=s;
	pair_bv_co=new double* [t];
	for(size_t i=0;i<t;i++){
		pair_bv_co[i]=new double[t];
	}
	for(size_t i=0;i<t;i++)
		for(size_t j=0;j<t;j++){
			pair_bv_co[i][j]=pairbvco[i][j];
		}
	virtatom=imageall(allatom,size,period,)
}
