#include "atom.h"
#include <stdio.h>
#include "image.h"
#include <new>
#include <math.h>
double distance(double* a,double* b){
	double s=0.0;
	for(size_t i=0;i<3;i++){
		s=s+(a[i]-b[i])*(a[i]-b[i]);
	}
	return sqrt(s);
}
box::box(atom* inputallatom,
		int t,
		int s,
		double* period,
		double cutoff
		){
	p=new double[3];
	for(size_t i=0;i<3;i++){
		p[i]=period[i];
	}
	allatom=new atom[s];
	type=t;
	size=s;
	int virt_size;
	virtatom=imageall(allatom,size,period,cutoff,virt_size);
	virtsize=virt_size;
}
void box::freezeforce(){
	for(size_t i=0;i<size;i++){
		allatom[i].force[0]=0.0;
		allatom[i].force[0]=0.0;
		allatom[i].force[0]=0.0;
	}
}
void box::updatelistbv(double rcut){
	double temp;
	for(size_t i=0;i<size;i++){
		allatom[i].neibv.clear();
		for(size_t j=0;j<virtsize;j++){
			temp=distance(allatom[i].position,virtatom[j].position);
			if(temp<rcut && temp>0.0000001){
				allatom[i].neibv.push_back(j);
			}
		}
	}
}
void box::updatelistbvv(double rcut){
	double temp;
	for(size_t i=0;i<size;i++){
		allatom[i].neibvv.clear();
		for(size_t j=0;j<virtsize;j++){
			temp=distance(allatom[i].position,virtatom[j].position);
			if(temp<rcut && temp>0.0000001){
				allatom[i].neibvv.push_back(j);
			}
		}
	}
}
/*starting a light version of bond valence with computing the force and bond valence energy*/
void box::computebv(double** pair_bv_co,double rcut){
	/*bond valence parameters.
	 *typeone typetwo r0 Nij S V0 rcut
	 * 1     1        #  #  # #  #   
	 * 1     1        #  #  # #  #
	 * 1     1        ###############
	 * 1     1        ###############
	 * 1     1        ###############
	 * 1     1        ###############
	 * 1     1        ###############
	 * 1     1        ###############
	 * 1     1        ###############
	 */

}
/*end define the bond-valence energy*/
