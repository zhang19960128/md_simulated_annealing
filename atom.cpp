#include "atom.h"
#include <stdio.h>
#include "image.h"
#include <new>
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
box::freezeforce(){
	for(size_t i=0;i<size;i++){
		allatom[i].force[0]=0.0;
		allatom[i].force[0]=0.0;
		allatom[i].force[0]=0.0;
	}
}
box::updatelistbv(){

}
/*starting a light version of bond valence with computing the force and bond valence energy*/
void box::computebv(double** pair_bv_co,double rcut){
	/*bond valence parameters.
	 *typeone typetwo r0 Nij S V0 rcut
	 * 1     1        ###############   
	 * 1     1        ###############
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
