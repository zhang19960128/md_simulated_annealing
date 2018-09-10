#include "atom.h"
#include <stdio.h>
#include "image.h"
#include <new>
#include <list>
#include <iomanip>
#include <iostream>
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
		double** pairbv_input,
		double** pairbvv_input,
        double** pairlj_input,
        double ljcut,
        double ewd_sigma,
        int ewd_nmax,
        int ewd_gmax
        ){
	/*this paribv_input should be similar to lammps input*/
	/*this paribvv_input should be similar to lammps input*/
	p=new double[3];
	for(size_t i=0;i<3;i++){
		p[i]=period[i];
	}
	allatom=new atom[s];
	std::copy(inputallatom,inputallatom+s,allatom);
	type=t;//specify how many type are in the simulation
	size=s;//specify how many atoms are in the simulation
	/*should have C(n,2) pair*/
	r0=new double* [t];/*equilibrium bond length*/
	v0=new double* [t];/*equilibrium valence zero only diagonal elements are considered*/
	cij=new double* [t];/*power law of bond valence*/
	sij=new double* [t];/*only contains diagonal elements*/
	svvij=new double* [t];/*only contains diagnoal elements*/
	bvrcut=new double* [t];/*cut-off for bond valence*/
	bvvrcut=new double* [t];/*cut-off for bond valence vector*/
	vv0=new double* [t];/*equlibrium bvv0*/
    coul = new double* [t];
    bij = new double* [t];
	for(size_t i=0;i<t;i++){
		r0[i]=new double[t];
		v0[i]=new double[t];
		cij[i]=new double[t];
		sij[i]=new double[t];
		svvij[i]=new double[t];
		bvrcut[i]=new double[t];
		bvvrcut[i]=new double[t];
		vv0[i]=new double[t];
        coul[i] = new double[t];
        bij[i] = new double[t];
	}
	size_t temp=0;
	double maxcutoff=0.0;
	for(size_t i=0;i<t;i++)
		for(size_t j=i;j<t;j++){
			r0[i][j]=pairbv_input[temp][0];
			r0[j][i]=r0[i][j];
			cij[i][j]=pairbv_input[temp][1];
			cij[j][i]=cij[i][j];
			sij[i][j]=pairbv_input[temp][2];
			sij[j][i]=sij[i][j];
			svvij[i][j]=pairbvv_input[temp][2];
			svvij[j][i]=svvij[i][j];
			v0[i][j]=pairbv_input[temp][3];
			v0[j][i]=v0[i][j];
			vv0[i][j]=pairbvv_input[temp][3];
			vv0[j][i]=vv0[i][j];
			bvrcut[i][j]=pairbv_input[temp][4];
			bvrcut[j][i]=bvrcut[i][j];
			bvvrcut[i][j]=pairbvv_input[temp][4];
			bvvrcut[j][i]=bvvrcut[i][j];
            coul[i][j] = pairlj_input[temp][0];
            coul[j][i] = coul[i][j];
            bij[i][j] = pairlj_input[temp][1];
            bij[j][i] = bij[i][j];
            maxcutoff=maxcutoff > bvrcut[i][j] ? maxcutoff : bvrcut[i][j];
			maxcutoff=maxcutoff > ljrcut ? maxcutoff : ljrcut;
            temp++;
		}
	int virt_size;
	virtatom=imageall(allatom,size,period,maxcutoff,virt_size);
	virtsize=virt_size;
    sigma = ewd_sigma;
    nmax = ewd_nmax;
    gmax = ewd_gmax;
    ljrcut = ljcut;
}
void box::freezeforce(){
	for(size_t i=0;i<size;i++){
		allatom[i].force[0]=0.0;
		allatom[i].force[1]=0.0;
		allatom[i].force[2]=0.0;
	}
}
void box::printnei(int i){
   for(std::list<int>::iterator a=allatom[i].neibv.begin();a!=allatom[i].neibv.end();a++){
      std::cout<<" "<<*a;
   }
}
/*finished computing bond valence force*/
/*end define the bond-valence energy*/
