#include "atom.h"
#include <stdio.h>
#include "image.h"
#include <new>
#include <list>
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
		double** pairbvv_input
		){
	/*this paribv_input should be similar to lammps input*/
	/*this paribvv_input should be similar to lammps input*/
	p=new double[3];
	for(size_t i=0;i<3;i++){
		p[i]=period[i];
	}
	allatom=new atom[s];
	type=t;//specify how many type are in the simulation
	size=s;//specify how many atoms are in the simulation
	/*should have C(n,2) pair*/
	r0=new double* [t];/*equilibrium bond length*/
	v0=new double* [t];/*equilibrium valence zero only diagonal elements are considered*/
	cij=new double* [t];/*power law of bond valence*/
	sij=new double* [t];/*only contains diagonal elements*/
	bvrcut=new double* [t];/*cut-off for bond valence*/
	for(size_t i=0;i<t;i++){
		r0[i]=new double[t];
		v0[i]=new double[t];
		cij[i]=new double[t];
		sij[i]=new double[t];
		bvrcut[i]=new double[t];
	}
	size_t temp=0;
	double maxcutoff=0.0;
	for(size_t i=0;i<t;i++)
		for(size_t j=i;j<t;j++){
			r0[i][j]=pairbv_input[temp][0];
			cij[i][j]=pairbv_input[temp][1];
			sij[i][j]=pairbv_input[temp][2];
			v0[i][j]=pairbv_input[temp][3];
			bvrcut[i][j]=pairbv_input[temp][4];
			maxcutoff=maxcutoff > bvrcut[i][j] ? maxcutoff : bvrcut[i][j];
			temp++;
		}
	int virt_size;
	virtatom=imageall(allatom,size,period,maxcutoff,virt_size);
	virtsize=virt_size;
}
void box::freezeforce(){
	for(size_t i=0;i<size;i++){
		allatom[i].force[0]=0.0;
		allatom[i].force[0]=0.0;
		allatom[i].force[0]=0.0;
	}
}
void box::updatelistbv(){
	double temp;
	double paircut;
	for(size_t i=0;i<size;i++){
		allatom[i].neibv.clear();
		for(size_t j=0;j<virtsize;j++){
			temp=distance(allatom[i].position,virtatom[j].position);
			paircut=bvrcut[allatom[i].type][virtatom[j].type];
			if(temp<paircut && temp>0.0000001){
				allatom[i].neibv.push_back(j);
			}
		}
	}
}
void box::updatelistbvv(){
	double temp;
	double paircut;
	for(size_t i=0;i<size;i++){
		allatom[i].neibvv.clear();
		for(size_t j=0;j<virtsize;j++){
			temp=distance(allatom[i].position,virtatom[j].position);
	    paircut=bvvrcut[allatom[i].type][virtatom[j].type];
			if(temp<paircut && temp>0.0000001){
				allatom[i].neibvv.push_back(j);
			}
		}
	}
}
/*starting a light version of bond valence with computing the force and bond valence energy*/
void box::computebv(){
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
	double delx,dely,delz,rsq,recip,r,s;
	bvenergy=0.00;
	for(size_t i=0;i<size;i++){
		allatom[i].s0=0;
		for(std::list<int>::iterator j=allatom[i].neibv.begin();j!=allatom[i].neibv.end();j++){
			delx=allatom[i].position[0]-virtatom[*j].position[0];
			dely=allatom[i].position[1]-virtatom[*j].position[1];
			delz=allatom[i].position[2]-virtatom[*j].position[2];
			rsq=delx*delx+dely*dely+delz*delz;
			r=sqrt(rsq);
			recip=1.0/r;
			allatom[i].s0+=pow(r0[allatom[i].type][virtatom[*j].type]/r,cij[allatom[i].type][virtatom[*j].type]);
		}
		s=allatom[i].s0-v0[allatom[i].type][allatom[i].type];
		bvenergy=sij[allatom[i].type][allatom[i].type]*(s*s)+bvenergy;
	}
	/*finished computing energy and started to compute force*/
	double Aij=0.0;
	for(size_t i=0;i<size;i++)
		for(std::list<int>::iterator j=allatom[i].neibv.begin();j!=allatom[i].neibv.end();j++){
			delx=allatom[i].position[0]-virtatom[*j].position[0];
			dely=allatom[i].position[1]-virtatom[*j].position[1];
			delz=allatom[i].position[2]-virtatom[*j].position[2];
			rsq=delx*delx+dely*dely+delz*delz;
			r=sqrt(rsq);
			recip=1.0/r;
			Aij=cij[allatom[i].type][virtatom[*j].type]*pow(r0[allatom[i].type][virtatom[*j].type]/r,cij[allatom[i].type][virtatom[*j].type])/r;
			allatom[i].force[0]+=2*sij[allatom[i].type][allatom[i].type]*(allatom[i].s0-v0[allatom[i].type][allatom[i].type])*Aij/r*delx;
			allatom[i].force[1]+=2*sij[allatom[i].type][allatom[i].type]*(allatom[i].s0-v0[allatom[i].type][allatom[i].type])*Aij/r*dely;
			allatom[i].force[2]+=2*sij[allatom[i].type][allatom[i].type]*(allatom[i].s0-v0[allatom[i].type][allatom[i].type])*Aij/r*delz;
	}
}
/*finished computing bond valence force*/
/*end define the bond-valence energy*/
