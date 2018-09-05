#ifndef atom_h
#define atom_h
#include <list>
#include <iostream>
/*since this is a light version of calculating the bond valence energy and force, we are not going to construct neighbor list for this atom*/
typedef struct Atom{
	double position[3];
	double force[3];
	double charge;//coulumb charge
	double bv0;//prefered bond valence_equilibrium.
	double s0;//real bond valence term
	double bvreal;//real valence according to bond valence formula
	double bvv0;//prefered bond valence vector equilibrium
	double bvvreal;//real valence according to bond valence vector formula.
	int type;
	std::list<int> neibv;
	std::list<int> neibvv;
}atom;
double distance(double* a,double* b);
class box{
	public:
		box()=default;//c++ 11 feature
		box(atom* inputallatom,
				int t,//type numbers
				int s,//total number of atoms
				double* period,
				double** pairbv_input,
				double** pairbvv_input
				);
		void freezeforce();/*freeze force for other people to calculate accumulative force*/
		void updatelistbv();/*update once and use forever, big trick*/
		void updatelistbvv();
		void updatebv(double** input);
		void computebv();
		void printnei();
		~box(){
			for(size_t i=0;i<type;i++){
				delete [] r0[i];
				delete [] v0[i];
				delete [] cij[i];
				delete [] sij[i];
				delete [] bvrcut[i];
				delete [] bvvrcut[i];
				delete [] vv0[i];
			}
			delete [] r0;
			delete [] v0;
			delete [] cij;
			delete [] sij;
			delete [] bvrcut;
			delete [] bvvrcut;
			delete [] vv0;
			delete [] allatom;
		};
	private:
		int virtsize;//store how many image atoms are there.
		double* p;//store the periodical boundary condition.
		atom* allatom;//store the atom array.
		int size;//store how many atoms are in the simulation box.
		int type;//store how many types of atoms are in the simulation box.
		atom* virtatom;//store the virtual atom image.
	  /*the name of the parameters refer to Shi Liu
		 *Reinterpretation of the bond-valence model with bond-valence model with bond-order formalism: An improved bond-valence-based interatomic potential for PbTiO3
		 * */
		double** r0;//bv r0
		double** v0;//bv v0
		double** cij;//bv powerlaw
		double** sij;//bv energy coeffiecient
		double** bvrcut;//cut-off for bv
		double bvenergy;//energy produced by bond valence.
		double** bvvrcut;//cut-off for bvv
		double** vv0;
};
#endif
