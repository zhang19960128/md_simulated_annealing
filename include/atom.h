#ifndef atom_h
#define atom_h
#include <list>
#include <iostream>
/*since this is a light version of calculating the bond valence energy and force, we are not going to construct neighbor list for this atom*/
typedef struct Atom{
	double position[3];
	double force[3];//MD force.
	double dftforce[3];//DFT force.
	double charge;//epsilonumb charge
	double bv0;//prefered bond valence_equilibrium.
	double s0;//real bond valence term
	double s0x;//bond valence vector x
	double s0y;//bond valence vector y
	double s0z;//bond valence vector z
	double bvreal;//real valence according to bond valence formula
	double bvv0;//prefered bond valence vector equilibrium
	double bvvreal;//real valence according to bond valence vector formula.
	int type;
	std::list<int> neibv;
	std::list<int> neibvv;
  std::list<int> neilj;
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
				double** pairbvv_input,
				double** pairlj_input,
                double ljcut=8.0
                );
		void init(atom* inputallatom,int s,double* period,double dft_energy,double** stress_dft,double w);
		void freezeforce();/*freeze force for other people to calculate accumulative force*/
		void updatelistbv();/*update once and use forever, big trick*/
		void updatelistbvv();
    void printlj();
		void updatebv(double** input);
    void updatebvv(double** input);
		void updatelistlj();
    void computebv();
		void computebvv();
		void computestress();
    void computelj();
    void computelong(double accuraccy);
		void printnei(int i);
		void printinfo();
		~box()
        {
			for(size_t i=0;i<type;i++){
				delete [] r0[i];
				delete [] v0[i];
				delete [] cij[i];
				delete [] sij[i];
                                delete [] svvij[i];
				delete [] bvrcut[i];
				delete [] bvvrcut[i];
				delete [] vv0[i];
                delete [] bij[i];
                delete [] epsilon[i];
			}
			delete [] r0;
			delete [] v0;
			delete [] cij;
			delete [] sij;
                        delete [] svvij;
			delete [] bvrcut;
			delete [] bvvrcut;
			delete [] vv0;
			delete [] allatom;
      delete [] bij;
      delete [] epsilon;
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
		double** svvij;//bvv energy coeffiecient
		double** bvrcut;//cut-off for bv
		double bvenergy;//energy produced by bond valence.
		double bvvenergy;//energy produced by bond valence energy.
		double ljenergy;//LJ 12 order energy
    double epsilonenergy;//Long range epsilonumb energy
    double** bvvrcut;//cut-off for bvv
		double** vv0;
    double ljrcut;
    double** bij;
    double** epsilon;
		double** stress;
		double** stressdft;
		double dftenergy;
		double weight;
};
#endif
