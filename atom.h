#ifndef atom_h
#define atom_h
/*since this is a light version of calculating the bond valence energy and force, we are not going to construct neighbor list for this atom*/
typedef struct Atom{
	double position[3];
	double force[3];
	double charge;//coulumb charge
	double bv0;//prefered bond valence_equilibrium.
	double bvreal;//real valence according to bond valence formula
	double bvv0;//prefered bond valence vector equilibrium
	double bvvreal;//real valence according to bond valence vector formula.
	int type;
	std::list<int> pairbv;
	std::list<int> pairbvv;
}atom;
class box{
	public:
		box()=default;//c++ 11 feature
		box(atom* inputallatom,
				int t,
				int s,
				double* period,
				double cutoff
				);
		void freezeforce();/*freeze force for other people to calculate accumulative force*/
		void updatelistforbv();/*update once and use forever, big trick*/
		void computebv(double** pair_bv_co,double rcut);
		~box(){
			delete p;
			delete allatom;
			delete virtatom;
		};
	private:
		int virtsize;
		double* p;//store the periodical boundary condition.
		atom* allatom;//store the atom array.
		int size;//store how many atoms are in the simulation box.
		int type;//store how many types of atoms are in the simulation box.
		atom* virtatom;//store the virtual atom image.
};
#endif
