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
}atom;
class box{
	public:
		box()=default;//c++ 11 feature
		box(atom* inputallatom,
				int t,
				int s,
				double* period,
				double** pairbvco,
				double cutoff
				);
		~box(){
			delete p;
			delete allatom;
			delete virtatom;
			delete pair_bv_co;
		};
	private:
		int virtsize;
		double* p;//store the periodical boundary condition.
		atom* allatom;//store the atom array.
		int size;//store how many atoms are in the simulation box.
		int type;//store how many types of atoms are in the simulation box.
		atom* virtatom;//store the virtual atom image.
		double** pair_bv_co;//store the pair_bv_co;
}
#endif
