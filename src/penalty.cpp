#include "penalty.h"
#include "atom.h"
#include "readpara.h"
/*
 * the main AIM of the two function is to map the list variable into bvv model parameters and charge parameters. 
 *
 * /
/*tell you which index you are going to map tick to (i,j) */
void indexbvvmap(int** bvvmatrixmap,int tick,int& i,int& j){
	size_t count=0;
	for(size_t m=0;m<control::pair_num;m++)
		for(size_t n=0;n<12;n++){
			if(bvvmatrixmap[m][n]==1){
				if(count==tick){
					i=m;
					j=n;
				}
				count++;
			}
		}
}
/*this tell you which index you are going to map tick to (i) in chargemap
*now we are dealing with only asite, bsite and osite charge change.
tick is the index of xp in the map function
*/
void indexchargemap(int* chargemap,int tick,int& i){
	int count=tick+1;
	int remain=count-control::paracount_bvv;/*the remaining parameters to be optimized, Assumed to be in the charge map*/
	size_t sum=0;
	for(size_t m=0;m<species::spe.size();m++){
		if(chargemap[m]==1){
			sum=sum+1;
			if(sum==remain){
				i=m;
			}
		}
	}
}

void map(double* xp){
	int m,n;
 	for(int i=0;i<control::paracount_bvv;i++){
		indexbvvmap(control::bvvmatrixmap,i,m,n);
		control::bvvmatrix[m][n]=xp[i];
	}
//    double accCharge = 0;
	for(int i=control::paracount_bvv;i<control::paracount_bvv+control::paracount_charge;i++){
		indexchargemap(control::chargemap,i,m);
        control::charge[m]=xp[i];
//        if (control::site[m] == 0 || control::site[m] == 1) accCharge += charge[m];
//        else accCharge += charge[m]*3;
    }
	/*perseve charge conservation law! Zhenbang*/
    for (int i=0; i<species::spe.size(); i++){
        if (chargemap[i] == 0){
            if (species::site[i] == 0 || species::site[i] == 1);
            control::charge[i] = (0-accCharge);

        }
    }
}

//Zhenbang
void mapToXp(double* xp){
    count = 0;
    for (size_t m=0; m<control::pair_num; m++){
        for (size_t n=0; n<12; n++){
            if (control::bvvmatrixmap[m][n] == 0) continue;
            xp[count] = control::bvvmatrix[m][n];
            count += 1;
        }
    }

    for (size_t i=0; i<species::spe.size(); i++){
        if (control::chargemap[i] == 1){
            xp[count] = control::charge[i];//the number of this charge array is Nspecies-1 
            count += 1;
        }
    }
}

int referenceStruct(box* system, int systemSize){
    int index = -1;
    reference = 1e20;
    for (size_t i=0; i<systemSize; i++){
        if (system[i].dftenergy < reference){
           index = i; 
        }
    }
    return i;
}

//Zhenbang
double PenaltyFunc(double* xp, box* system){
    int indexRef = saconst::sa_refindex;
    map(xp);
    double penalty = 0.0;
    box* ionall = system; 

    double* totalEnergy = new double[number];

    for (size_t i=0; i<number; i++){
        totalEnergy[i] = 0;
        ionall[i].updatebvparameter(control::bvmatrix);
        for (size_t j=0; j<ionall[i].size; j++){
            for (size_t k=0; k<species::spe.size()){
                if (ionall[i].allatom[j].type == k){
                    ionall[i].allatom[j].charge == control::charge[k];
                }
            }
        }
        ionall[i].computeAll();
        totalEnergy[i] = ionall[i].bvenergy + ionall[i].bvvenergy + ionall[i].ljenergy + ionall[i].epsilonenergy;
        
    /*Calculate the penalty*/
    double PenaltyE = 0;
    double PenaltyF = 0;
    for (size_t i=0; i<number; i++){
        PenaltyE += fabs((totalEnergy[i]-totalEnergy[indexRef]) - (ionall[i].dftenergy-ionall[indexRef].dftenergy))*ionall[i].weight;
        for (size_t j=0; j<ionall[i].size; j++){
            for (size_t k=0; k<3; k++){
                PenaltyF += fabs((ionall[i].allatom[j].force[k]-ionall[indexRef].allatom[j].force[k]) - (ionall[i].allatom[j].dftforce[k]-ionall[indexRef].allatom[j].dftforce[k]))*ionall[i].weight; 
            }
        }
    }
    PenaltyE = PenaltyE/number*saconst::sa_eweight;
    PenaltyF = PenaltyF/(number*3*ionall[0].size)*saconst::sa_fweight;
    delete[] totalEnergy;

    return penalty = PenaltyE + PenaltyF;
}
