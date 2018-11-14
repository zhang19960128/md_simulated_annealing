#include "penalty.h"
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
void indexchargemap(int* chargemap,int tick,int& i){
	int count=tick+1;
	int remain=count-control::paracount_bvv;/*the remaining parameters to be optimized, Assumed to be in the charge map*/
	size_t sum=0;
	for(size_t m=0;m<species::num.size();m++){
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
	for(int i=control::paracount_bvv;i<control::paracount_bvv+control::paracount_charge;i++){
		indexchargemap(control::chargemap,i,m);
		control::charge[m]=xp[i];
	}
	/*perseve charge conservation law!*/
}
double PenaltyFunc(double* ){
}
