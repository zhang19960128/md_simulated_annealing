#include "atom.h"
#include <new>
#include "image.h"
#include "stdio.h"
atom* imageall(atom* input,int size,double* p,double cutoff,int& virt_size){
	double min=0.0;
	for(size_t i=0;i<3;i++){
		if(p[i]>min)
			min=p[i];
	}
	int nimage=ceil(cutoff/min);
	virt_size=nimage*2+1;
	atom* allimage=new atom[virt*size*virt*virt];
	/*copy all the memory to different images*/
	for(size_t i=0;i<virt*virt*virt;i++){
		memcpy(allimage+i*size,input,sizeof(atom)*size);
	}
	int* shiftv[3]={0,0,0};
	int imagetick=0;/*this specify the image you need to shift*/
	for(int i=-1*virt;i<=virt;i++)
		for(int j=-1*virt;j<=virt;j++)
			for(int k=-1*virt;k<=virt;k++){
				shiftv[0]=i;
				shiftv[1]=j;
				shiftv[2]=k;
				imgetick=(i+virt)+(j+virt)*virt+(k+virt)*virt*virt;
				shift(allimage+imagetick,size,p,shiftv);
	}
	return allimage;
}
void shift(atom* input,int size,double* p,int* shiftv){
	for(int i=0;i<size;i++){
		for(size_t j=0;j<3;j++){
			(input+i)->position[j]=(input+i)->position[j]+shiftv[j]*p[i];
		}
	}
}
