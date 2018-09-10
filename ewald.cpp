#include "atom.h"
#include <iostream>
#include <math.h>
#include <Eigen/Dense>
#include <stdlib.h>
#include <complex>
#include <list>
#define pi 3.14159265359
using namespace Eigen;
using namespace std;

/*Perform the Ewald summation. The lattice paramters should be given*/
/*Lattice should be a flattened array*/
void box::computelong(){
    double ShortRange = 0;
    double LongRange = 0;
    double selfe = 0;
    double *lattice = p;
    double epsil = 0;
    atom *input = allatom;
    MatrixXd mtx_lattice(3,3);
    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            if(i==j){
                mtx_lattice(i,j) = lattice[i];
            }
            else{
                mtx_lattice(i,j) = 0;
            }
        }
    }

    double volume = mtx_lattice.determinant();
    Matrix3d gvecs = 2*pi*mtx_lattice.inverse(); /*By doing this, we actually got the transpose of the g lattice, which
    is however preferred due to following matrices multiplication*/ 
    /*Calculate the short-range energy*/
    double delx, dely, delz,rsq,r;
    for (size_t i=0; i<size;i++){
        for (std::list<int>::iterator j=allatom[i].neilj.begin(); j!=allatom[i].neilj.end(); j++){ 
            delx=allatom[i].position[0]-virtatom[*j].position[0];
			dely=allatom[i].position[1]-virtatom[*j].position[1];
			delz=allatom[i].position[2]-virtatom[*j].position[2];
			rsq=delx*delx+dely*dely+delz*delz;
			r=sqrt(rsq);
            epsil = coul[allatom[i].type][allatom[*j].type]; 
            ShortRange += 1/pi/epsil/4*allatom[i].charge*virtatom[*j].charge/r;
            allatom[i].force[0] += -1/pi/epsil/4*allatom[i].charge*virtatom[*j].charge/(r*r*r)*delx; 
            allatom[i].force[1] += -1/pi/epsil/4*allatom[i].charge*virtatom[*j].charge/(r*r*r)*dely; 
            allatom[i].force[2] += -1/pi/epsil/4*allatom[i].charge*virtatom[*j].charge/(r*r*r)*delz; 
        }
    }
    ShortRange = ShortRange/2;



    /*Calculate the self energy*/
    for (int i=0; i<size; i++){
        selfe += 1/pow(2*pi,0.5)/sigma/4/pi/epsil*pow((input+i)->charge,2);
    }
    /*Calculate the long-range energy*/
    Vector3d rij;
    for (int gx=-gmax; gx<gmax; gx++){
        for (int gy=-gmax; gy<gmax; gy++){
            for (int gz=-gmax; gz<gmax; gz++){
                if (gx==0 && gy==0 && gz == 0){
                    continue;
                }
                else{
                    std::complex<double> sk(0,0);
                    double dgx =gx, dgy=gy, dgz=gz;
                    Vector3d ng;
                    ng << dgx, dgy, dgz;
                    Vector3d gvec;
                    gvec = gvecs*ng;
                    for (int i=0; i<size; i++){
                        std::complex<double> unit_i(0,1);
                        Vector3d vpos;
                        for (int j=0; j<3; j++){
                            vpos(j) = (input+i)->position[j];
                        }
                        sk += (input+i)->charge*exp(unit_i*gvec.dot(vpos));
                    }
                    double k = gvec.norm();
                    LongRange += 1/epsil/2*exp(-sigma*sigma*k*k/2)*pow(norm(sk),2)/volume/k/k;
                    /*Finish the calculation of energy, continue to calculate the force*/
                    for (size_t i=0; i<size; i++){
                        for (size_t j=0; j<size; j++){
                            Vector3d vi;
                            Vector3d vj;
                            for (size_t kk=0; kk<3; kk++){
                                vi(kk) = (input+i)->position[kk];
                                vj(kk) = (input+j)->position[kk];
                            }
                            rij = vi-vj;
                            for (size_t kk=0; kk<3; kk++){
                                (input+i)->force[kk] +=
                                -1/volume/epsil*(input+i)->charge*(input+j)->charge*sin(gvec.dot(rij))*exp(-sigma*sigma*k*k/2)*gvec[kk]/k/k; 
                            }
                        }
                    }
                }
            }
        }
    }
    
    coulenergy =  ShortRange - selfe + LongRange;    
    std::cout<<ShortRange<<"  "<<selfe<<"  "<<LongRange<<std::endl;
}
