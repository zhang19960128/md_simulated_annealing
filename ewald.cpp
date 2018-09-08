#include "atom.h"
#include <iostream>
#include <math.h>
#include <Eigen/Dense>
#include <stdlib.h>
#include <complex>
#define pi 3.14159265359
using namespace Eigen;

/*Perform the Ewald summation. The lattice paramters should be given*/
/*Lattice should be a flattened array*/
void box::computelong(){
    double ShortRange = 0;
    double LongRange = 0;
    double selfe = 0;
    double *lattice = p;
    double epsil = 0; 
    atom *input = allatom;
    Matrix3d mtx_lattice;
    for (int i=0; i<9; i++){
        mtx_lattice << lattice[i];
    }
    double volume = mtx_lattice.determinant();
    Matrix3d gvecs = 2*pi*mtx_lattice.inverse(); /*By doing this, we actually got the transpose of the g lattice, which
    is however preferred due to following matrices multiplication*/
    
    /*Calculate the short-range energy*/
    Vector3d rijn;
    double deno = 0;
    double alpha = 1/pow(1/2,1/2)/sigma;/*Used to calculate the force*/
    for (int lx=-nmax; lx<nmax; lx++){
        for (int ly=-nmax; ly<nmax; ly++){
            for (int lz=-nmax; lz<nmax; lz++){
                for (size_t i=0; i<size; i++){
                    for (size_t j=0; j<size; j++){
                        if (i==j && lx==0 && ly==0 && lz==0){
                            continue;
                        }
                        else{
                            Vector3d vi;
                            Vector3d vj;
                            Vector3d nreal;
                            nreal << lx, ly, lz;
                            for (size_t k=0; k<3; k++){
                                vi << (input+i)->position[k];
                                vj << (input+j)->position[k];
                            }
                            rijn = (vi - vj + mtx_lattice.transpose()*nreal); 
                            deno = (rijn).norm();
                            epsil = coul[(input+i)->type][(input+j)->type];
                            ShortRange += 1/8/pi/epsil*(input+i)->charge*(input+j)->charge*erfc(deno/pow(2,0.5)/sigma)/deno;
                            for (size_t k=0; k<3; k++){
                                (input+i)->force[k] += 1/4/pi/epsil*(input+i)->charge*(input+j)->charge*rijn[k]/pow(deno,3)*(erfc(alpha*deno) + 2*alpha/pow(pi,1/2)*deno*exp(-pow(alpha*deno, 2))); 
                            } 
                        }
                    }
                }
            }
        }
    }
    
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
                    Vector3d ng;
                    ng << gx, gy, gz;
                    Vector3d gvec = gvecs*ng;
                    for (int i=0; i<size; i++){
                        std::complex<double> unit_i(0,1);
                        Vector3d vpos;
                        for (int j=0; j<3; j++){
                            vpos << (input+i)->position[j];
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
                                vi << (input+i)->position[kk];
                                vj << (input+j)->position[kk];
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
}
