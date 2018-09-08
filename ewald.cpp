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
    atom *input = allatom;
    Matrix3d mtx_lattice;
    for (int i=0; i<9; i++){
        mtx_lattice << lattice[i];
    }
    double volume = mtx_lattice.determinant();
    Matrix3d gvecs = 2*pi*mtx_lattice.inverse(); /*By doing this, we actually got the transpose of the g lattice, which
    is however preferred due to following matrices multiplication*/
    
    /*Calculate the short-range energy*/
    for (int lx=-nmax; lx<nmax; lx++){
        for (int ly=-nmax; ly<nmax; ly++){
            for (int lz=-nmax; lz<nmax; lz++){
                for (int i=0; i<size; i++){
                    for (int j=0; j<size; j++){
                        if (i==j && lx==0 && ly==0 && lz==0){
                            continue;
                        }
                        else{
                            Vector3d vi;
                            Vector3d vj;
                            Vector3d nreal;
                            nreal << lx, ly, lz;
                            for (int k=0; k<3; k++){
                                vi << (input+i)->position[k];
                                vj << (input+j)->position[k];
                            }
                            double deno = (vi - vj + mtx_lattice.transpose()*nreal).norm();
                            ShortRange += 0.5*(input+i)->charge*(input+j)->charge*erfc(deno/pow(2,0.5)/sigma)/deno;
                        }
                    }
                }
            }
        }
    }
    
    /*Calculate the self energy*/
    for (int i=0; i<size; i++){
        selfe += 1/pow(2*pi,0.5)/sigma*pow((input+i)->charge,2);
    }
    
    /*Calculate the long-range energy*/
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
                    LongRange += 2*pi*exp(-sigma*sigma*k*k/2)*pow(norm(sk),2)/volume/k/k;
                }
            }
        }
    }
    
    coulenergy =  ShortRange - selfe + LongRange;    
}
