/*************************************************************
* Contents: Tri Diagonol Matrix Algorithm
* Author  : G Ravi Kumar
* Dated   : 19-Feb-2019
**************************************************************
* Refer tdma() for 1D and tdma2R() for 2D
*************************************************************/
#include "matrixtools.h"
using namespace std;

// Solves 1D N TDMA equations into solution array T.
// "matrix" (N by 4) has N rows of [ a_w    a_e    a_p    S_u].
double* tdma(double* matrix, size_t N){
    double* T = new double[N];
    double* A = new double[N];
    double* B = new double[N];
    // Forward Elimination
    A[0] = getval2d(matrix, 4, N, 0, 1)/ getval2d(matrix, 4, N, 0, 2);
    B[0] = getval2d(matrix, 4, N, 0, 3)/ getval2d(matrix, 4, N, 0, 2);
    for(size_t k = 1; k < N; k++ ){
        A[k] = getval2d(matrix, 4, N, k, 1)
               /(getval2d(matrix, 4, N, k, 2)-(getval2d(matrix, 4, N, k, 0)*A[k-1]));
        B[k] = ((getval2d(matrix, 4, N, k, 0)*B[k-1])+getval2d(matrix, 4, N, k, 3))
              /(getval2d(matrix, 4, N, k, 2)-(getval2d(matrix, 4, N, k, 0)*A[k-1]));
    }
    // Back substitution
    T[N-1] = B[N-1];
    for(size_t k=N-2;  ;k--){
        T[k]=(A[k]*T[k+1])+B[k];
        if(k==0) break;
    }
    // Free up memory
    delete[] A,B;
    return T;
}

double* trans(double* x, size_t Nx, size_t Ny){
    const size_t Nz = 6;
    double* y = new double[Nx*Ny*Nz];
    for(size_t i = 0;  i < Nx; i++){
        for(size_t j = 0;  j < Ny; j++){
            for(size_t k = 0; k < Nz; k++){   
                if(k==0 || k==1){
                    y[(Ny*i + j)*Nz + k] = getval3d(x,Nx,Ny,Nz,j,i,k+2);
                }
                else if(k==2 || k==3){
                    y[(Ny*i + j)*Nz + k] = getval3d(x,Nx,Ny,Nz,j,i,k-2);
                }
                else{
                    y[(Ny*i + j)*Nz + k] = getval3d(x,Nx,Ny,Nz,j,i,k);
                }
            }
        }
    }
    return y;
}

// Solves 2D TDMA with default sweep along +X
//"matrix" (Ny by Nx by 6) 3D pointer array 1D array of:
// [ a_n    a_s    a_e    a_w    a_p    b ]
// at each point matrix[Ny*j + Nx*i].
// default sweep direction and axis : +X
void tdma2d(double* matrix, double* T, size_t Nx, size_t Ny, bool change_sweep_axis){
    if(change_sweep_axis) {
        matrix = trans(matrix, Nx, Ny);
        size_t dum = Nx; Nx = Ny; Ny = dum;
    }
    const size_t Nz = 6;    
    double* Tx;
    double* coeff = new double[Ny*4] {0};

    for(size_t j = 0; j < Nx ; j++){         // Sweep loop [Y]
        for(size_t i = 0; i < Ny; i++){      // TDMA loop [X]
            coeff[i*4] = getval3d(matrix,Nx,Ny,Nz,i,j,0);    // a_n
            coeff[i*4+1] = getval3d(matrix,Nx,Ny,Nz,i,j,1);  // a_s
            coeff[i*4+2] = getval3d(matrix,Nx,Ny,Nz,i,j,4);  // a_p
            coeff[i*4+3] = getval3d(matrix,Nx,Ny,Nz,i,j,5) + // S + a_w*T_w + a_e*T_e
                            (getval3d(matrix,Nx,Ny,Nz,i,j,3)*(j==0?0:T[i*Nx+j-1])) + 
                            (getval3d(matrix,Nx,Ny,Nz,i,j,2)*(j==Nx-1?0:T[i*Nx+j+1]));
        }
        Tx = tdma(coeff, Ny);
        for(size_t p = 0; p < Ny; p++){
            T[p*Nx + j]=Tx[p];
        }
    }
    
    delete[] Tx;
    delete[] coeff;
    if(change_sweep_axis) T = trans2d(T, Ny, Nx);
    return ;
}

double* tdma2R(double* matrix, double* Ti, size_t Nx, size_t Ny, double err, bool change_sweep_axis){
    double* Ti_old = new double[Nx*Ny];
    copymatrix(Ti, Ti_old, Nx*Ny);
    tdma2d(matrix, Ti, Nx, Ny, change_sweep_axis);

    size_t iter = 1;
    while(sumofdiff(Ti, Ti_old, Nx*Ny) >= err){
        iter++;
        copymatrix(Ti, Ti_old, Nx*Ny);
        tdma2d(matrix, Ti, Nx, Ny, change_sweep_axis);
    }
    cout<<"\nSystem converged in "<<iter<<" Iterations.\n\n";
    return Ti;
}