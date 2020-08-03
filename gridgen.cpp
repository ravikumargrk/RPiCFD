/*************************************************************
* Contents: Grid Generator and Solver
* Author  : G Ravi Kumar
* Dated   : 19-April-2019
*************************************************************/
// Job: Generates grid for given time step and append it to 
// gridfile.grd

// .grd file contents (each line):
// a_n    a_s    a_e    a_w    a_p    b  
// rows seperated by empty line and time frame sperated by 2 empty 
// lines.

#include <string>
#include <iostream>
#include "parameters.cpp"
#include "tdma.h"
#include <cmath>
using namespace std;

// Material Properties

// Dimenstions

double x(size_t i,size_t j) {return DX*(0.5 + double(j));}
double y(size_t i,size_t j) {return DY*(0.5 + double(i));}

double Diff_G(double x, double y){
    // Can this vary ? let it be a constant for now.
    return 0.45;
}

double rhoC(size_t i, size_t j){
    // Can this vary ? let it be a constant for now.
    return 2610*950;
}

double coeff(size_t dir,size_t i, size_t j){
    if(dir==0){ // North
        if(j==0) return (Diff_G(x(i,j),y(i,j)+DY/2))/(DY/2);
        else return (Diff_G(x(i,j),y(i,j)+DY/2))/(DY);
    }
    if(dir==1){ // South
        if(j==Nx-1) return (Diff_G(x(i,j),y(i,j)-DY/2))/(DY/2);
        else return (Diff_G(x(i,j),y(i,j)-DY/2))/(DY);
    }
    if(dir==2){ // EAST
        if(i==Ny-1) return (Diff_G(x(i,j)-DX/2,y(i,j)))/(DX/2);
        else return (Diff_G(x(i,j)-DX/2,y(i,j)))/(DX);
    }
    if(dir==3){ // WEST
        if(i==0) return (Diff_G(x(i,j)+DX/2,y(i,j)))/(DX/2);
        else return (Diff_G(x(i,j)+DX/2,y(i,j)))/(DX);
    }
}

bool isbound(size_t dir,size_t i, size_t j){
    if(dir==0){ // North
        if(j!=0) return false;
        else return true;
    }
    if(dir==1){ // South
        if(j!=Nx-1) return false;
        else return true;
    }
    if(dir==2){ // EAST
        if(i!=Ny-1) return false;
        else return true;
    }
    if(dir==3){ // WEST
        if(i!=0) return false;
        else return true;
    }
}

double getphi(double* matrix, size_t Nx, size_t Ny, size_t i, size_t j ) 
{return matrix[(j+1) + Nx*(i+1)];}


double source_term(size_t i, size_t j, double time_stamp){
    // Specifying source location here.
    // Note that this term is time varying !
    if((i >= S1X_begin/DX && i <= S1X_end/DX)&&((j >= S1Y_begin/DY && j <= S1Y_end/DY))){
        return 200000*cos(time_stamp);
    }
    if((i >= S2X_begin/DX && i <= S2X_end/DX)&&((j >= S2Y_begin/DY && j <= S2Y_end/DY))){
        return 200000*cos(time_stamp);
    }
    return 0;
}


// (i,j) is grid ID. the coeficients are assigned with respect to
// their grid IDs. all element lengths dxPE, dyPN.. etc are defined
// as functions of grid IDs 

void gridgen(double* grid_coeff, double* prev_phi, double time_stamp){
    for(size_t i = 0; i < Ny; i++){
        for(size_t j = 0; j < Nx; j++){
            for(size_t dir = 0; dir < 4 ; dir++){
                if(!isbound(dir,i,j)){
                    grid_coeff[(Nx*i + j)*6 + dir] = theta*coeff(dir,i,j);
                }
            }
            grid_coeff[(Nx*i + j)*6 + 4] = (2610*950*DX*DY/DT
                                            + (theta)*( coeff(0,i,j)+coeff(1,i,j)+coeff(2,i,j)+coeff(3,i,j) )
                                            );
            if(isbound(2,i,j)){
                grid_coeff[(Nx*i + j)*6 + 4] -= (coeff(2,i,j)/(1+(he/coeff(2,i,j))));
            }
            grid_coeff[(Nx*i + j)*6 + 5] = source_term(i, j, time_stamp)
                                         + (2610*950*DX*DY/DT
                                            - (1-theta)*( coeff(0,i,j)+coeff(1,i,j)+coeff(2,i,j)+coeff(3,i,j) )
                                            )*getphi(prev_phi,Nx,Ny,i,j);
                                         
            if(isbound(0,i,j)){
                grid_coeff[(Nx*i + j)*6 + 5] += coeff(0,i,j)*getphi(prev_phi,Nx,Ny,i,j+1);
            }
            else{
                grid_coeff[(Nx*i + j)*6 + 5] += (1-theta)*coeff(0,i,j)*getphi(prev_phi,Nx,Ny,i,j+1);
            }
            if(isbound(1,i,j)){
                grid_coeff[(Nx*i + j)*6 + 5] += coeff(1,i,j)*getphi(prev_phi,Nx,Ny,i,j-1);
            }
            else{
                grid_coeff[(Nx*i + j)*6 + 5] += (1-theta)*coeff(1,i,j)*getphi(prev_phi,Nx,Ny,i,j-1);
            }
            if(isbound(3,i,j)){
                grid_coeff[(Nx*i + j)*6 + 5] += coeff(3,i,j)*getphi(prev_phi,Nx,Ny,i-1,j);
            }
            else{
                grid_coeff[(Nx*i + j)*6 + 5] += (1-theta)*coeff(3,i,j)*getphi(prev_phi,Nx,Ny,i-1,j);
            }
            
            if(isbound(2,i,j)){
                grid_coeff[(Nx*i + j)*6 + 5] += (1-theta)*coeff(2,i,j)*getphi(prev_phi,Nx,Ny,i+1,j) + (Tinf/(1+(coeff(2,i,j)/he)));
            }
            else{
                grid_coeff[(Nx*i + j)*6 + 5] += (1-theta)*coeff(2,i,j)*getphi(prev_phi,Nx,Ny,i-1,j);
            }                                        
            // update convective boundary
            if(isbound(2,i,j)){
                prev_phi[(j+1) + Nx*(i+2)] =( he*Tinf + getphi(prev_phi,Nx,Ny,i,j)*coeff(2,i,j) )/ ( he + coeff(2,i,j) ); 
            }
        }
    }
    
    return ;
}

void solver(double* grid, double* guessphi){
    //k is same everywhere, but there are sources at two places.
    double* noboundguessphi = new double [Nx*Ny];
    for(size_t i = 0; i < Ny; i++){
        for(size_t j = 0; j < Nx; j++){
            noboundguessphi[(Nx*i + j)] = guessphi[(Nx*(i+1)+(j+1))];
        }
    }
    double* noboundnewphi = tdma2R(grid,noboundguessphi,Nx,Ny,1,false);
    for(size_t i = 0; i < Ny; i++){
        for(size_t j = 0; j < Nx; j++){
            guessphi[(Nx*(i+1)+(j+1))] = noboundnewphi[(Nx*i + j)];
        }
    }
    return ;
}