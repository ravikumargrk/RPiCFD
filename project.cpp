/*************************************************************
* Contents: Assembler
* Author  : G Ravi Kumar
* Dated   : 19-April-2019
*************************************************************/
// Job: calls gridgen and solver in correct order, 
// saves grid and solution data.

#include <fstream>
#include "gridgen.cpp"

using namespace std;
void writephi(double* phi){
    ofstream solution;
    solution.open("result",ios::app);
    for(size_t i = 0; i < Ny; i++){
        for(size_t j = 0; j < Nx; j++){
            solution<<y(i,j)<<"\t"<<x(i,j)<<"\t"<<phi[(Nx*(i+1) + (j+1))]<<"\n";
        }
        solution<<"\n";
    }
    solution<<"\n";
    solution.close();
}

void writegrid(double* grid){
    ofstream gridfile;
    gridfile.open("grid",ios::app);
    for(size_t i = 0; i < Ny; i++){
        for(size_t j = 0; j < Nx; j++){
            gridfile<<x(i,j)<<"\t"<<y(i,j)<<"\t";
            for(size_t k = 0; k < 6; k++){
                gridfile<<grid[(Nx*i + j)*6 + k]<<"\t";
            }
            if(grid[(Nx*i + j)*6 + 4] > grid[(Nx*i + j)*6 + 0] + grid[(Nx*i + j)*6 + 1] + grid[(Nx*i + j)*6 + 2] + grid[(Nx*i + j)*6 + 3] ){
                gridfile<<"True V";
            }else {gridfile<<"False X";}
            gridfile<<"\n";
        }
        gridfile<<"\n";
    }
    gridfile<<"\n";
    gridfile.close();
}

int main(){
    remove("grid");
    remove("result");
    // Initialise grid and guess solution
    double* phi = new double [(Nx+2)*(Ny+2)] {0};
    for(size_t i = 0; i < (Nx+2)*(Ny+2); i++){
        phi[i] = 300;
    }
    // Iterate
        // Call gridgen
    double* grid = new double[Nx*Ny*6] {0};
    for(size_t t = 0; t<2; t++){
        gridgen(grid,phi,double(t)*DT); // grid updates
        // Write grid to file
        writegrid(grid);
        solver(grid,phi);
    }
    delete[] grid;
    writephi(phi);
    delete[] phi;
    system("gnuplot \"plotscript.plt\"");
}

