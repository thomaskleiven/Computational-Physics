#include "Schrodinger.hpp"
#include <armadillo>
#include <cstdlib>
#include <cassert>

extern "C" void dpttrf_(int* n, double* diag, double* e, int* info);
extern "C" void dpttrs_(int* n, int* nrhs, double* diag, double* subdiag, double* rhs, int* nx, int* info);
extern "C" void dstevd_(char* JOBZ, int* N, double* D, double* E, double* Z,int* LDZ,double* WORK,int* LWORK, int* IWORK, int* LIWORK,int* INFO );



using namespace std;

Schrodinger::Schrodinger(int nx):nx(nx){
  diagonal.set_size(nx);
  sub_diagonal.set_size(nx-1);
  diagonal.fill(0);
  sub_diagonal.fill(0);
  dx = 1.0/nx;
}


void Schrodinger::buildSubDiag(){
  for (int i = 0; i < nx-1; i++){
    sub_diagonal(i) = - 1.0/(dx*dx);
  }
}

void Schrodinger::eigenvalueSolver(){

  char mode[] = "V";
  arma::mat orthonormal_eigenvectors(nx,nx);
  arma::vec liwork( nx );
  int LIWORK = 3+5*nx;
  int LWORK = 1 + 4*nx + nx*nx;
  arma::vec work ( LWORK );
  int iwork[LIWORK];
  int info;

  dstevd_(mode, &nx, diagonal.memptr(), sub_diagonal.memptr(), orthonormal_eigenvectors.memptr(), &nx, work.memptr(), &LWORK, iwork, &LIWORK, &info);

  diagonal.save("eigenvalues.csv", arma::csv_ascii);
}
