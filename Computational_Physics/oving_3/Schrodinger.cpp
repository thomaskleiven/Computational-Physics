#include "Schrodinger.hpp"
#include <armadillo>
#include <cstdlib>
#include <cassert>
#include <sstream>
#include <cmath>
#include <complex>

extern "C" void dstevd_(char* JOBZ, int* N, double* D, double* E, double* Z,int* LDZ,double* WORK,int* LWORK, int* IWORK, int* LIWORK,int* INFO );



using namespace std;

Schrodinger::Schrodinger(int nx):nx(nx), IMUNIT (0.0, 1.0){
  nx = nx;
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

double Schrodinger::trapezoidal(const arma::vec eigenvector){
  double integral = 0;
  for (int i = 1; i < eigenvector.n_elem; i++){
    integral += 2*eigenvector(i);
  }
  integral += eigenvector(0);
  integral += eigenvector(eigenvector.n_elem-1);
  return integral/(2*eigenvector.n_elem);
}

void Schrodinger::setEigenvectorWithTime(arma::vec &eigenvector, double t, int eigenvalue){
  complex_eigenvector = arma::conv_to<arma::cx_vec>::from (eigenvector);
  complex_eigenvector*=exp(IMUNIT*diagonal(eigenvalue)*t);
}


void Schrodinger::eigenvalueSolver(){

  char mode[] = "V";
  eigenvectors.set_size(nx,nx);
  arma::vec liwork( nx );
  int LIWORK = 3+5*nx;
  int LWORK = 1 + 4*nx + nx*nx;
  arma::vec work ( LWORK );
  int iwork[LIWORK];
  int info;

  dstevd_(mode, &nx, diagonal.memptr(), sub_diagonal.memptr(), eigenvectors.memptr(), &nx, work.memptr(), &LWORK, iwork, &LIWORK, &info);
  //save(eigenvectors, diagonal);
}

void Schrodinger::save(arma::mat &eigenvectors, arma::vec &eigenvalues){

    stringstream fname;
    fname << "eigenvalues/eigenvales_" << nx << ".csv";
    eigenvalues.save(fname.str().c_str(), arma::csv_ascii);
    eigenvectors = eigenvectors.t();
    fname.str("");
    fname << "eigenvectors/eigenvector_" << nx << ".csv";
    eigenvectors.save(fname.str().c_str(), arma::csv_ascii);
}


void Schrodinger::checkOrtogonality(){
    for (int i = 0; i<eigenvectors.n_cols; i++){
      for (int j = 0; j<eigenvectors.n_cols; j++){
        cout << arma::dot(eigenvectors.col(i), eigenvectors.col(j)) << endl;
      }
    }
}
