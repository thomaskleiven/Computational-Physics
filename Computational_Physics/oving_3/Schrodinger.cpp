#include "Schrodinger.hpp"
#include <armadillo>
#include <cstdlib>
#include <cassert>
#include <sstream>
#include <cmath>
#include <complex>

extern "C" void dstevd_(char* JOBZ, int* N, double* D, double* E, double* Z,int* LDZ,double* WORK,int* LWORK, int* IWORK, int* LIWORK,int* INFO );

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>


using namespace std;

Schrodinger::Schrodinger(int nx): IMUNIT (0.0, 1.0){
  diagonal.set_size(nx);
  sub_diagonal.set_size(nx-1);
  diagonal.fill(0);
  sub_diagonal.fill(0);
  dx = 1.0/(nx+1);
}


void Schrodinger::setInitialCondition(){
  arma::mat eigvecs;
  bool status = eigvecs.load("eigenvectors/eigenvector_with_potential_100.csv", arma::csv_ascii);
  if(status == true){ cout << "Eigenvectors loaded ok" << endl;}else{cout << "Could not load eigenvectors" << endl;}
  arma::vec eigenvector = eigvecs.col(0);
  last_psi = arma::conv_to<arma::cx_vec>::from (eigenvector);
}




void Schrodinger::buildSubDiag(){
  for (int i = 0; i < sub_diagonal.n_elem; i++){
    sub_diagonal(i) = - 1.0/(dx*dx);
  }
}

double Schrodinger::trapezoidal(const arma::vec &eigenvector){
  double integral = 0;
  for (int i = 1; i < eigenvector.n_elem; i++){
    integral += 2*eigenvector(i);
  }
  integral += eigenvector(0);
  integral += eigenvector(eigenvector.n_elem-1);
  return (double) integral/(2*eigenvector.n_elem);
}

arma::cx_vec& Schrodinger::getEigenvectorWithTime(arma::vec &eigenvector, double t, int eigenvalue){
  complex_eigenvector = arma::conv_to<arma::cx_vec>::from (eigenvector);
  complex_eigenvector *=exp(IMUNIT*diagonal(eigenvalue)*t);
  return complex_eigenvector;
}

void Schrodinger::eigenvalueSolver(){

  char mode[] = "V";
  eigenvectors.set_size(diagonal.n_elem, diagonal.n_elem);
  arma::vec liwork( diagonal.n_elem );
  int LIWORK = 3+5*diagonal.n_elem;
  int LWORK = 1 + 4*diagonal.n_elem + diagonal.n_elem*diagonal.n_elem;
  arma::vec work ( LWORK );
  int iwork[LIWORK];
  int info;
  int v = diagonal.n_elem;

  dstevd_(mode, &v, diagonal.memptr(), sub_diagonal.memptr(), eigenvectors.memptr(), &v, work.memptr(), &LWORK, iwork, &LIWORK, &info);
}

double Schrodinger::getMaxEigenvalue(){
  return arma::max(diagonal);
}

double Schrodinger::getMinEigenvalue(){
  return arma::min(diagonal);
}

void Schrodinger::normalizeEigenvectors(){
  for (int i = 0; i<eigenvectors.n_rows; i++){
    arma::vec eigenvector = arma::pow(eigenvectors.col(i),2);
    double n_factor = trapezoidal(eigenvector);
    eigenvectors.col(i) /= sqrt(n_factor);
  }
  //save(eigenvectors, diagonal);;
}

double Schrodinger::interpolation(arma:: vec &eigenvector){
  arma::vec x = arma::linspace(0.0, 1.0, eigenvector.n_elem);
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, eigenvector.n_elem);
  gsl_spline_init (spline, x.memptr(), eigenvector.memptr(), eigenvector.n_elem);
  double y{0};
  for(double i = 0; i<x.n_elem; i++){
    y = gsl_spline_eval(spline, x(i), acc);
  }
  return gsl_spline_eval_integ(spline, x(0), x(x.n_elem-1), acc);
  //gsl_spline_free(spline);
  //gsl_interp_accel_free(acc);


}





void Schrodinger::save(arma::mat &eigenvectors, arma::vec &eigenvalues){
    stringstream fname;
    fname << "eigenvalues/eigenvales_with_potential_" << diagonal.n_elem << ".csv";
    eigenvalues.save(fname.str().c_str(), arma::csv_ascii);
    fname.str("");
    fname << "eigenvectors/eigenvector_with_potential_" << diagonal.n_elem << ".csv";
    eigenvectors.save(fname.str().c_str(), arma::csv_ascii);
}
