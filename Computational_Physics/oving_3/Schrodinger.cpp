#include "Schrodinger.hpp"
#include <armadillo>
#include <cstdlib>
#include <cassert>
#include <sstream>
#include <cmath>
#include <complex>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

//Include LAPACK
extern "C" void dstevd_(char* JOBZ, int* N, double* D, double* E, double* Z,int* LDZ,double* WORK,int* LWORK, int* IWORK, int* LIWORK,int* INFO );
extern "C" void zgttrf_(int* n, double* subdiag, double* diag, double* super_diag, double* DU2, int* IPIV, int* info);
extern "C" void zgttrs_(char* TRANS, int* N, int* NRHS, double* DL, double* D, double* DU, double* DU2, int* IPIV, double* B, int* LDB, int* INFO);


using namespace std;


//Constructor
Schrodinger::Schrodinger(int nx): IMUNIT (0.0, 1.0){
  diagonal.set_size(nx);
  sub_diagonal.set_size(nx-1);
  diagonal.fill(0);
  sub_diagonal.fill(0);
  dx = 1.0/(nx+1);
}

//Build sub diagonal
void Schrodinger::buildSubDiag(){
  for (int i = 0; i < sub_diagonal.n_elem; i++){
    sub_diagonal(i) = - 1.0/(dx*dx);
  }
}

//Step-by-step evolution by Crank-Nicolson
void Schrodinger::CrankNicolsonScheme(){
  int info;
  int IPIV[ u_last.n_elem ];
  int nt = 30000;
  int NRHS = 1;
  int n = u_last.n_elem;
  arma::cx_vec DU2( u_last.n_elem );
  arma::cx_vec u = u_last;
  arma::mat results(u_last.n_elem,0);
  results.fill(0);
  int counter = 0;

  //Symmetric tridiagonal system, sub equals super
  crank_super_diagonal_A = crank_sub_diag_A;

  //Cast complex vectors in order to use LAPACK
  double *super_diag = reinterpret_cast<double*>(crank_super_diagonal_A.memptr());
  double *subdiag = reinterpret_cast<double*>(crank_sub_diag_A.memptr());
  double *diag = reinterpret_cast<double*>(crank_diag_A.memptr());
  double *DU2_cast = reinterpret_cast<double*>(DU2.memptr());


  //Compute factorisation
  zgttrf_(&n, subdiag, diag, super_diag, DU2_cast, IPIV, &info);
  for(int i = 0; i<nt; i++){
    u_last = u;
    for(int j=0; j<u_last.n_elem-2; j++){
      //Compute the left hand side of the equation
      u(j+1) = crank_sub_diag_B(j)*u_last(j) + crank_diag_B(j+1)*u_last(j+1) + crank_sub_diag_B(j+1)*u_last(j+2);
    }
    u(0) = crank_diag_B(0)*u_last(0) + crank_sub_diag_B(0)*u_last(1);
    u(u.n_elem-1) = crank_diag_B(crank_diag_B.n_elem-1)*u_last(u.n_elem-1) + u_last(u.n_elem-2)*crank_sub_diag_B(crank_sub_diag_B.n_elem-2);

    //Solve the system Ax = Bb
    char TRANS[] = "N";
    double *u_cast = reinterpret_cast<double*>(u.memptr());
    zgttrs_(TRANS, &n, &NRHS, subdiag, diag, super_diag, DU2_cast, IPIV, u_cast, &n, &info);

    //Save to file
    if(i%100 == 0){
      arma::vec result = arma::pow(arma::abs(u),2);
      results.insert_cols(counter++, result);
    }

  }
  results.save("crankNicolsonScheme.csv", arma::csv_ascii);
}


//Compute alpha-coeffs
double Schrodinger::trapezoidal(const arma::vec &eigenvector){
  double integral = 0;
  for (int i = 1; i < eigenvector.n_elem; i++){
    integral += 2*eigenvector(i);
  }
  integral += eigenvector(0);
  integral += eigenvector(eigenvector.n_elem-1);
  return (double) integral/(2*eigenvector.n_elem);
}

//Get eigenvectors with time
arma::cx_vec& Schrodinger::getEigenvectorWithTime(arma::vec &eigenvector, double t, int eigenvalue){
  complex_eigenvector = arma::conv_to<arma::cx_vec>::from (eigenvector);
  complex_eigenvector *=exp(IMUNIT*diagonal(eigenvalue)*t);
  return complex_eigenvector;
}

//Get eigenvalues and eigenvectors
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

  //Eigensolver
  dstevd_(mode, &v, diagonal.memptr(), sub_diagonal.memptr(), eigenvectors.memptr(), &v, work.memptr(), &LWORK, iwork, &LIWORK, &info);
}

double Schrodinger::getMaxEigenvalue(){
  return arma::max(diagonal);
}

double Schrodinger::getMinEigenvalue(){
  return arma::min(diagonal);
}

//Normalization
void Schrodinger::normalizeEigenvectors(){
  for (int i = 0; i<eigenvectors.n_rows; i++){
    arma::vec eigenvector = arma::pow(eigenvectors.col(i),2);
    double n_factor = trapezoidal(eigenvector);
    eigenvectors.col(i) /= sqrt(n_factor);
  }
  //save(eigenvectors, diagonal);;
}

//Compute alpha-coeffs
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
}


//Save to file
void Schrodinger::save(arma::mat &eigenvectors, arma::vec &eigenvalues){
    stringstream fname;
    fname << "eigenvalues/eigenvales_with_potential_" << diagonal.n_elem << ".csv";
    eigenvalues.save(fname.str().c_str(), arma::csv_ascii);
    fname.str("");
    fname << "eigenvectors/eigenvector_with_potential_" << diagonal.n_elem << ".csv";
    eigenvectors.save(fname.str().c_str(), arma::csv_ascii);
}
