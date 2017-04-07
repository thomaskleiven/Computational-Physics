#ifndef SCHRODINGER_TPP
#define SCHRODINGER_TPP

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

using namespace std;

template<class V>
void Schrodinger::initDiagonals(const V &potential){
  buildDiag(potential);
  buildSubDiag();
}


template<class V>
void Schrodinger::buildDiag( const V &potential ){
  for (int i = 0; i < diagonal.n_elem; i++){
      diagonal(i) = 2.0/( dx*dx ) + potential( i*dx );
  }
}

template<class V>
void Schrodinger::setDiagForCrank( const V& potential ){
  setInitialCondition();
  cdouble im (0.0,1.0);
  crank_diag_A.set_size( u_last.n_elem );
  crank_diag_B.set_size( u_last.n_elem );
  crank_sub_diag_A.set_size( u_last.n_elem - 1 );
  crank_sub_diag_B.set_size( u_last.n_elem - 1 );
  double dt = 1E-6;
  double del_x = 1.0/ u_last.n_elem;

  for( int i=0; i<u_last.n_elem; i++ ){
    crank_diag_A (i) =  1.0 + im*(dt/2.0) * (-2.0/(del_x*del_x) + potential(del_x*(i)));
    crank_diag_B (i) =  1.0 - im*(dt/2.0) * (-2.0/(del_x*del_x) + potential(del_x*(i)));
    if(i>0){
      crank_sub_diag_A (i-1) = -im*(dt/2.0) * 1.0/(del_x*del_x);
      crank_sub_diag_B (i-1) = im*(dt/2.0) * 1.0/(del_x*del_x);
    }
  }
}




template<class V>
void Schrodinger::euler(const V& potential){
  //setInitialCondition();
  arma::cx_vec psi ( u_last.n_elem );
  arma::vec x = arma::linspace(0.0,1.0, u_last.n_elem);
  unsigned int nt = 80;
  double dt = 1E-5;
  double del_x = 1.0 / x.n_elem;
  cdouble F (dt/(del_x*del_x),0.0);
  cdouble im (0.0,1.0);

  for ( int i = 0; i<nt; i++ ){
    for ( int j = 1; j<u_last.n_elem-1; j++ ){
      psi(j) = u_last(j) + im*F*(u_last(j+1) - 2.0*u_last(j) + u_last(j-1)) + im*u_last(j)*potential(x(j))*dt;
    }
    psi(0) = u_last(0);
    psi(u_last.n_elem-1) = u_last(u_last.n_elem-1);
    u_last = psi;
  }
  arma::vec result = arma::pow(arma::abs(psi),2);
  result.save("euler_scheme.csv", arma::csv_ascii);
}




template<class Function>
void Schrodinger::project( const Function &condition ){
  alpha_coeff.set_size( eigenvectors.n_rows );
  alpha_coeff.fill(0);
  arma::vec eigenvector ( eigenvectors.n_rows );
  eigenvector.fill(0);
  for (int i = 0; i < eigenvectors.n_cols; i++){
    eigenvector = eigenvectors.col(i);
    for (int j = 0; j < eigenvectors.n_rows; j++ ){
      double x = (j+1)*dx;
      eigenvector(j) *= condition(x);
    }
    alpha_coeff(i) = trapezoidal(eigenvector);
  }
}

#endif
