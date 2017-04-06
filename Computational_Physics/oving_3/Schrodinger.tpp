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
void Schrodinger::euler(const V& potential){
  setInitialCondition();
  arma::cx_vec psi ( last_psi.n_elem );
  arma::vec x = arma::linspace(0.0,1.0, last_psi.n_elem);
  unsigned int nt = 80;
  double dt = 1E-5;
  double del_x = 1.0 / x.n_elem;
  cdouble F (dt/(del_x*del_x),0.0);
  cdouble im (0.0,1);

  for ( int i = 0; i<nt; i++ ){
    for ( int j = 1; j<last_psi.n_elem-1; j++ ){
      psi(j) = last_psi(j) + im*F*(last_psi(j+1) - 2.0*last_psi(j) + last_psi(j-1)) + im*last_psi(j)*potential(x(j))*dt;
    }
    psi(0) = last_psi(0);
    psi(last_psi.n_elem-1) = last_psi(last_psi.n_elem-1);
    last_psi = psi;
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
    //alpha_coeff(i) = interpolation(eigenvector);
    //cout << alpha_coeff(i) << endl;
  }
}

#endif
