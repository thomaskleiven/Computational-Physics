#ifndef SCHRODINGER_TPP
#define SCHRODINGER_TPP

using namespace std;

template<class V>
void Schrodinger::initDiagonals(const V &potential){
  buildDiag(potential);
  buildSubDiag();
}


template<class V>
void Schrodinger::buildDiag( const V &potential ){
  for (int i = 0; i < nx; i++){
      diagonal(i) = 2.0/( dx*dx ) + potential( i*dx );
  }
}

template<class Function>
void Schrodinger::project( const Function &condition ){
  double dx = 1.0/( eigenvectors.n_cols+1 );
  alpha_coeff.set_size( eigenvectors.n_cols );
  alpha_coeff.fill(0);
  for (int i = 0; i < eigenvectors.n_cols; i++){
    double x = (i+1)*dx;
    arma::vec eigenvector = eigenvectors.col(i);
    eigenvector *= condition(x);
    alpha_coeff(i) = trapezoidal(eigenvector);
  }
}


#endif