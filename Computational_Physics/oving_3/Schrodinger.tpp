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
  double dx = 1.0/( eigenvectors.n_rows + 1);
  alpha_coeff.set_size( eigenvectors.n_rows );
  alpha_coeff.fill(0);
  arma::vec initial_condition ( eigenvectors.n_rows );
  arma::vec eigenvector ( eigenvectors.n_rows );
  eigenvector.fill(0);
  initial_condition.fill(0);
  //cout << eigenvectors.col(0) << endl;
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
