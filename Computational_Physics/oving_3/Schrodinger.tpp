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
  double dx = 1.0/( eigenvectors.n_cols);
  alpha_coeff.set_size( eigenvectors.n_cols );
  alpha_coeff.fill(0);
  for (int i = 0; i < eigenvectors.n_cols; i++){
    double x = (i+1)*dx;
    arma::vec eigenvector = eigenvectors.col(i);
    eigenvector *= condition(x);
    alpha_coeff(i) = trapezoidal(eigenvector);
<<<<<<< HEAD
    if(i==0){
      cout << eigenvector << endl;
      cout << condition(x) << endl;
    }
=======
    //cout << alpha_coeff(i) << endl;
>>>>>>> 101b471f55c9ec883116e581b6261792ddfbda1e
  }
}


#endif
