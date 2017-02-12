#include "HelmHoltz2D.hpp"
#include<iostream>
#include<sstream>


using namespace std;



HelmHoltz2D::HelmHoltz2D(unsigned int N_Discretization):N(N_Discretization){
  unsigned int M = 5*N*N - 4*N;
  locations.set_size(2, M);
  values.set_size(M);
}

HelmHoltz2D::~HelmHoltz2D(){
  if(matrix != NULL){
    delete matrix;
  }
}


void HelmHoltz2D::diagonal(){
  for (int i = 0; i < N*N; i++){
    locations(0,i) = i;
    locations(1,i) = i;
    values(i) = 4;
  }
}


void HelmHoltz2D::subDiagonal(){
  int size_subDiag = N*N;
  int counter = 1;
  for (int i = 1; i < size_subDiag; i++){
    int index_locations = N*N+counter-1;

    if(i%N == 0){
      continue;
    }

    locations(0,index_locations) = i;
    locations(1,index_locations) = i-1;
    values(index_locations) = -1;
    counter++;
  }
}

void HelmHoltz2D::superDiagonal(){
  int size_superDiag = N*N-1;
  int counter = 0;
  for (int i = 0; i< size_superDiag; i++){
    int index_locations = 2*N*N - N + counter;

    if((i+1)%N == 0){
      continue;
    }

    locations(0, index_locations) = i;
    locations(1, index_locations) = i+1;
    values(index_locations) = -1;
    counter++;
  }
}

void HelmHoltz2D::subsubDiag(){
  int size_subsub = N*N-N;
  for(int i = 0; i < size_subsub; i++){
    int index_locations = 3*N*N - 2*N +i;

    locations(0, index_locations) = i + N;
    locations(1, index_locations) = i;
    values(index_locations) = -1;
  }
}

void HelmHoltz2D::supsupDiag(){
  int size_supsup = N*N-N;
  for(int i = 0; i< size_supsup; i++){
    int index_locations = 4*N*N-3*N + i;

    locations(0, index_locations) = i;
    locations(1, index_locations) = i + N;
    values(index_locations) = -1;
  }
}

void HelmHoltz2D::initDiagonals(){
  diagonal();
  subDiagonal();
  superDiagonal();
  subsubDiag();
  supsupDiag();
}

void HelmHoltz2D::buildMatrix(){
  initDiagonals();
  if(matrix != NULL){delete matrix;}
  matrix = new arma::sp_mat(locations, values);
}

void HelmHoltz2D::solve(double numberOfEigenvalues){
  arma::eigs_sym(eigval, eigvec, *matrix, numberOfEigenvalues, "sm");
}

void HelmHoltz2D::save(){
  eigval *= (N*N);
  eigval.save("data/eigenvalues.csv", arma::csv_ascii);
  for (int i = 0; i < eigval.n_elem; i++){

    cout << eigvec.n_rows << endl;
    cout << eigvec.n_cols << endl;

    arma::mat A = eigvec.col(i);
    A.reshape(N,N);
    stringstream fname;
    fname << "data/eigenvector_" << static_cast<int>(i+1) << ".csv";
    A.save(fname.str().c_str(), arma::csv_ascii);
  }
}


////////////////////////////////////////////////////////////////////////////////
void HelmHoltz2DDebug::printMatrix(){
  buildMatrix();
  arma::mat denseMatrix(*matrix);
  cout << denseMatrix << endl;

  arma::eig_sym(eigval, eigvec, denseMatrix);
}

void HelmHoltz2DDebug::printLocations(){
  initDiagonals();
  cout << locations << endl;
}

void HelmHoltz2DDebug::checkOrtogonality(double numberOfEigenvalues){
  solve(numberOfEigenvalues);
  for (int i = 0; i<eigvec.n_cols;i++){
    for (int j = 0; j<eigvec.n_cols; j++){
      cout << i << " " << j << " " << arma::dot(eigvec.col(i), eigvec.col(j)) << endl;
    }
  }
}
