#include "TimeEvolution.hpp"
#include<iostream>
#include<sstream>
#include<armadillo>
#include<cassert>
#include<iostream>
#include<cmath>


using namespace std;

void TimeEvolution::loadMatrix(const string &filename){
  matrix.load(filename);
  matrix *= matrix.n_rows;
  //cout << matrix << endl;

}

void TimeEvolution::loadEigenvalue(const string &filename, unsigned int modenumber){
  assert(modenumber >= 1);
  arma::vec eigenvlues;
  eigenvlues.load(filename);
  eigenvalue = eigenvlues(modenumber-1);
}

double TimeEvolution::trap(const arma::vec &value){
  double integral = 0;
  double N = value.n_elem;
  for(int i = 1; i < value.n_elem-1; i++){
    integral += 2*value(i);
  }
  integral += value(0);
  integral += value(value.n_elem-1);

  return integral/(2*N);
}

void TimeEvolution::getMatrixWithTime(double t, arma::mat &result){
  result = matrix*equationTime(t);
}

////////////////////////////////////////////////////
double WaveEquation::equationTime(double t){
  return cos(sqrt(eigenvalue)*t);
}

double DiffusionEquation::equationTime(double t){
  return exp(-eigenvalue*t);
}












//////////////////////////////////////////////////////////////////
void TimeEvolutionDebug::checkIntegral(){
  arma::vec vector(6);
  vector.fill(1);
  cout << trap(vector) << endl;
}
