#ifndef SCHRODINGER
#define SCHRODINGER
#include <armadillo>
#include <cmath>
#include <complex>
typedef std::complex<double> cdouble;

class Schrodinger{
public:
  Schrodinger(int nx);
  template<class V>
  void initDiagonals(const V &potential);
  void eigenvalueSolver();
  template<class Function>
  void project(const Function &condition);
  void setEigenvectorWithTime(arma::vec &eigenvector, double t, int eigenvalue);
  void checkOrtogonality();
  arma::cx_vec getComplexEigenvector();
  arma::vec alpha_coeff;
  arma::mat eigenvectors;
  arma::cx_vec complex_eigenvector;

private:
  arma::vec diagonal;
  arma::vec sub_diagonal;
  arma::vec total_solution;
  double dx;
  double trapezoidal(const arma::vec eigenvector);
  template<class V>
  void buildDiag(const V &potential);
  void buildSubDiag();
  int nx{0};
  cdouble IMUNIT;

};

#include "Schrodinger.tpp"
#endif
