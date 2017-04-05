#ifndef SCHRODINGER
#define SCHRODINGER
#include <armadillo>
#include <cmath>
#include <complex>
typedef std::complex<double> cdouble;

enum class InitialConditions_t{
  EIGENMODE, DRIAC
};

class Schrodinger{
public:
  Schrodinger(int nx);
  template<class V>
  void initDiagonals(const V &potential);
  void eigenvalueSolver();
  template<class Function>
  void project(const Function &condition);
  arma::cx_vec& getEigenvectorWithTime(arma::vec &eigenvector, double t, int eigenvalue);
  void save(arma::mat &eigenvectors, arma::vec &eigenvalues);
  arma::cx_vec getComplexEigenvector();
  arma::vec alpha_coeff;
  arma::mat eigenvectors;
  arma::cx_vec complex_eigenvector;
  double getMaxEigenvalue();
  double getMinEigenvalue();
  double trapezoidal(const arma::vec &eigenvector);
  double getNormalizationFactor(arma::vec &eigenvector, int size, int n);
  void normalizeEigenvectors();
  arma::mat normalized_eigvecs;

private:
  arma::vec diagonal;
  double normalization_factor{1};
  arma::vec sub_diagonal;
  arma::vec total_solution;
  double dx;
  template<class V>
  void buildDiag(const V &potential);
  void buildSubDiag();
  int nx{0};
  cdouble IMUNIT;

};

#include "Schrodinger.tpp"
#endif
