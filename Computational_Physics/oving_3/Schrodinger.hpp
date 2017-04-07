#ifndef SCHRODINGER
#define SCHRODINGER
#include <armadillo>
#include <cmath>
#include <complex>
#include <sstream>
typedef std::complex<double> cdouble;

enum class InitialConditions_t{
  EIGENMODE, DRIAC
};

class Schrodinger{
public:
  Schrodinger(int nx);
  Schrodinger(){};
  template<class V>
  void initDiagonals(const V &potential);
  void eigenvalueSolver();
  template<class Function>
  void project(const Function &condition);
  void normalizeEigenvectors();
  template<class V>
  void euler(const V &potential);
  template<class V>
  void setDiagForCrank(const V& potential);
  void CrankNicolsonScheme();
  template<class init>
  void setInitialCondition(const init& condition);
  double getMinEigenvalue();
  double getMaxEigenvalue();
  double dx;
  arma::vec alpha_coeff;
  arma::vec diagonal;
  arma::mat eigenvectors;
  arma::cx_vec& getEigenvectorWithTime(arma::vec &eigenvector, double t, int eigenvalue);

private:
  template<class V>
  void buildDiag(const V &potential);
  void save(arma::mat &eigenvectors, arma::vec &eigenvalues);
  void buildSubDiag();
  void setSubDiagForCrank();
  double normalization_factor{1};
  double interpolation(arma::vec &eigenvector);
  double trapezoidal(const arma::vec &eigenvector);
  double *next_step_cast;
  arma::cx_vec getComplexEigenvector();
  arma::cx_vec u_last;
  arma::vec sub_diagonal;
  arma::vec total_solution;
  arma::cx_vec complex_eigenvector;
  arma::mat normalized_eigvecs;
  arma::cx_vec crank_diag_A;
  arma::cx_vec crank_diag_B;
  arma::cx_vec crank_sub_diag_A;
  arma::cx_vec crank_sub_diag_B;
  arma::cx_vec crank_super_diagonal_A;
  cdouble IMUNIT;

};

#include "Schrodinger.tpp"
#endif
