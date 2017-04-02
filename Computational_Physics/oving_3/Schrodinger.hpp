#ifndef SCHRODINGER
#define SCHRODINGER
#include <armadillo>

class Schrodinger{
public:
  Schrodinger(int nx);
  template<class V>
  void initDiagonals(const V &potential);
  void eigenvalueSolver();
protected:
  arma::vec diagonal;
  arma::vec sub_diagonal;
  template<class V>
  void buildDiag(const V &potential);
  void buildSubDiag();

  int nx{0};
  double dx;
};

#include "Schrodinger.tpp"
#endif
