#ifndef HELMHOLTZ2D_H
#define HELMHOLTZ2D_H
#include<armadillo>

class HelmHoltz2D{
public:
  HelmHoltz2D(unsigned int N_Discretization);

  ~HelmHoltz2D();
  void buildMatrix();

  void save();

  void solve(double number_of_eigenvalues);

protected:
  unsigned int N{0};                   //Number of discretication points in each direction

  arma::umat locations;
  arma::vec values;
  arma::sp_mat *matrix{NULL};
  arma::vec eigval;
  arma::mat eigvec;


  /** Sets the diagonal element */
  void diagonal();

  /**Sets the subDiagonal element */
  void subDiagonal();

  /**Sets the superDiagonal element*/
  void superDiagonal();

  /**Sets the supersuperDiagonal element*/
  void supsupDiag();

  /**Sets the subsubDiagonal element*/
  void subsubDiag();

  void initSparseMatrix();

  void initDiagonals();
};


class HelmHoltz2DDebug: public HelmHoltz2D{
public:
  HelmHoltz2DDebug(unsigned int N): HelmHoltz2D(N){};
  void printMatrix();

  void printLocations();

  void checkOrtogonality(double numberOfEigenvalues);

};



#endif
