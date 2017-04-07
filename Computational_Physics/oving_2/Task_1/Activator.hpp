#ifndef ACTIVATOR
#define ACTIVATOR

#define ARMA_NO_DEBUG
#include<armadillo>

enum class LatticeType_t{
  SQUARE, TRIANGULAR, HONEYCOMB
};

class Lattice;

class Activator{
public:
  Activator(int N, LatticeType_t lattice);
  ~Activator();
  void run_loops(int n_loops);
  std::string uid;

private:

  bool check{true};
  unsigned int rand_identifier{0};
  int num_of_bonds{0};
  int n_sites;
  int N{0};
  int count{0};
  void checkOutput();
  void pushBinomialCoeff();
  void calculateConvolution();
  void calculateChi();
  void calculateChi(arma::vec &convolution_p, arma::vec &convolution_p_inf_squared);
  void activateBonds(int i);
  double average_p_inf_value{0};
  double average_p_inf_value_squared{0};
  double lnFacBond;
  arma::vec p_val;
  LatticeType_t lattice{LatticeType_t::SQUARE};
  Lattice *grid{NULL};
  std::vector<double> binomialCoeff;
  arma::vec p_inf_values;
  arma::vec p_inf_sq_values;
  arma::vec avg_clusterSize;
  arma::vec chi_values;
  arma::vec binomial;

};
#endif
