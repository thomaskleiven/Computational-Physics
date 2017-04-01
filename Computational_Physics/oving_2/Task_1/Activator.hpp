#ifndef ACTIVATOR
#define ACTIVATOR

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

  std::vector<double> binomialCoeff;
  arma::vec p_inf_values;
  arma::vec p_inf_sq_values;
  arma::vec avg_clusterSize;
  arma::vec chi_values;
  arma::vec binomial;
  int num_of_bonds{0};
  double lnFacBond;
  void pushBinomialCoeff();
  void calculateConvolution();
  void calculateChi();
  int n_sites;
  int N{0};
  void activateBonds(int i);
  int count{0};
  double average_p_inf_value{0};
  double average_p_inf_value_squared{0};
  void calculateChi(arma::vec &convolution_p, arma::vec &convolution_p_inf_squared);
  arma::vec p_val;
  LatticeType_t lattice{LatticeType_t::SQUARE};
  bool check{true};
  void checkOutput();
  Lattice *grid{NULL};

};





#endif
