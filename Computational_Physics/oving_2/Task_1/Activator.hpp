#ifndef ACTIVATOR
#define ACTIVATOR

#include<armadillo>

class Activator{
public:
  Activator(int N);
  std::vector<double> binomialCoeff;
  arma::vec averaged_p_inf_values;
  arma::vec averaged_p_inf_sq_values;
  arma::vec averaged_avg_clusterSize;
  arma::vec averaged_chi_values;
  double lnFacBond;
  void pushBinomialCoeff();
  void calculateConvolution();
  void run_loops(int n_loops);
  int n_sites;
  int N{0};
};





#endif
