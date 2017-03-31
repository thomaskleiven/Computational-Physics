#include "Activator.hpp"
#include "Lattice.hpp"
#include<armadillo>

extern "C" double gsl_sf_lnfact(unsigned int n);
extern "C" double gsl_sf_fact(unsigned int n);


using namespace std;

Activator::Activator(int N):N(N){
  n_sites = N*N;
  averaged_p_inf_values.set_size(2*n_sites);
  averaged_chi_values.set_size(2*n_sites);
  averaged_avg_clusterSize.set_size(2*n_sites);
  averaged_p_inf_sq_values.set_size(2*n_sites);

  averaged_p_inf_values.fill(0);
}


void Activator::run_loops(int n_loops){
  #pragma omp parallel for
  for (int i = 0; i < n_loops; i++){
    SquareLattice* sq = new SquareLattice(N);
    sq->generateNeighbors();
    sq->activateSites();
    for (int k=0; k<averaged_p_inf_values.n_elem; k++){
    averaged_p_inf_values[k] += sq->p_inf_values[k];
    averaged_chi_values[k] += sq->chi_values[k];
    averaged_p_inf_sq_values[k] += sq->p_inf_sq_values[k];
    averaged_avg_clusterSize[k] += sq->avg_clusterSize[k];}
    delete sq;
  }
  averaged_p_inf_values /= n_loops;
  averaged_chi_values /= n_loops;
  averaged_avg_clusterSize /= n_loops;
  calculateConvolution();
}

void Activator::calculateConvolution(){
  arma::vec p = arma::linspace(0.0, 1.0, 1E4);
  arma::vec convolution_p( p.n_elem );
  arma::vec convolution_chi( p.n_elem );
  arma::vec convolution_avg( p.n_elem );
  convolution_p.fill(0);
  convolution_avg.fill(0);
  convolution_chi.fill(0);
  unsigned int M = averaged_p_inf_values.n_elem;
  double lnFacBond = gsl_sf_lnfact(M);
  float percentage = 0;

  #pragma omp parallel for
  for(int i = 0; i<p.n_elem;i++){
    percentage += (1.0/(p.n_elem));
    if(omp_get_thread_num() == 0){cout << "\r Percentage: " << percentage;}
    for(int n=0;n<M; n++){
      double lnfac_n = gsl_sf_lnfact(n);
      double lnfac_M_n = gsl_sf_lnfact(M-n);
      //if(i == p.n_elem/10){binomialproba(n) = exp(lnFacBond - lnfac_n - lnfac_M_n + n*log(p(i))+(M - n)*log(1-p(i)));}
      convolution_p(i) += averaged_p_inf_values(n)*(exp(lnFacBond - lnfac_n - lnfac_M_n + n*log(p(i)) + (M - n)*log(1-p(i))));
      //convolution_avg(i) += averaged_avg_clusterSize(n)*exp(lnFacBond - lnfac_n - lnfac_M_n + n*log(p(i))+(M - n)*log(1-p(i)));
      //convolution_chi(i) += averaged_chi_values(n)*(exp(lnFacBond - lnfac_n - lnfac_M_n + n*log(p(i))+(M - n)*log(1-p(i))));
    }
  }
  stringstream fname;
  fname << "p" << ".csv";
  convolution_p.save(fname.str().c_str(), arma::csv_ascii);
  fname.str("");
  fname << "chi" << ".csv";
  convolution_chi.save(fname.str().c_str(), arma::csv_ascii);
  fname.str("");
  fname <<"avg" << ".csv";
  convolution_avg.save(fname.str().c_str(), arma::csv_ascii);
}
