#include "Activator.hpp"
#include "Lattice.hpp"
#include<armadillo>

extern "C" double gsl_sf_lnfact(unsigned int n);
extern "C" double gsl_sf_fact(unsigned int n);


using namespace std;

Activator::Activator(int N):N(N){
  n_sites = N*N;
  p_inf_values.set_size(1.5*n_sites);
  chi_values.set_size(1.5*n_sites);
  avg_clusterSize.set_size(1.5*n_sites);
  p_inf_sq_values.set_size(1.5*n_sites);
  binomial.set_size(1.5*n_sites);
  num_of_bonds = 1.5*n_sites;
}


void Activator::run_loops(int n_loops){
  #pragma omp parallel for
  for (int k = 0; k < n_loops; k++){                                            //Number of times to average over
    HoneycombLattice* sq = new HoneycombLattice(N);
    sq->generateNeighbors();
    sq->shuffleBonds();
    for(int i= 0; i<sq->bonds.size(); i++){                                     //Activating one by one bond
      if(k==0){binomial(i) = sq->binomial_coeff(i);};
      avg_clusterSize(i) += sq->calcAverageClusterSize(sq->bonds[i]);
      p_inf_values(i) += sq->getPvalue();
      p_inf_sq_values(i) += pow(sq->getPvalue(), 2);
      chi_values(i) += sq->getChi(i, p_inf_values, p_inf_sq_values);
    }
    delete sq;
  }
  p_inf_values /= n_loops;
  chi_values /= n_loops;
  avg_clusterSize /= n_loops;
  calculateConvolution();                                                       //Calculate convolution
}


void Activator::calculateConvolution(){
  arma::vec p = arma::linspace(0.0, 1.0, 1E4);
  arma::vec convolution_p( p.n_elem );
  arma::vec convolution_chi( p.n_elem );
  arma::vec convolution_avg( p.n_elem );
  convolution_p.fill(0);
  convolution_avg.fill(0);
  convolution_chi.fill(0);
  float percentage = 0;
  #pragma omp parallel for
  for(int i = 0; i<p.n_elem;i++){
    percentage += (1.0/(p.n_elem));
    if(omp_get_thread_num() == 0){cout << "\r Percentage: " << percentage;}
    for(int n=0;n<num_of_bonds; n++){
      convolution_p(i) += p_inf_values(n)*(exp(binomial(n) + n*log(p(i)) + (num_of_bonds - n)*log(1-p(i))));
      convolution_avg(i) += avg_clusterSize(n)*exp(binomial(n) + n*log(p(i))+(num_of_bonds - n)*log(1-p(i)));
      convolution_chi(i) += chi_values(n)*(exp(binomial(n) + n*log(p(i))+(num_of_bonds - n)*log(1-p(i))));
    }
  }
  stringstream fname;
  fname << folder << "/" << "p" << uid << ".csv";
  convolution_p.save(fname.str().c_str(), arma::csv_ascii);
  fname.str("");
  fname << folder<< "/" << "chi" << uid << ".csv";
  convolution_chi.save(fname.str().c_str(), arma::csv_ascii);
  fname.str("");
  fname << folder << "/"<<"avg" << uid << ".csv";
  convolution_avg.save(fname.str().c_str(), arma::csv_ascii);
}
