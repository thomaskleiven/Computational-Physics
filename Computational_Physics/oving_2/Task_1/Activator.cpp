#include "Activator.hpp"
#include "Lattice.hpp"
#include <armadillo>
#include <cstdlib>
#include <cassert>

extern "C" double gsl_sf_lnfact(unsigned int n);
extern "C" double gsl_sf_fact(unsigned int n);


using namespace std;

Activator::Activator(int N, LatticeType_t lattice):N(N), lattice(lattice){
  n_sites = N*N;
  rand_identifier = rand()%10000000;
  double n_bonds = 0;
  switch(lattice){
    case LatticeType_t::SQUARE:
        n_bonds = SquareLattice::bonds_per_site*N*N;
        break;
    case LatticeType_t::TRIANGULAR:
        n_bonds = TriangularLattice::bonds_per_site*N*N;
        break;
    case LatticeType_t::HONEYCOMB:
        n_bonds = HoneycombLattice::bonds_per_site*N*N;
        break;
  }
  p_inf_values.set_size(n_bonds);
  chi_values.set_size(n_bonds);
  chi_values.fill(0);
  avg_clusterSize.set_size(n_bonds);
  p_inf_sq_values.set_size(n_bonds);
  binomial.set_size(n_bonds);

  p_inf_sq_values.fill(0.0);
  p_inf_values.fill(0.0);
  avg_clusterSize.fill(0.0);

  num_of_bonds = n_bonds;
}

Activator::~Activator(){
  delete grid;
}


void Activator::run_loops(int n_loops){
  #pragma omp parallel for
  for (int k = 0; k < n_loops; k++){                                            //Number of times to average over
    Lattice *sq = NULL;
    unsigned int id = omp_get_thread_num();
    if(id==0){delete grid;}

    switch(lattice){
      case LatticeType_t::SQUARE:
        sq = new SquareLattice(N);
        break;
      case LatticeType_t::TRIANGULAR:
        sq = new TriangularLattice(N);
        break;
      case LatticeType_t::HONEYCOMB:
        sq = new HoneycombLattice(N);
        break;
      default: cout << "Unknown" << endl;
    }

    if(id == 0){grid = sq;}


    //if((id == 0) && check){checkOutput();}
    sq->generateNeighbors();
    sq->shuffleBonds();
    for(int i= 0; i<sq->bonds.size(); i++){                                     //Activating one by one bond
      if(k==0){binomial(i) = sq->binomial_coeff(i);};
      avg_clusterSize(i) += sq->calcAverageClusterSize(sq->bonds[i]);
      p_inf_values(i) += sq->getPvalue();
    }
    if(id != 0){delete sq; sq = NULL;}
    count++;
  }

  p_inf_values /= n_loops;
  p_inf_sq_values /= n_loops;
  avg_clusterSize /= n_loops;
  calculateConvolution();                                                     //Calculate convolution
}

void Activator::checkOutput(){
  assert(grid!=NULL);
  check = false;
  string fname = grid->folder;
  fname += "test.csv";
  ofstream out;
  out.open(fname.c_str());
  if(!out.good()){throw logic_error("Could not open output");}
  out.close();
}

void Activator::calculateChi(arma::vec &convolution_p, arma::vec &convolution_p_inf_squared){
  arma::vec chi ( convolution_p.n_elem );

  for(int i = 0; i<chi.n_elem; i++){
    if(convolution_p_inf_squared(i)-pow(convolution_p(i),2) < 0.0){
      chi(i) = 0;
    }else{
      double sub =  convolution_p_inf_squared(i)-pow(convolution_p(i),2);
      sub = sqrt(sub);
      chi(i) = sub*n_sites;}
  }
  stringstream filename;
  filename << grid->folder<< "/" << "chi" << uid << "_" << rand_identifier << ".csv";
  chi.save(filename.str().c_str(), arma::csv_ascii);
  filename.str("");
}


void Activator::calculateConvolution(){
  arma::vec p = arma::linspace(0.0, 1.0, 3E4);
  arma::vec convolution_p( p.n_elem );
  arma::vec convolution_p_inf_squared( p.n_elem );
  arma::vec convolution_avg( p.n_elem );
  arma::vec p_sq = arma::pow(p_inf_values, 2);
  convolution_p.fill(0);
  convolution_avg.fill(0);
  convolution_p_inf_squared.fill(0);
  float percentage = 0;
  #pragma omp parallel for
  for(int i = 0; i<p.n_elem;i++){
    percentage += (1.0/(p.n_elem));
    if(omp_get_thread_num() == 0){cout << "\r Percentage: " << percentage;}
    for(int n=0;n<num_of_bonds; n++){
      convolution_p(i) += p_inf_values(n)*(exp(binomial(n) + n*log(p(i)) + (num_of_bonds - n)*log(1-p(i))));
      convolution_avg(i) += avg_clusterSize(n)*exp(binomial(n) + n*log(p(i))+(num_of_bonds - n)*log(1-p(i)));
      convolution_p_inf_squared(i) += p_sq(n)*(exp(binomial(n) + n*log(p(i))+(num_of_bonds - n)*log(1-p(i))));
    }
  }

  calculateChi(convolution_p, convolution_p_inf_squared);

  stringstream fname;
  fname << grid->folder << "/" << "p" << uid << "_" << rand_identifier <<".csv";
  convolution_p_inf_squared.save(fname.str().c_str(), arma::csv_ascii);
  fname.str("");
  fname << grid->folder << "/" << "p" << uid << "_" << rand_identifier <<".csv";
  convolution_p.save(fname.str().c_str(), arma::csv_ascii);
  fname.str("");
  fname << grid->folder << "/"<<"avg" << uid << "_" << rand_identifier <<".csv";
  convolution_avg.save(fname.str().c_str(), arma::csv_ascii);
}
