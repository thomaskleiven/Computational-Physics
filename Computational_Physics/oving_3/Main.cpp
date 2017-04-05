#include "Schrodinger.hpp"
#include "Potentials.hpp"
#include "InitialConditions.hpp"
#include <stdlib.h>
#include <armadillo>
#include <cmath>

using namespace std;

int main(int argc, char* argv[]){
  if(argc != 3){cout << "Wrong arguments" << endl; return 1;}
  FreeParticalPotential freePartical;
  Dirac dirac;
  FirstEigenmode mode;
  InitialCondition* impuls = 0;

  string initial(argv[2]);
  if(initial == "dirac"){
    impuls = &dirac;
  }else if(initial == "first_mode"){
    impuls = &mode;
  }else{
    cout << "Unknown initial condition" << endl; return 1;
  }

  Schrodinger schrodinger(atoi(argv[1]));
  schrodinger.initDiagonals(freePartical);
  schrodinger.eigenvalueSolver();
  schrodinger.normalizeEigenvectors();
  schrodinger.project(*impuls);
  arma::cx_vec solution;
  solution.set_size( schrodinger.alpha_coeff.n_elem );
  double nt = 100;
  double dt = 5000*(1.0/schrodinger.getMaxEigenvalue());
  arma::mat time_evolution(solution.n_elem, nt);
  time_evolution.fill(0);
  double sum_alpha_squared = arma::accu(arma::pow(schrodinger.alpha_coeff,2));

  for(int i = 0; i<nt; i++){
    double t = i*dt;
    solution.fill(0.0);
    for(int j = 0; j < schrodinger.alpha_coeff.n_elem; j++){
      arma::vec eigenvector = schrodinger.eigenvectors.col(j);
      solution += schrodinger.getEigenvectorWithTime(eigenvector, t, j)*schrodinger.alpha_coeff(j);
    }
    time_evolution.col(i) = arma::pow(arma::abs(solution), 2)/sum_alpha_squared;
  }

  time_evolution = time_evolution.t();
  time_evolution.save("time_evolution.csv", arma::csv_ascii);
  arma::vec psi = arma::pow(arma::abs(solution), 2) / sum_alpha_squared;
  psi.save("psi.csv", arma::csv_ascii);
};
