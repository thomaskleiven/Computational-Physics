#include "Schrodinger.hpp"
#include "Potentials.hpp"
#include "InitialConditions.hpp"
#include <stdlib.h>
#include <armadillo>

using namespace std;

int main(int argc, char* argv[]){
  FreeParticalPotential freePartical;
  InitialCondition init;

  if(argc != 2){cout << "Wrong arguments" << endl; return 1;}

  Schrodinger schrodinger(atoi(argv[1]));
  schrodinger.initDiagonals(freePartical);
  schrodinger.eigenvalueSolver();
  schrodinger.project(init);
  schrodinger.checkOrtogonality();
  std::vector<arma::cx_vec> modes;
  arma::cx_vec solution;
  solution.set_size( schrodinger.alpha_coeff.n_elem );
  int nt = 7000;
  double dt = 1.0/nt;
  for(int i = 0; i<nt; i++){
    double t = i*dt;
    for(int j = 0; j < schrodinger.alpha_coeff.n_elem; j++){
      arma::vec eigenvector = schrodinger.eigenvectors.col(j);
      schrodinger.setEigenvectorWithTime(eigenvector, t, j);
      modes.push_back(schrodinger.complex_eigenvector);
      solution += modes[j]*schrodinger.alpha_coeff(j);
    }
  }
  arma::vec real_solution = arma::conv_to<arma::vec>::from(arma::pow(arma::abs(solution),2));
  real_solution.save("solution.csv", arma::csv_ascii);
};
