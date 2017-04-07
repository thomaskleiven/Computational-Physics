#include "Schrodinger.hpp"
#include "Potentials.hpp"
#include <armadillo>
#include <stdlib.h>

int main(){
  BarrierPotential barrier;

  Schrodinger step_evolution;
  //step_evolution.euler( barrier );
  step_evolution.setDiagForCrank( barrier );
  step_evolution.CrankNicolsonScheme();

  return 0;
}
