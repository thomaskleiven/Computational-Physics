#include "Schrodinger.hpp"
#include "Potentials.hpp"
#include <armadillo>
#include <stdlib.h>

int main(){
  BarrierPotential barrier;

  Schrodinger step_evolution;
  step_evolution.euler( barrier );

  return 0;
}
