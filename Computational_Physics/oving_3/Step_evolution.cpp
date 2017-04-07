#include "Schrodinger.hpp"
#include "Potentials.hpp"
#include "InitialConditions.hpp"
#include <armadillo>
#include <stdlib.h>

int main(int argc, char* argv[]){

  if(argc != 2){cout << "Wrong arguments" << endl; return 1;}

  //Potentials
  BarrierPotential barrier;

  //Initial conidtions
  FirstEigenmode mode;
  FirstTwoEigenmode modes;
  Dirac dirac;

  //Create an instance
  Schrodinger step_evolution;


  InitialCondition *impuls = 0;
  //Activate given intial condition
  string initial_condition(argv[1]);
  if(initial_condition == "1mode"){
    impuls = &mode;
  }else if(initial_condition == "2modes"){
    impuls = &modes;
  }else if(initial_condition == "dirac"){
    impuls = &dirac;
  }else{
    cout << "Unknown initial condition" << endl; return 1;
  }

  //Initialize system
  step_evolution.setInitialCondition( *impuls );
  step_evolution.setDiagForCrank( barrier );
  //Solve system.
  step_evolution.CrankNicolsonScheme();

  //step_evolution.euler( barrier );

  return 0;
}
