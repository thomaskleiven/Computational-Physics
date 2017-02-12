#include "TimeEvolution.hpp"
#include<string>
#include<iostream>
//#define DEBUG

using namespace std;


int main(){

  TimeEvolution evolution;
  GaussianSpike spike;
  spike.setSigma(0.001);
  spike.setCenter(0.5,0.5);

  #ifdef DEBUG
    TimeEvolutionDebug debug;
    debug.checkIntegral();
    return 0;
  #endif
  string fname = "data/eigenvector_4.csv";
  string eigenvalueFile = "data/eigenvalues.csv";
  evolution.loadMatrix(fname);
  evolution.loadEigenvalue(eigenvalueFile, 4);
  cout << evolution.project(spike) << endl;

}
