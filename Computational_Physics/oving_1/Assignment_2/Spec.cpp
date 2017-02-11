#include<iostream>
#include<vector>
#include<armadillo>
#include"HelmHoltz2D.hpp"
#include<pei/dialogBox.hpp>
#include<map>
#include<string>
#define DEBUG

using namespace std;


int main(){




  map<string, double> params;
  params["size"] = 3.0;
  pei::DialogBox box(params);
  #ifdef DEBUG
    HelmHoltz2DDebug debugSolver(4);
    debugSolver.printMatrix();
    //debugSolver.printLocations();
    return 0;
  #endif
  box.show();
  HelmHoltz2D solver(params.at("size"));
  solver.buildMatrix();



  return 0;

}
;
