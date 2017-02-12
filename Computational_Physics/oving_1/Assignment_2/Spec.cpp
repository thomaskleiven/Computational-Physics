#include<iostream>
#include<vector>
#include<armadillo>
#include"HelmHoltz2D.hpp"
#include<pei/dialogBox.hpp>
#include<map>
#include<string>
//#define DEBUG

using namespace std;


int main(){
  map<string, double> params;
  params["size"] = 3.0;
  params["number_of_eigenvalues"] = 2.0;
  pei::DialogBox box(params);
  #ifdef DEBUG
    HelmHoltz2DDebug debugSolver(3);
    debugSolver.printMatrix();
    //debugSolver.printLocations();
    debugSolver.checkOrtogonality(3);
    debugSolver.save();
    return 0;
  #endif
  box.show();
  HelmHoltz2D solver(params.at("size"));
  solver.buildMatrix();
  solver.solve(params.at("number_of_eigenvalues"));
  solver.save();
  return 0;

}
;
