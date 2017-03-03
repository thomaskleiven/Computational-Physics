#include "TimeEvolution.hpp"
#include<string>
#include<iostream>
#include<vector>
#include<sstream>
#include<visa/visa.hpp>
//#define DEBUG

using namespace std;


int main(){
  visa::WindowHandler plots;
  plots.addPlot("Eigenmode");
  plots.get("Eigenmode").setCmap( visa::Colormaps::Colormap_t::NIPY_SPECTRAL );
  plots.get("Eigenmode").setColorLim( 0, 100 );

  bool updateColorLimits = false;


  vector<TimeEvolution*> modes;
  vector<double> projectionsCoeff;
  GaussianSpike spike;
  spike.setSigma(0.1);
  spike.setCenter(0.5,0.5);

  for(int i = 1; i<30; i++){
    DiffusionEquation *evolution = new DiffusionEquation;
    stringstream fname;
    fname << "data/eigenvector_" << static_cast<int>(i) << ".csv";
    evolution->loadMatrix(fname.str());
    evolution->loadEigenvalue("data/eigenvalues.csv", i);
    modes.push_back(evolution);
    projectionsCoeff.push_back(evolution->project(spike));
  }

  arma::mat totalSolution;
  arma::mat modeSolution;
  unsigned int nt = 7000;
  double dt = 0.00005;
  for(int i = 0; i<nt; i++){
    double t = i*dt;
    for(int j = 0; j < modes.size(); j++){
      modes[j]->getMatrixWithTime(t, modeSolution);
      if(j==0){
        totalSolution = modeSolution*projectionsCoeff[j];
      }else{
        totalSolution += modeSolution*projectionsCoeff[j];}
    }
    if(updateColorLimits){
      double min = arma::min(arma::min(totalSolution));
      double max = arma::max(arma::max(totalSolution));
      plots.get("Eigenmode").setColorLim(min,max);
      updateColorLimits = false;
    }
    plots.get("Eigenmode").fillVertexArray(totalSolution);
    plots.show();
    plots.get("Eigenmode").clear();
  }

  for(int i = 0; i < modes.size(); i++){
    delete modes[i];
  }


  #ifdef DEBUG
    TimeEvolutionDebug debug;
    debug.checkIntegral();
    return 0;
  #endif

}
