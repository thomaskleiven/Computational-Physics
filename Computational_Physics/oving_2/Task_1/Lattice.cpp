#include "Lattice.hpp"
#include "Activator.hpp"
#include <armadillo>
#include <cmath>
#include <iostream>
#include <vector>
#include <iterator>
#include <stdexcept>
#include <sstream>
#include <fstream>
extern "C" double gsl_sf_lnfact(unsigned int n);
extern "C" double gsl_sf_fact(unsigned int n);

//#define DEBUG

using namespace std;

Lattice::Lattice(int N):N(N), sites(N*N){
  n_sites = sites.size();
  fill(sites.begin(), sites.end(), -1);
  p_inf_values.set_size(2*n_sites);
  p_inf_sq_values.set_size(2*n_sites);
  avg_clusterSize.set_size(2*n_sites);
  chi_values.set_size(2*n_sites);

  p_inf_values.fill(0);
  avg_clusterSize.fill(0);
  p_inf_sq_values.fill(0);
  chi_values.fill(0);
  average_s = n_sites;
}

void Lattice::generateNeighbors(){
  for(int i = 0; i < sites.size(); i++){
    findNeighbor(i);
  }
}

void Lattice::calcAverageClusterSize(Bond &bond){
  unsigned int index_rootnode_start = getRootNode(bond.startPos);
  unsigned int index_rootnode_end = getRootNode(bond.neighbor);
  if(index_rootnode_start != index_rootnode_end){

    unsigned int largest = (sites[index_rootnode_end] < sites[index_rootnode_start]) ? index_rootnode_end : index_rootnode_start;
    unsigned int smallest = (largest == index_rootnode_end) ? index_rootnode_start : index_rootnode_end;

    average_s -= pow(sites[largest], 2);                                             //Subtract merged cluster size squared
    average_s -= pow(sites[smallest], 2);

    sites[largest] += sites[smallest];                                               //Expand largest cluster with the size of the smallest
    sites[smallest] = largest;                                                       //Change rootnode of the smallest cluster

    if(sites[largest] < sites[largestCluster]){largestCluster = largest;}            //Check if new cluster now is larger than the previous largest cluster

    //Check if largest cluster found and if smallest cluster is merged to largest
    #ifdef DEBUG
      for(int i = 0; i<sites.size(); i++){
        if(sites[i] < sites[largestCluster]){
          throw logic_error("Not largest cluster");}}

      if(sites[largest] > sites[smallest]){throw logic_error("largest smaller than smallest");}

      if(largest != getRootNode(largest)){throw logic_error("Wrong rootnode");}

      if(largest == smallest){throw logic_error("Smallest and largest cannot be equal");}
    #endif

    average_s += pow(sites[largest], 2);                                                          //Add new cluster size to average clustersize squared

    if(getPvalue() < 1.0){expected_s = (double)(average_s - pow(n_sites*getPvalue(), 2))/(n_sites*(1-getPvalue()));}
    else{expected_s = 0;}
  }
  avg_clusterSize(num_activatedBonds) += expected_s;
}

void Lattice::run_loops(int n_loops){
  for (int i = 0; i < n_loops; i++){
    num_activatedBonds = 0;
    average_s = n_sites;
    fill(sites.begin(), sites.end(), -1);
    activateSites();
  }
  p_inf_values /= n_loops;
  avg_clusterSize /= n_loops;
  chi_values /= n_loops;
  //calculateConvolution();
}

void Lattice::activateSites(){
  shuffleBonds();
  lnFacBond = gsl_sf_lnfact(bonds.size());
  #pragma omp parallel for
  for(int i = 0; i<bonds.size(); i++){
    activateBond(bonds[i]);
  }
}

void Lattice::activateBond(Bond &bond){
  pushBinomialCoeff();                                                        //Calculate binomial coefficient after activating a bond
  calcAverageClusterSize(bond);                                               //Calc avg cluster size
  p_inf_values(num_activatedBonds) = getPvalue();                            //Calc p_inf value
  p_inf_sq_values(num_activatedBonds) = pow(getPvalue(), 2);                  //Calc p_inf squared value
  chi_values(num_activatedBonds) = getChi(num_activatedBonds);                //Calc chi value
  num_activatedBonds++;
}



int Lattice::getRootNode(int site){
  if(sites[site] < 0){
    return site;
  }else{
    getRootNode(sites[site]);
  }
}

double Lattice::getPvalue(){
  return (double)abs(sites[largestCluster])/sites.size();
}

double Lattice::getChi(int i){
  return sites.size()*sqrt(p_inf_values(i) - p_inf_sq_values(i));

}

void Lattice::getMainCluster(){
  for (int i = 0; i<sites.size(); i++){
    if(sites[i] == largestCluster){
      mainCluster.push_back(i);
    }else if(sites[i] >= 0 && getRootNode(i) == largestCluster){
      mainCluster.push_back(i);
    }
  }
  mainCluster.push_back(largestCluster);
}

void Lattice::save(vector<int> &vector){
  ofstream output_file("sites.csv");
  ostream_iterator<int> output_iterator(output_file, "\n");
  copy(vector.begin(), vector.end(), output_iterator);
}

void Lattice::shuffleBonds(){
  random_shuffle (bonds.begin(), bonds.end());
}

void Lattice::pushBinomialCoeff(){
  binomialCoeff.push_back(lnFacBond - gsl_sf_lnfact(num_activatedBonds) - gsl_sf_lnfact(bonds.size()-num_activatedBonds));
}

void Lattice::calculateConvolution(){
  arma::vec p = arma::linspace(0.0, 1.0, 1E4);
  arma::vec convolution_p( p.n_elem );
  arma::vec convolution_chi( p.n_elem );
  arma::vec convolution_avg( p.n_elem );
  convolution_p.fill(0);
  convolution_avg.fill(0);
  convolution_chi.fill(0);
  arma::vec binomialproba( p.n_elem );
  binomialproba.set_size( bonds.size() );
  unsigned int M = bonds.size();
  float percentage = 0;

  #pragma omp parallel for
  for(int i = 0; i<p.n_elem;i++){
    percentage += (1.0/(p.n_elem));
    if(omp_get_thread_num() == 0){cout << "\r Percentage: " << percentage;}
    for(int n=0;n<bonds.size(); n++){
      double lnfac_n = gsl_sf_lnfact(n);
      double lnfac_M_n = gsl_sf_lnfact(M-n);
      //if(i == p.n_elem/10){binomialproba(n) = exp(lnFacBond - lnfac_n - lnfac_M_n + n*log(p(i))+(M - n)*log(1-p(i)));}
      convolution_p(i) += p_inf_values(n)*(exp(lnFacBond - lnfac_n - lnfac_M_n + n*log(p(i)) + (M - n)*log(1-p(i))));
      convolution_avg(i) += avg_clusterSize(n)*exp(lnFacBond - lnfac_n - lnfac_M_n + n*log(p(i))+(M - n)*log(1-p(i)));
      convolution_chi(i) += chi_values(n)*(exp(lnFacBond - lnfac_n - lnfac_M_n + n*log(p(i))+(M - n)*log(1-p(i))));
    }
  }
  stringstream fname;
  fname <<folder << "p" << uid << ".csv";
  convolution_p.save(fname.str().c_str(), arma::csv_ascii);
  fname.str("");
  fname << folder<< "chi" << uid << ".csv";
  convolution_chi.save(fname.str().c_str(), arma::csv_ascii);
  fname.str("");
  fname <<folder<< "avg" << uid << ".csv";
  convolution_avg.save(fname.str().c_str(), arma::csv_ascii);
}





void TriangularLattice::findNeighbor(int position){
  setCoordinates();
  Bond bond;
  Bond bond1;
  Bond bond2;
  bond.startPos = position;
  bond1.startPos = position;
  bond2.startPos = position;
  bond.startPos = position;
  bond.neighbor = (position+1)%N + (position/N)*N;
  bond1.startPos = position;
  bond1.neighbor = (position+N)%(n_sites);
  if(position == n_sites-1){
    bond2.neighbor = 0;
  }else if((position+1) % N == 0){
    bond2.neighbor = position+1;
  }else if((position/N) == N-1){
    bond2.neighbor = position%N+1;
  }else{
   bond2.neighbor = position + N +1;
  }
  bonds.push_back(bond);
  bonds.push_back(bond1);
  bonds.push_back(bond2);
}

void TriangularLattice::setCoordinates(){
  crd.resize(N*N);
  for (int i = 0; i<N*N; i++){
        crd[i].x = i%N;
        crd[i].y = i/N;
        cout << crd[i].x << " " << crd[i].y << endl;
    }
}



void SquareLattice::findNeighbor(int position){
  setCoordinates();
  Bond bond;
  Bond bond1;
  bond.startPos = position;
  bond.neighbor = (position+1)%N + (position/N)*N;
  bond1.startPos = position;
  bond1.neighbor = (position+N)%(N*N);

  bonds.push_back(bond);
  bonds.push_back(bond1);
}

void SquareLattice::setCoordinates(){
  crd.resize(N*N);
  for (int i = 0; i<N*N; i++){
    crd[i].x = i%N;
    crd[i].y = i/N;
  }
}

void Lattice::saveGrid(const char* fname){
  ofstream out;
  out.open(fname);
  if(!out.good()){
    throw runtime_error("Could not open file in saveGrid");
  }
  for (int i= 0; i<bonds.size(); i++){
    out << crd[bonds[i].startPos].x << "," << crd[bonds[i].startPos].y << "," << crd[bonds[i].neighbor].x << "," << crd[bonds[i].neighbor].y <<"\n";


    //Square
    //out << crd[bonds[i].startPos].x << "," << crd[bonds[i].startPos].y << "," << crd[bonds[i].neighbor].x << "," << crd[bonds[i].neighbor].y << "\n";
  }
  out.close();
  cout << "Grid saved to file " << fname << endl;
}



void HoneycombLattice::findNeighbor(int position){
  setCoordinates();
  Bond bond;
  Bond bond1;
  bond.startPos = position;
  bond1.startPos = position;

  if(position == n_sites-1){
    bond.neighbor = 0;
    bond1.neighbor = position - 2*N +1;
    bonds.push_back(bond);
    bonds.push_back(bond1);
  }else if(checkIfFirstRow(position) && !checkIfEvenNumber(position)){
    bond.neighbor = position + N - 1;
    bond1.neighbor = position + n_sites - N - 1;
    bonds.push_back(bond);
    bonds.push_back(bond1);
  }else if(checkIfLastColumn(position) && !checkIfEvenNumber(position) && !checkIfEvenRow(position)){
    bond.neighbor = position - 2*N +1;
    bond1.neighbor = position + 1;
    bonds.push_back(bond);
    bonds.push_back(bond1);
  }else if(checkIfEvenNumber(position)){
    bond.neighbor = position + 1;
    bonds.push_back(bond);
  }else if(!checkIfEvenNumber(position) && !checkIfEvenRow(position)){
    bond.neighbor = position - N + 1;
    bond1.neighbor = position + N + 1;
    if(checkIfLastRow(position)){
      bond1.neighbor = position - n_sites + N + 1;
    }
    bonds.push_back(bond);
    bonds.push_back(bond1);
  }else if(!checkIfEvenNumber(position) && checkIfEvenRow(position)){
    bond.neighbor = position - N - 1;
    bond1.neighbor = position + N - 1;
    bonds.push_back(bond);
    bonds.push_back(bond1);
  }
}

void HoneycombLattice::setCoordinates(){
  crd.resize(N*N);
  for (int i = 0; i<N*N; i++){
    crd[i].x = i%N;
    crd[i].y = i/N;
  }
}

bool HoneycombLattice::checkIfEvenNumber(int position){if(position%2 == 0){return true;} return false;}
bool HoneycombLattice::checkIfEvenRow(int position){if((position/N)%2 == 0){return true;} return false;}
bool HoneycombLattice::checkIfLastRow(int position){if((position/N) == N-1){return true;} return false;}
bool HoneycombLattice::checkIfLastColumn(int position){if((position+1)%N == 0){return true;} return false;}
bool HoneycombLattice::checkIfFirstRow(int position){if(position/N == 0){return true;} return false;}







/////////////////////////////////////////////////////////
void DebugLattice::makeLattice(){
  generateNeighbors();
  saveGrid("HoneyLattice.csv");
}

void DebugLattice::activate(){
    activateSites();
}


void DebugLattice::printBonds(){
  for (int i=0; i<bonds.size(); i++){
    cout << bonds[i].startPos << " " << bonds[i].neighbor << endl;
  }
}
