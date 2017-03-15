#include "Lattice.hpp"
#include <armadillo>
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <iterator>
#include <cstdint>
#include <fstream>
extern "C" double gsl_sf_lnfact(unsigned int n);
extern "C" double gsl_sf_fact(unsigned int n);


using namespace std;

Lattice::Lattice(int N):N(N), sites(N*N){
  fill(sites.begin(), sites.end(), -1);
  binomialCoeff.set_size(2*N*N);
  binomialCoeff.fill(0);
  pValue += getPvalue();
  pValueSquared += pow(getPvalue(),2);
  chi = getChi(0);
}


void Lattice::generateNeighbors(){
  for(int i = 0; i < sites.size(); i++){
    findNeighbor(i);
  }
}

void Lattice::activateBond(Bond &bond){
  pushBinomialCoeff();
  if(getRootNode(bond.neighbor) != getRootNode(bond.startPos)){
    int largest = (sites[getRootNode(bond.neighbor)] < sites[getRootNode(bond.startPos)]) ? getRootNode(bond.neighbor) : getRootNode(bond.startPos);
    int smallest = (sites[getRootNode(bond.neighbor)] > sites[getRootNode(bond.startPos)]) ? getRootNode(bond.neighbor) : getRootNode(bond.startPos);

    average_s -= pow(sites[largest], 2);                                        //Subtract merged cluster size squared
    average_s -= pow(sites[smallest], 2);

    if(smallest == largest){smallest = getRootNode(bond.neighbor);}             //If cluster size is equal
    sites[largest] += sites[smallest];                                          //Expand largest cluster with the size of the smallest
    sites[smallest] = getRootNode(largest);                                     //Change rootnode of the smallest cluster
    if(sites[largest] < sites[largestCluster]){largestCluster = largest;}       //Check if new cluster now is larger than the previous largest cluster

    average_s += pow(sites[largest], 2);                                        //Add new cluster size to average clustersize squared

    pValue += getPvalue();
    pValueSquared += pow(getPvalue(),2);

    chi = getChi(num_activatedBonds);
  }
  num_activatedBonds++;
}

void Lattice::activateSites(){
  shuffleBonds();
  shuffleBonds();
  shuffleBonds();
  lnFacBond = gsl_sf_lnfact(bonds.size());
  for(int i = 0; i<bonds.size(); i++){
    activateBond(bonds[i]);
    if(i == 9990*9990){
      getMainCluster();
      save(mainCluster);
    }
  }
}

double Lattice::getAverageClusterSize(){
  if(getPvalue() != 1){return (average_s - pow(sites.size()*getPvalue(),2))/(sites.size()*(1-getPvalue()));}
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
  return sites.size()*sqrt((pValueSquared/(i+1)) - pow(pValue/(i+1),2));

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
  ofstream output_file("./res.csv");
  ostream_iterator<int> output_iterator(output_file, "\n");
  copy(vector.begin(), vector.end(), output_iterator);
}

void Lattice::shuffleBonds(){
  random_shuffle (bonds.begin(), bonds.end());
}

void Lattice::pushBinomialCoeff(){
  double coeff =  lnFacBond - gsl_sf_lnfact(num_activatedBonds) - gsl_sf_lnfact(bonds.size()-num_activatedBonds);
}



void TriangularLattice::findNeighbor(int position){
  Bond bond;
  Bond bond1;
  Bond bond2;
  bond.startPos = position;
  bond1.startPos = position;
  bond2.startPos = position;
  bond.startPos = position;
  bond.neighbor = (position+1)%N + (position/N)*N;
  bond1.startPos = position;
  bond1.neighbor = (position+N)%(N*N);
  if(position == N*N-1){
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

void SquareLattice::findNeighbor(int position){
  Bond bond;
  Bond bond1;
  bond.startPos = position;
  bond.neighbor = (position+1)%N + (position/N)*N;
  bond1.startPos = position;
  bond1.neighbor = (position+N)%(N*N);

  bonds.push_back(bond);
  bonds.push_back(bond1);
}



void HoneycombLattice::findNeighbor(int position){
  Bond bond;
  Bond bond1;
  bond.startPos = position;
  bond1.startPos = position;

  if(position == N*N-1){
    bond.neighbor = 0;
    bond1.neighbor = position - 2*N +1;
    bonds.push_back(bond);
    bonds.push_back(bond1);
  }else if(checkIfFirstRow(position) && !checkIfEvenNumber(position)){
    bond.neighbor = position + N - 1;
    bond1.neighbor = position + N*N - N - 1;
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
      bond1.neighbor = position - N*N + N + 1;
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

bool HoneycombLattice::checkIfEvenNumber(int position){if(position%2 == 0){return true;} return false;}
bool HoneycombLattice::checkIfEvenRow(int position){if((position/N)%2 == 0){return true;} return false;}
bool HoneycombLattice::checkIfLastRow(int position){if((position/N) == N-1){return true;} return false;}
bool HoneycombLattice::checkIfLastColumn(int position){if((position+1)%N == 0){return true;} return false;}
bool HoneycombLattice::checkIfFirstRow(int position){if(position/N == 0){return true;} return false;}







/////////////////////////////////////////////////////////
void DebugLattice::printBonds(){
  generateNeighbors();

  //for (int i = 0; i<bonds.size(); i++){
    ////cout << "Startposition: " <<bonds[i].startPos << " Bond-to: " << bonds[i].neighbor << endl;
  //}
}

void DebugLattice::printSites(){
  activateSites();
  ////cout << sites.size() << endl;
  ////cout << getAverageClusterSize() << endl;
  for (int i = 0; i < sites.size(); i++){
    ////cout << "Site " << i << ": " << sites[i] << endl;
    ////cout << "Pvalue: " << pValues(i) << endl;
    ////cout << "PValueSquared: " << pValuesSquared(i) << endl;
    ////cout << "Chi: " << chi(i) << endl;
    ////cout << averageClusterSize[i] << endl;
  }
  /*for (int i = 0; i < mainCluster.size(); i++){
    ////cout << "MainSite " << i << ": " << mainCluster[i] << endl;
  }*/
}
