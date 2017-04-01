#include "Lattice.hpp"
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

using namespace std;

double HoneycombLattice::bonds_per_site = 1.5;


Lattice::Lattice(int N):N(N), sites(N*N){
  n_sites = sites.size();
  fill(sites.begin(), sites.end(), -1);
  average_s = n_sites;
}

void Lattice::generateNeighbors(){
  for(int i = 0; i < sites.size(); i++){
    findNeighbor(i);
  }
}


double Lattice::calcAverageClusterSize(Bond &bond){
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

    average_s += pow(sites[largest], 2);                                                          //Add new cluster size to average clustersize squared

    if(getPvalue() < 1.0){expected_s = (double)(average_s - pow(n_sites*getPvalue(), 2))/(n_sites*(1-getPvalue()));}
    else{expected_s = 0;}
  }
  return expected_s;
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

double Lattice::getChi(int i, arma::vec &a, arma::vec &b){
  return sites.size()*sqrt(a(i) - b(i));

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

double Lattice::binomial_coeff(int num_activatedBonds){
  return gsl_sf_lnfact(bonds.size()) - gsl_sf_lnfact(num_activatedBonds) - gsl_sf_lnfact(bonds.size()-num_activatedBonds);
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
    }
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

void SquareLattice::setCoordinates(){
  crd.resize(N*N);
  for (int i = 0; i<N*N; i++){
    cout << "here" << endl;
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
}



void HoneycombLattice::findNeighbor(int position){
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

void DebugLattice::printBonds(){
  for (int i=0; i<bonds.size(); i++){
    cout << bonds[i].startPos << " " << bonds[i].neighbor << endl;
  }
}
