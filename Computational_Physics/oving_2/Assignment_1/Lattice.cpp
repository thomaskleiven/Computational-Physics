#include"Lattice.hpp"
#include<armadillo>
#include<cmath>
#include<iostream>
#include<vector>
#include<algorithm>
#include<random>
#include<iterator>
#include <cstdint>


using namespace std;

Lattice::Lattice(int N):N(N), sites(N*N){
  fill(sites.begin(), sites.end(), -1);
}


void Lattice::generateNeighbors(){
  for(int i = 0; i < sites.size(); i++){
    findNeighbor(i);
  }
  //shuffleBonds();
}

void Lattice::activateBond(Bond &bond){
  if(getRootNode(bond.neighbor) != getRootNode(bond.startPos)){
    int largest = (sites[getRootNode(bond.neighbor)] < sites[getRootNode(bond.startPos)]) ? getRootNode(bond.neighbor) : getRootNode(bond.startPos) ;
    int smallest = (sites[getRootNode(bond.neighbor)] > sites[getRootNode(bond.startPos)]) ? getRootNode(bond.neighbor) : getRootNode(bond.startPos);
    if(smallest == largest){
      smallest = getRootNode(bond.neighbor);
    }
    sites[largest] += sites[smallest];
    sites[smallest] = getRootNode(largest);
  }
}
int Lattice::getRootNode(int site){
  if(sites[site] < 0){
    return site;
  }else{
    getRootNode(sites[site]);
  }
}

int Lattice::getBiggestCluster(int startPosition, int neighbor){
  return (startPosition > neighbor) ? startPosition : neighbor;
}

int Lattice::getSmallesCluster(int startPosition, int neighbor){
  return (startPosition < neighbor) ? startPosition : neighbor;
}


void Lattice::shuffleBonds(){
  std::random_shuffle ( bonds.begin(), bonds.end() );
}


bool Lattice::checkIfLastRow(int position){
  if((position+1) % N == 0){
    return true;
  }
}

bool Lattice::checkIfLastColumn(int position){
  if(position > pow(N,2) - N - 1){
    return true;
  }
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
  for (int i = 0; i<bonds.size(); i++){
    cout << "Startposition: " <<bonds[i].startPos << " Bond-to: " << bonds[i].neighbor << endl;
  }
}

void DebugLattice::printSites(){
  for (int i = 0; i < sites.size(); i++){
    cout << "Site " << i << ": " <<sites[i] << endl;
  }
}
