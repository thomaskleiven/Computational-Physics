#include"Lattice.hpp"
#include<armadillo>
#include<cmath>
#include<iostream>
#include<vector>


using namespace std;

Lattice::Lattice(int N):N(N){
  sites.set_size(pow(N,2));
}


void Lattice::generateNeighbors(){
  for(int i = 0; i < sites.size(); i++){
    findNeighbor(i);
  }
}

void SquareLattice::findNeighbor(int position){
  Bond *bond = new Bond;
  Bond *bond1 = new Bond;
  bond->startPos = position;
  bond1->startPos = position;
  if((position > pow(N,2) - N - 1) && (position+1)%N == 0){
    bond->neighbor = position - N + 1;
    bond1->neighbor = position - N*N + N;
  }else if(position > pow(N,2) - N - 1){
    bond->neighbor = position - N*N + N;
    bond1->neighbor = position + 1;
  }else if((position+1)%N == 0){
    bond->neighbor = position - N + 1;
    bond1->neighbor = position + N;
  }else{
    bond->neighbor = position + 1;
    bond1->neighbor = position + sqrt(sites.n_elem);
  };

  bonds.push_back(bond);
  bonds.push_back(bond1);
}












/////////////////////////////////////////////////////////
void DebugLattice::printBonds(){
  generateNeighbors();
  for (int i = 0; i<bonds.size(); i++){
    cout << "Startposition: " <<bonds[i]->startPos << " Bond-to: " << bonds[i]->neighbor << endl;
  }
}
