#include"Lattice.hpp"
#include<iostream>
#define DEBUG


using namespace std;

int main(){
  //Lattice solver(3);
  #ifdef DEBUG
    DebugLattice debug(3);
    debug.printBonds();
  #endif
};
