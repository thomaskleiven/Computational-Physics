#include"Lattice.hpp"
#include<iostream>
#include<cstdlib>
#include<ctime>

#define DEBUG


using namespace std;

int main(){
  srand(time(0));
  //Lattice solver(3);
  #ifdef DEBUG
    DebugLattice debug(10);
    debug.printBonds();
    debug.printSites();
  #endif
};
