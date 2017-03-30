#include"Lattice.hpp"
#include<iostream>
#include<cstdlib>
#include<ctime>
#include <stdlib.h>

//#define DEBUG
using namespace std;


int main(int argc, char* argv[]){
  srand(time(0));
  #ifdef DEBUG
    DebugLattice debug(40);
    debug.makeLattice();
  #else
    if(argc != 2){cout << "Wrong number of arguemts" << endl;return 1;}
    SquareLattice solver(atoi(argv[1]));
    solver.uid = argv[1];
    solver.generateNeighbors();
    solver.run_loops(90);
  #endif
  return 0;
};
