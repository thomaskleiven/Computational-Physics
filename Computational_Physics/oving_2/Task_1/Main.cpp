#include"Lattice.hpp"
#include<iostream>
#include<cstdlib>
#include<ctime>
#include <stdlib.h>

#define DEBUG
using namespace std;


int main(int argc, char* argv[]){
  srand(time(0));
  #ifdef DEBUG
    DebugLattice debug(50);
    debug.makeLattice();
    //debug.printBonds();
    debug.activateSites();
  #else
    if(argc != 3){cout << "Wrong number of arguemts" << endl;return 1;}
    //SquareLattice solver(atoi(argv[1]));
    //TriangularLattice solver(atoi(argv[1]));
    //HoneycombLattice solver(atoi(argv[1]));
    solver.uid = argv[1];
    solver.folder = argv[2];
    solver.generateNeighbors();
    solver.run_loops(90);
  #endif
  return 0;
};
