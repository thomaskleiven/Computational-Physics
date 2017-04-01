#include"Lattice.hpp"
#include"Activator.hpp"
#include<iostream>
#include<cstdlib>
#include<ctime>
#include <stdlib.h>
#include<stdexcept>

//#define DEBUG
using namespace std;


int main(int argc, char* argv[]){
  srand(time(0));
  #ifdef DEBUG
    //DebugLattice debug(100);
    //debug.makeLattice();
    //debug.printBonds();
    //debug.activateSites();
  #else
  try{
    if(argc != 3){cout << "Wrong number of arguemts" << endl;return 1;}
    string latticeType(argv[2]);
    LatticeType_t lt;
    if(latticeType == "square"){lt = LatticeType_t::SQUARE;}
    else if(latticeType == "triangular"){lt = LatticeType_t::TRIANGULAR;}
    else if(latticeType == "honey"){lt = LatticeType_t::HONEYCOMB;}
    else{
      cout << "Unknown lattice" << endl;
      return 1;
    }
    Activator act(atoi(argv[1]), lt);
    act.uid = argv[1];
    act.run_loops(1000);}
    catch (exception &ex){ cout << ex.what() << endl; return 1;}

    //SquareLattice solver(atoi(argv[1]));
    //TriangularLattice solver(atoi(argv[1]));
    //HoneycombLattice solver(atoi(argv[1]));
    //solver.uid = argv[1];
    //solver.folder = argv[2];
    //solver.generateNeighbors();
    //solver.run_loops(90);
  #endif
  return 0;
};
