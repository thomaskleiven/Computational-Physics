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

  try{
    if(argc != 3){cout << "Wrong number of arguemts" << endl;return 1;}

    //Decide lattice fron user input
    string latticeType(argv[2]);
    LatticeType_t lt;
    if(latticeType == "square"){lt = LatticeType_t::SQUARE;}
    else if(latticeType == "triangular"){lt = LatticeType_t::TRIANGULAR;}
    else if(latticeType == "honey"){lt = LatticeType_t::HONEYCOMB;}
    else{
      cout << "Unknown lattice" << endl;
      return 1;
    }
    //Create an instance
    Activator act(atoi(argv[1]), lt);

    //Folder ID
    act.uid = argv[1];

    //Number of averages
    act.run_loops(1000);}
    
    catch (exception &ex){ cout << ex.what() << endl; return 1;}
  return 0;
};
