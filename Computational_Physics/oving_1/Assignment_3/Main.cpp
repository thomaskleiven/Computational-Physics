#include"Hopf.hpp"

using namespace std;


int main(){
  ConservativeLaxWendroff solver(1000);
  solver.setInitialCondition();
  solver.transport();

  return 0;
};
