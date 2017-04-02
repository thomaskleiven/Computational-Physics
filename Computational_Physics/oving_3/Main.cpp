#include "Schrodinger.hpp"
#include "Potentials.hpp"

using namespace std;

int main(){
  FreeParticalPotential freePartical;
  Schrodinger schrodinger(100);
  schrodinger.initDiagonals(freePartical);
  schrodinger.eigenvalueSolver();
}
