#ifndef LATTICE
#define LATTICE

#include<armadillo>
#include<vector>

class Lattice{
public:
  Lattice(int N);

protected:
  int N{0};
  virtual void findNeighbor(int position) = 0;
  void shuffleBonds();
  void generateNeighbors();
  bool checkIfLastRow(int position);
  bool checkIfLastColumn(int position);
  struct Bond{int startPos; int neighbor; int neighbor1;};
  std::vector<Bond*> bonds;
  arma::vec sites;

};

class SquareLattice: public Lattice{
public:
  SquareLattice(int N): Lattice(N){};

protected:
  void findNeighbor(int position);
};

class TriangularLattice: public Lattice{
public:
  TriangularLattice(int N): Lattice(N){};

protected:
  void findNeighbor(int position);
};

////////////////////////////////////////////////////////////////////////////////


class DebugLattice: public TriangularLattice{
public:
  DebugLattice(int N):TriangularLattice(N){};
  void printBonds();
};

#endif
