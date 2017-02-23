#ifndef LATTICE
#define LATTICE

#include<armadillo>
#include<vector>

class Lattice{
public:
  Lattice(int N);

protected:
  virtual void findNeighbor(int position) = 0;
  void generateNeighbors();
  int N{0};
  arma::vec sites;
  bool checkIfLastRow(int position);
  bool checkIfLastColumn(int position);

};

class SquareLattice: public Lattice{
public:
  SquareLattice(int N): Lattice(N){};
  struct Bond{int startPos; int neighbor; int neighbor1;};
  std::vector<Bond*> bonds;
private:
  void findNeighbor(int position);
};

class TriangularLattice: public Lattice{
public:
  TriangularLattice(int N): Lattice(N){};
protected:
  struct Bond{int startPos; int neighbor; int neighbor1;};
  std::vector<Bond*> bonds;
  void findNeighbor(int position);
};


class DebugLattice: public TriangularLattice{
public:
  DebugLattice(int N):TriangularLattice(N){};
  void printBonds();
};

#endif
