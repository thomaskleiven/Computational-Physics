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
  struct Bond{int startPos; int neighbor;};
  void isCluster(Bond bond);
  std::vector<Bond*> bonds;

};

class SquareLattice: public Lattice{
public:
  SquareLattice(int N): Lattice(N){};
private:
  void findNeighbor(int position);
};

class TriangularLattice: public Lattice{
public:
  TriangularLattice(int N): Lattice(N){};
private:
  void findNeighbor(int position);
}


class DebugLattice: public SquareLattice{
public:
  DebugLattice(int N):SquareLattice(N){};
  void printBonds();
};

#endif
