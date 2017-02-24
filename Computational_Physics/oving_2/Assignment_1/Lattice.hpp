#ifndef LATTICE
#define LATTICE

#include<armadillo>
#include<vector>


struct Bond{
  int startPos;
  int neighbor;
};


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
  int getRootNode(int site);
  void activateBond(Bond &bond);
  int getBiggestCluster(int position, int neighbor);
  int getSmallesCluster(int position, int neighbor);
  std::vector<Bond> bonds;
  std::vector<int> sites;

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

class HoneycombLattice: public Lattice{
public:
  HoneycombLattice(int N): Lattice(N){};
protected:
  void findNeighbor(int position);
  bool checkIfEvenNumber(int position);
  bool checkIfEvenRow(int position);
  bool checkIfLastRow(int position);
  bool checkIfLastColumn(int position);
  bool checkIfFirstRow(int pos);
};

////////////////////////////////////////////////////////////////////////////////


class DebugLattice: public HoneycombLattice{
public:
  DebugLattice(int N):HoneycombLattice(N){};
  void printBonds();
};

#endif
