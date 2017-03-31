
#ifndef LATTICE
#define LATTICE

#include<armadillo>
#include<vector>


struct Bond{
  int startPos;
  int neighbor;
};

struct Coordinate{
  int x;
  int y;
};


class Lattice{
public:
  Lattice(int N);
  void generateNeighbors();
  std::vector<Bond> bonds;
  double lnFacBond;
  void shuffleBonds();
  double getPvalue();
  double calcAverageClusterSize(Bond &bond);
  double getAverageClusterSize();
  double getChi(int i, arma::vec &a, arma::vec &b);
  double binomial_coeff(int num_activatedBonds);

protected:
  virtual void findNeighbor(int position) = 0;
  void save(std::vector<int> &vector);
  void getMainCluster();
  void buildNodeTable();
  void saveGrid(const char* fname);
  int getRootNode(int site);
  int num_activatedBonds{0};
  int largestCluster{0};
  int N{0};
  double average_s{0};
  double expected_s{0};
  int n_sites{0};
  std::vector<Coordinate> crd;
  std::vector<int> sites;
  std::vector<int> mainCluster;
};

class SquareLattice: public Lattice{
public:
  SquareLattice(int N): Lattice(N){};
protected:
  void findNeighbor(int position);
  void setCoordinates();
};

class TriangularLattice: public Lattice{
public:
  TriangularLattice(int N): Lattice(N){};
protected:
  void findNeighbor(int position);
  void setCoordinates();
};

class HoneycombLattice: public Lattice{
public:
  HoneycombLattice(int N): Lattice(N){};
protected:
  void setCoordinates();
  void findNeighbor(int position);
  bool checkIfEvenNumber(int position);
  bool checkIfEvenRow(int position);
  bool checkIfLastRow(int position);
  bool checkIfLastColumn(int position);
  bool checkIfFirstRow(int pos);
};


////////////////////////////////////////////////////////////////////////////////


class DebugLattice: public TriangularLattice{
public:
  DebugLattice(int N):TriangularLattice(N){};
  void makeLattice();
  void printBonds();
  void activate();
  //void printSites();
  //void checkIfLargestClusterFound();
};

#endif
