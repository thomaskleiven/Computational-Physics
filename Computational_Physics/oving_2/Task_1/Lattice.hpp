
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
  void activateSites();
  std::string uid;
  std::string folder;
  void run_loops(int n_loops);
  arma::vec p_inf_values;
  arma::vec p_inf_sq_values;
  arma::vec chi_values;
  arma::vec avg_clusterSize;

protected:
  virtual void findNeighbor(int position) = 0;
  void shuffleBonds();
  void calculateAverageClusterSize(int i);
  void activateBond(Bond &bond);
  void save(std::vector<int> &vector);
  void saveCoeff(std::vector<double> &vector);
  void calcAverageClusterSize(Bond &bond);
  void getMainCluster();
  void buildNodeTable();
  void saveGrid(const char* fname);
  int getRootNode(int site);
  int num_activatedBonds{0};
  int largestCluster{0};
  int N{0};
  double average_s{0};
  double lastValue{0};
  int num{0};
  double getPvalue();
  double getAverageClusterSize();
  double getChi(int i);
  double pValue;
  double pValueSquared;
  double chi;
  double expected_s{0};
  int n_sites{0};
  int bondsActivated{0};
  std::vector<Coordinate> crd;
  std::vector<Bond> bonds;
  std::vector<int> sites;
  std::vector<int> mainCluster;
private:
  unsigned int count{0};
  std::vector<double> binomialCoeff;
  double lnFacBond;
  void pushBinomialCoeff();
  void calculateConvolution();
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
