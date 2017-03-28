#include "Lattice.hpp"
#include <armadillo>
#include <cmath>
#include <iostream>
#include <vector>
#include <iterator>
extern "C" double gsl_sf_lnfact(unsigned int n);
extern "C" double gsl_sf_fact(unsigned int n);


using namespace std;

Lattice::Lattice(int N):N(N), sites(N*N){
  fill(sites.begin(), sites.end(), -1);
  p_inf_values.set_size(2*N*N);
  p_inf_sq_values.set_size(2*N*N);
  avg_clusterSize.set_size(2*N*N);
  chi_values.set_size(2*N*N);
  average_s = N*N;
}


void Lattice::generateNeighbors(){
  for(int i = 0; i < sites.size(); i++){
    findNeighbor(i);
  }
}


void Lattice::activateBond(Bond &bond){
  pushBinomialCoeff();                                                        //Calculate binomial coefficient after activating a bond
  calcAverageClusterSize(bond);                                               //Calc avg cluster size
  p_inf_values(num_activatedBonds-1) = getPvalue();                           //Calc p_inf value
  p_inf_sq_values(num_activatedBonds-1) = pow(getPvalue(), 2);                //Calc p_inf squared value
  chi_values(num_activatedBonds-1) = getChi(num_activatedBonds-1);            //Calc chi value
  num_activatedBonds++;
}

void Lattice::calcAverageClusterSize(Bond &bond){
  unsigned int index_rootnode_start = getRootNode(bond.startPos);
  unsigned int index_rootnode_end = getRootNode(bond.neighbor);                 //Root node of end node
  if(index_rootnode_start != index_rootnode_end){

    unsigned int largest = (sites[index_rootnode_end] < sites[index_rootnode_start]) ? index_rootnode_end : index_rootnode_start;
    unsigned int smallest = (largest == index_rootnode_end) ? index_rootnode_start : index_rootnode_end;

    average_s -= pow(sites[largest], 2);                                                          //Subtract merged cluster size squared
    average_s -= pow(sites[smallest], 2);

    //cout << "Start position: " << bond.startPos << endl;
    //cout << "End position: " << bond.neighbor << endl;
    //cout << "Largest cluster: " << largest << "        SIZE: " << sites[largest] <<endl;
    //cout << "Smallest cluster: " << smallest << "        SIZE: " << sites[smallest] <<endl;
    //cout << "Size of new merged cluster: " << sites[largest] << endl;
    //cin.get();


    sites[largest] += sites[smallest];                                                            //Expand largest cluster with the size of the smallest
    sites[smallest] = largest;                                                                    //Change rootnode of the smallest cluster

    if(sites[largest] < sites[largestCluster]){largestCluster = largest;}                         //Check if new cluster now is larger than the previous largest cluster

    average_s += pow(sites[largest], 2);                                                          //Add new cluster size to average clustersize squared

    if(getPvalue() < 1.0){expected_s = (double)(average_s - pow(N*N*getPvalue(), 2))/(N*N*(1-getPvalue()));}
    else{expected_s = 0;}
  }
  avg_clusterSize(num_activatedBonds-1) = expected_s;
}

void Lattice::activateSites(){
  shuffleBonds();
  lnFacBond = gsl_sf_lnfact(bonds.size());
  for(int i = 0; i<bonds.size(); i++){
    activateBond(bonds[i]);
  }
  calculateConvolution();
}

int Lattice::getRootNode(int site){
  if(sites[site] < 0){
    return site;
  }else{
    getRootNode(sites[site]);
  }
}

double Lattice::getPvalue(){
  return (double)abs(sites[largestCluster])/sites.size();
}

double Lattice::getChi(int i){
  return sites.size()*sqrt(p_inf_values(i) - p_inf_sq_values(i));

}

void Lattice::getMainCluster(){
  for (int i = 0; i<sites.size(); i++){
    if(sites[i] == largestCluster){
      mainCluster.push_back(i);
    }else if(sites[i] >= 0 && getRootNode(i) == largestCluster){
      mainCluster.push_back(i);
    }
  }
  mainCluster.push_back(largestCluster);
}

/*void Lattice::save(vector<int> &vector){
  ofstream output_file("./res.csv");
  ostream_iterator<int> output_iterator(output_file, "\n");
  copy(vector.begin(), vector.end(), output_iterator);
}*/

void Lattice::shuffleBonds(){
  random_shuffle (bonds.begin(), bonds.end());
}

void Lattice::pushBinomialCoeff(){
  binomialCoeff.push_back(lnFacBond - gsl_sf_lnfact(num_activatedBonds) - gsl_sf_lnfact(bonds.size()-num_activatedBonds));
}

void Lattice::calculateConvolution(){
  arma::vec p = arma::linspace(0.0, 1.0, 1E4);
  arma::vec convolution_p( p.n_elem );
  arma::vec convolution_chi( p.n_elem );
  arma::vec convolution_avg( p.n_elem );
  unsigned int M = bonds.size();
  float percentage = 0;
  double lnfac_n;
  double lnfac_M_n;

  #pragma omp parallel for private(lnfac_n, lnfac_M_n)
  for(int i = 0; i<p.n_elem;i++){
    percentage += (1.0/(p.n_elem));
    if(omp_get_thread_num() == 0){
    cout << "\r Percentage: " << percentage;}
    convolution_p(i) = 0.0;
    for(int n=0;n<bonds.size(); n++){
      lnfac_n = gsl_sf_lnfact(n);
      lnfac_M_n = gsl_sf_lnfact(M-n);
      convolution_p(i) += p_inf_values(n)*(exp(lnFacBond - lnfac_n - lnfac_M_n + n*log(p(i))+(M - n)*log(1-p(i))));
      convolution_avg(i) += avg_clusterSize(n)*exp(lnFacBond - lnfac_n - lnfac_M_n + n*log(p(i))+(M - n)*log(1-p(i)));
      //convolution_chi(i) += chi_values[n]*(exp(lnFacBond - lnfac_n - lnfac_M_n + n*log(p(i))+(M - n)*log(1-p(i))));
    }
  }
  convolution_p.save("p.csv", arma::csv_ascii);
  //convolution_chi.save("chi.csv", arma::csv_ascii);
  convolution_avg.save("avg.csv", arma::csv_ascii);
  //avg_clusterSize.save("avg.csv", arma::csv_ascii);
}





void TriangularLattice::findNeighbor(int position){
  Bond bond;
  Bond bond1;
  Bond bond2;
  bond.startPos = position;
  bond1.startPos = position;
  bond2.startPos = position;
  bond.startPos = position;
  bond.neighbor = (position+1)%N + (position/N)*N;
  bond1.startPos = position;
  bond1.neighbor = (position+N)%(N*N);
  if(position == N*N-1){
    bond2.neighbor = 0;
  }else if((position+1) % N == 0){
    bond2.neighbor = position+1;
  }else if((position/N) == N-1){
    bond2.neighbor = position%N+1;
  }else{
   bond2.neighbor = position + N +1;
  }
  bonds.push_back(bond);
  bonds.push_back(bond1);
  bonds.push_back(bond2);
}

void SquareLattice::findNeighbor(int position){
  Bond bond;
  Bond bond1;
  bond.startPos = position;
  bond.neighbor = (position+1)%N + (position/N)*N;
  bond1.startPos = position;
  bond1.neighbor = (position+N)%(N*N);

  bonds.push_back(bond);
  bonds.push_back(bond1);
}



void HoneycombLattice::findNeighbor(int position){
  Bond bond;
  Bond bond1;
  bond.startPos = position;
  bond1.startPos = position;

  if(position == N*N-1){
    bond.neighbor = 0;
    bond1.neighbor = position - 2*N +1;
    bonds.push_back(bond);
    bonds.push_back(bond1);
  }else if(checkIfFirstRow(position) && !checkIfEvenNumber(position)){
    bond.neighbor = position + N - 1;
    bond1.neighbor = position + N*N - N - 1;
    bonds.push_back(bond);
    bonds.push_back(bond1);
  }else if(checkIfLastColumn(position) && !checkIfEvenNumber(position) && !checkIfEvenRow(position)){
    bond.neighbor = position - 2*N +1;
    bond1.neighbor = position + 1;
    bonds.push_back(bond);
    bonds.push_back(bond1);
  }else if(checkIfEvenNumber(position)){
    bond.neighbor = position + 1;
    bonds.push_back(bond);
  }else if(!checkIfEvenNumber(position) && !checkIfEvenRow(position)){
    bond.neighbor = position - N + 1;
    bond1.neighbor = position + N + 1;
    if(checkIfLastRow(position)){
      bond1.neighbor = position - N*N + N + 1;
    }
    bonds.push_back(bond);
    bonds.push_back(bond1);
  }else if(!checkIfEvenNumber(position) && checkIfEvenRow(position)){
    bond.neighbor = position - N - 1;
    bond1.neighbor = position + N - 1;
    bonds.push_back(bond);
    bonds.push_back(bond1);
  }
}

bool HoneycombLattice::checkIfEvenNumber(int position){if(position%2 == 0){return true;} return false;}
bool HoneycombLattice::checkIfEvenRow(int position){if((position/N)%2 == 0){return true;} return false;}
bool HoneycombLattice::checkIfLastRow(int position){if((position/N) == N-1){return true;} return false;}
bool HoneycombLattice::checkIfLastColumn(int position){if((position+1)%N == 0){return true;} return false;}
bool HoneycombLattice::checkIfFirstRow(int position){if(position/N == 0){return true;} return false;}







/////////////////////////////////////////////////////////
void DebugLattice::printBonds(){
  generateNeighbors();
}

void DebugLattice::printSites(){
  activateSites();
}
