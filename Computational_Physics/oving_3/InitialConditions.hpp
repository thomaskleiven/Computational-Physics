#ifndef INITIALCONDITIONS
#define INITIALCONDITIONS
#include<armadillo>

class InitialCondition{
public:
  virtual double operator() (double x) const = 0;
};

class FirstEigenmode: public InitialCondition{
public:
  double operator()(double x) const override{
    return sqrt(2)*sin(4*atan(1)*x);
  }
};

class Dirac: public InitialCondition{
public:
  double operator()(double x) const override{
    return exp(-pow(x-x0,2)/sigma)/0.4;
  };
  void setSigma(double newSigma){sigma = newSigma;};
  void setCenter(double x){
    x0 = x;
  };
private:
  double sigma{0.01};
  double x0{0.2};
};

class FirstTwoEigenmode: public InitialCondition{
public:
  arma::mat eigvecs;
  bool status = eigvecs.load("eigenvectors/eigenvector_with_potential_1000.csv", arma::csv_ascii);
  arma::vec eigvec_n_one = eigvecs.col(0);
  arma::vec eigvec_n_two = eigvecs.col(2);
  double operator()(double x) const override{
    return (1.0)/sqrt(2)*( eigvec_n_one(x*eigvec_n_one.n_elem) + eigvec_n_two(x*eigvec_n_two.n_elem ));
  }
};

#endif
