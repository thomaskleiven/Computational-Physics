#ifndef TIMEEVOLUTION
#define TIMEEVOLUTION
#include<armadillo>
#include<string>
#include<cmath>
#include<sstream>


class TimeEvolution{
public:
  TimeEvolution(){};
  void loadMatrix(const std::string &filename);
  void loadEigenvalue(const std::string &filename, unsigned int eigenmode);
  template<class Function>
  double project(Function &initialCondition);
  void getMatrixWithTime(double time, arma::mat &matrix);
  virtual double equationTime(double t) = 0;
protected:
  arma::mat matrix;
  double eigenvalue;
  static double trap(const arma::vec &value);

};

class WaveEquation: public TimeEvolution{
public:
  virtual double equationTime(double t) override;
};

class DiffusionEquation: public TimeEvolution{
public:
  virtual double equationTime(double t) override;
};

class GaussianSpike{
public:
  double operator()(double x, double y){
    return exp(-pow(x-x0,2)/sigma - pow(y-y0, 2)/sigma)/sigma;
  };
  void setSigma(double newSigma){sigma = newSigma;};
  void setCenter(double x, double y){
    x0 = x;
    y0 = y;
  };
private:
  double sigma{0.001};
  double x0{0.5};
  double y0{0.5};
};


class TimeEvolutionDebug: public TimeEvolution{
public:
  void checkIntegral();
  virtual double equationTime(double t) override{};
};

#include "TimeEvolution.tpp"

#endif
