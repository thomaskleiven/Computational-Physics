#ifndef INITIALCONDITIONS
#define INITIALCONDITIONS

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
    return exp(-pow(x-x0,2)/sigma)/sigma;
  };
  void setSigma(double newSigma){sigma = newSigma;};
  void setCenter(double x){
    x0 = x;
  };
private:
  double sigma{0.01};
  double x0{0.5};
};

#endif
