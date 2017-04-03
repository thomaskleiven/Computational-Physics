#ifndef INITIALCONDITIONS
#define INITIALCONDITIONS

class InitialCondition{
public:
  double operator()(double x) const{
    return sin(4*atan(1)*x);
  }
};

#endif
