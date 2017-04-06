#ifndef POTENTIALS
#define POTENTIALS

class FreeParticalPotential{
public:
  double operator()(double x) const{
    return 0;
  }
};

class BarrierPotential{
public:
  double operator()(double x) const{
    if(((1.0/3) < x) && (x < (2.0/3))){return 1E3;}else{return 0;}
  }
};


#endif
