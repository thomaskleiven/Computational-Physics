#ifndef HOPF
#define HOPF
#include<armadillo>
#include<string>


class Hopf{
public:
  explicit Hopf(int nx);
  void setDiagonals();
  void save(arma::vec &sol, std::string const &name);
  void setInitialCondition();
  void transport();
  virtual void applyScheme() = 0;

protected:
  Hopf(int nx, const char* newname);
  static double dx;
  static double dt;
  static int nt;

  double getDxDt(){
    return dt/dx;
  }

  unsigned int n_save{5};

  arma::vec solution;

  std::string name;
};

class ImplicitCentered: public Hopf{
public:
  explicit ImplicitCentered(int nx):Hopf(nx){};
  void applyScheme();
private:
  arma::vec diagonal();

  arma::vec superDiagonal();

  arma::vec subDiagonal();

  arma::vec top();

  arma::vec bottom();
};

class Upwind: public Hopf{
public:
  explicit Upwind(int nx):Hopf(nx){};
  void applyScheme();

};

class Downwind: public Hopf{
public:
  explicit Downwind(int nx):Hopf(nx){};
  void applyScheme();
};

class ExplicitCentered: public Hopf{
public:
  explicit ExplicitCentered(int nx):Hopf(nx){};
  void applyScheme();
};

class LaxFriedrich: public Hopf{
public:
  explicit LaxFriedrich(int nx):Hopf(nx, "LaxFriedrich"){};
  void applyScheme() override;
};

class LaxWendroff: public Hopf{
public:
  explicit LaxWendroff(int nx):Hopf(nx, "laxWendroff"){};
  void applyScheme() override;
};

class ConservativeLaxWendroff: public Hopf{
public:
  explicit ConservativeLaxWendroff(int nx):Hopf(nx, "ConservativeLaxWendroff"){};
  void applyScheme() override;
};






#endif
