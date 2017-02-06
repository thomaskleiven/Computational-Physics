#include<iostream>
#include<vector>
#include<armadillo>
extern "C" void dpttrf_(int* n, double* diag, double* e, int* info);
extern "C" void dpttrs_(int* n, int* nrhs, double* diag, double* subdiag, double* rhs, int* nx, int* info);

using namespace std;

//Declare variables
static int nx = 40;                       //Number of length steps
static int nt = 5000;                     //Number of time steps

static double dx = 1.0/nx;                //Size of length steps
static double dt = 0.1*dx*dx*1;           //Size of time steps

static arma::vec u(nx);                   //unknown u at new time level
static arma::vec u_l(nx);                 //u at previous time levelS

//static double F = dt/(dx*dx);

static arma::vec diagonal(nx);             //Diagonal
static arma::vec subDiagonal(nx);          //SubDiagonal

static arma::vec diagonalB(nx);           //Diag and Subdiag for matrix B in Crank Nicoloson
static arma::vec subDiagonalB(nx);

static int counter = 0;

static int xMin = -0.5;


enum class Scheme_t{
  FORWARD_EULER, BACKWARDS_EULER, CRANK_NICOLSON
};

class ConstantDiffusion{
public:
  double operator()(double x){
    return 1.0;
  }
};

class StepDiffusion{
public:
  double operator()(double x){
    if(x > 0){return 1.0;}
    else{
      return 0.5;
    }
  }
};


double getX(int i){
  return xMin + i*dx;
}


template<class Diffusion>
void generateA(Scheme_t scheme, Diffusion &D){
  u_l.fill(0);
  u.fill(0);
  for(int i = 0; i<nx; i++){
    double fPlus = D(getX(i+0.5))*dt/(dx*dx);
    double fMinus = D(getX(i-0.5))*dt/(dx*dx);

    switch (scheme) {

      case Scheme_t::FORWARD_EULER:
        break;

      case Scheme_t::BACKWARDS_EULER:
        diagonal(i) = 1 + fPlus + fMinus;        //Set diagonal and subdiagonal
        subDiagonal(i) = -fMinus;
        break;

      case Scheme_t::CRANK_NICOLSON:
        diagonal(i) = (1+0.5*fPlus+0.5*fMinus);
        subDiagonal(i) = (-fMinus/2);
        diagonalB(i) = (1-0.5*fPlus-0.5*fMinus);
        subDiagonalB(i) = fMinus/2;
        break;
        }

  }
}


void setInitialCondition(arma::vec &u){
  u(nx/2) = 1;
}






void explicitEuler(bool reflective){

  ConstantDiffusion D;

  double F = D(0)*dt/(dx*dx);

  generateA(Scheme_t::FORWARD_EULER, D);

  arma::mat results(nx, nt);

  for (int i=0; i<nt;i++){
    for (int j=1; j<nx-1; j++){
        u(j) = u_l(j) + F*(u_l(j-1) - 2*u_l(j) + u_l(j+1));
    }
    if(reflective){
      u(0) = u(1);
      u(nx-1) = u(nx-2);
    }

    //Write to file every 5th iteration
    //if((i*5)%nt == 0){
      results.col(counter++) = u;
    //}

    //Update u_l before next step
    u_l = u;
  }
  results.save("explicitEuler.csv", arma::csv_ascii);
}

template<class Diffusion>
void setReflectiveBoundriesForImplicitEuler(arma::vec &diag, Diffusion &D){

  double F = D(0)*dt/(dx*dx);                               //Get F when X equals zero
  diag(0) =1.0 + F;

  F = D(getX(diag.n_elem-1))*dt/(dx*dx);                    //Get F when X euals nx-1
  diag(diag.n_elem-1) = 1.0 + F;
}

void setReflectiveForCrank(arma::vec &diag, arma::vec &diagB, arma::vec &subdiagB){

}


void implicitEuler(bool reflective){
  int info;

  ConstantDiffusion D;                                                                                 //Has to be initialized in order to use LAPACK, returns SUCCESS or FAILURE

  generateA(Scheme_t::BACKWARDS_EULER, D);                                                               //Generate matrix

  setInitialCondition(u);

  arma::mat results(nx, 5);

  if(reflective){
    setReflectiveBoundriesForImplicitEuler(diagonal, D);
  }

  dpttrf_(&nx, diagonal.memptr(), subDiagonal.memptr(), &info);                               //Compute the symmetric system Ax = B using LAPACK

  //Solve the tridiag problem
  for(int i = 0; i<nt; i++){
    int nrhs = 1;
    dpttrs_(&nx, &nrhs, diagonal.memptr(), subDiagonal.memptr(), u.memptr(), &nx, &info);     // Solve the system Ax


    if((i*5)%nt == 0){
      results.col(counter++) = u;
    }

  }

  results.save("implicitEuler.csv", arma::csv_ascii);
}



  void crankNicolson(bool reflective){
    int info;
    arma::mat results(nx, 5);

    ConstantDiffusion D;

    StepDiffusion S;

    generateA(Scheme_t::CRANK_NICOLSON, S);

    setInitialCondition(u);


    //Compute factorisation
    dpttrf_(&nx, diagonal.memptr(), subDiagonal.memptr(), &info);

    for(int i = 0; i < nt; ++i){
      int nrhs = 1;
      for(int i=0; i < nx-2; ++i){
        u(i+1) = subDiagonalB(i)*u(i) + diagonalB(i+1)*u(i+1) + subDiagonalB(i+1)*u(i+2);
      }

      //Solve the system Ax = Bb
      dpttrs_(&nx, &nrhs, diagonal.memptr(), subDiagonal.memptr(), u.memptr(), &nx,  &info);

      if((i*5)%nt == 0){
        results.col(counter++) = u;
      }
    }

    results.save("crankNicolson.csv", arma::csv_ascii);



  }






int main(){
  //explicitEuler(true);
  //implicitEuler(true);
  crankNicolson(false);
  return 0;
}
;
