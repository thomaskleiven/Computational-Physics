#include<iostream>
#include<vector>
#include<armadillo>
extern "C" void dpttrf_(int* n, double* diag, double* e, int* info);
extern "C" void dpttrs_(int* n, int* nrhs, double* diag, double* subdiag, double* rhs, int* nx, int* info);
#include <math.h>
#include<sstream>

using namespace std;



//Declare variables
static double endtime = 0.02;

static int nx = 401;                       //Number of length steps

static double dx = 1.0/(nx+1);                //Size of length steps
static double dt = 1.0/2E7;           //Size of time steps

static int nt = endtime/dt;

static arma::vec u(nx);                   //unknown u at new time level
static arma::vec u_l(nx);                 //u at previous time levelS2
static double F = dt/(dx*dx);

static arma::vec diagonal(nx);             //Diagonal
static arma::vec subDiagonal(nx-1);          //SubDiagonal

static arma::vec diagonalB(nx);           //Diag and Subdiag for matrix B in Crank Nicoloson
static arma::vec subDiagonalB(nx-1);

static int counter = 0;

static float xMin = -0.5;


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
      return 0.1;
    }
  }
};


double getX(double i){
  return xMin + (i+1)*dx;
}


template<class Diffusion>
void generateA(Scheme_t scheme, Diffusion &D){
  u_l.fill(0);
  u.fill(0);
  for(int i = 0; i<nx; i++){
    double fPlus = D(getX(i+0.5))*dt/(dx*dx);
    double fMinus = D(getX(i-0.5))*dt/(dx*dx);

    //cout << i*dx << " " <<  D(getX(i+0.5)) << " " << xMin << endl;
    //cout << D(getX(i+0.5)) << endl;

    switch (scheme) {

      case Scheme_t::FORWARD_EULER:
        break;

      case Scheme_t::BACKWARDS_EULER:
        diagonal(i) = 1 + fPlus + fMinus;        //Set diagonal and subdiagonal
        if(i < nx-1){
        subDiagonal(i) = -fMinus;}
        break;

      case Scheme_t::CRANK_NICOLSON:
        diagonal(i) = (1+0.5*fPlus+0.5*fMinus);
        diagonalB(i) = (1-0.5*fPlus-0.5*fMinus);
        if(i > 0){
          subDiagonal(i-1) = (-fMinus/2);
          subDiagonalB(i-1) = fMinus/2;
          }
        break;
        }

  }
}


void setInitialCondition(arma::vec &u){
  u(nx/2) = 1/dx;
}






void explicitEuler(bool reflective){

  ConstantDiffusion D;

  double F = D(0)*dt/(dx*dx);

  generateA(Scheme_t::FORWARD_EULER, D);

  setInitialCondition(u_l);

  arma::mat results(nx, 1);

  for (int i=0; i<nt;i++){
    for (int j=1; j<nx-1; j++){
        u(j) = u_l(j) + F*(u_l(j-1) - 2*u_l(j) + u_l(j+1));
    }
    if(reflective){
      u(0) = u(1);
      u(nx-1) = u(nx-2);
    }

    //Write to file every 5th iteration
    if(abs(0.07-i*dt) <1E-10){
      results.col(counter++) = u;
    }

    //Update u_l before next step
    u_l = u;
  }
  stringstream fname;
  fname << "explicitEuler_" << static_cast<int>(dt*1E6) << ".csv";
  results.save(fname.str().c_str(), arma::csv_ascii);
}

template<class Diffusion>
void setReflectiveBoundriesForImplicitEuler(arma::vec &diag, Diffusion &D){

  double F = D(0)*dt/(dx*dx);                               //Get F when X equals zero
  diag(0) =1.0 + F;

  F = D(getX(diag.n_elem-1))*dt/(dx*dx);                    //Get F when X euals nx-1
  diag(diag.n_elem-1) = 1.0 + F;
}

void setReflectiveForCrank(arma::vec &diag, arma::vec &diagB){

  double F = dt/(dx*dx);
  diag(0) = 1.0 + F/2;                         //Diagonal A

  diag(diag.n_elem-1) = 1.0 + F/2;         //Diagonal A

  diagB(0) = 1.0 - F/2;

  diagB(diagB.n_elem-1) = 1.0 - F/2;
}


void implicitEuler(bool reflective){
  int info;

  ConstantDiffusion D;                                                                                 //Has to be initialized in order to use LAPACK, returns SUCCESS or FAILURE

  generateA(Scheme_t::BACKWARDS_EULER, D);                                                               //Generate matrix

  setInitialCondition(u);

  arma::mat results(nx, 1);

  if(reflective){
    setReflectiveBoundriesForImplicitEuler(diagonal, D);
  }

  dpttrf_(&nx, diagonal.memptr(), subDiagonal.memptr(), &info);                               //Compute the symmetric system Ax = B using LAPACK

  //Solve the tridiag problem
  for(int i = 0; i<nt; i++){
    int nrhs = 1;
    dpttrs_(&nx, &nrhs, diagonal.memptr(), subDiagonal.memptr(), u.memptr(), &nx, &info);     // Solve the system Ax


    if(abs(0.07-i*dt) <1E-20){
      results.col(counter++) = u;
    }

  }
  stringstream fname;
  fname << "implicitEuler_" << static_cast<int>(dt*1E6) << ".csv";
  results.save(fname.str().c_str(), arma::csv_ascii);
}



  void crankNicolson(bool reflective){
    int info;
    arma::mat results(nx, 1);

    ConstantDiffusion D;

    StepDiffusion S;

    generateA(Scheme_t::CRANK_NICOLSON, S);

    setInitialCondition(u);

    //cout << diagonalB(0) << " " << diagonalB(diagonalB.n_elem -1) << endl;

    if(reflective){
      setReflectiveForCrank(diagonal, diagonalB);
    }

    //cout << diagonalB(0) << " " << diagonalB(diagonalB.n_elem -1) << endl;


    //cout << subDiagonal << endl;
    //Compute factorisation
    dpttrf_(&nx, diagonal.memptr(), subDiagonal.memptr(), &info);



    for(int i = 0; i < nt; ++i){
      int nrhs = 1;
      u_l = u;
      for(int j=0; j < nx-2; ++j){
        u(j+1) = subDiagonalB(j)*u_l(j) + diagonalB(j+1)*u_l(j+1) + subDiagonalB(j+1)*u_l(j+2);
      }
      u(0) = u_l(0)*diagonalB(0) + subDiagonalB(0)*u_l(1);
      u(u.n_elem-1) = u_l(u.n_elem-1)*diagonalB(diagonalB.n_elem-1) + u_l(u.n_elem-2)*subDiagonalB(subDiagonalB.n_elem-1);

      //Solve the system Ax = Bb
      dpttrs_(&nx, &nrhs, diagonal.memptr(), subDiagonal.memptr(), u.memptr(), &nx,  &info);

      if(abs(0.01-i*dt) <1E-20){
        results.col(counter++) = u;
      }
    }
    //stringstream fname;
    //fname << "crankNicolson_" << static_cast<int>(dt*1E6) << ".csv";
    results.save("crankNicolson_dx.csv", arma::csv_ascii);
}






int main(){
  //explicitEuler(false);
  //implicitEuler(false);
  crankNicolson(false);
  return 0;
}
;
