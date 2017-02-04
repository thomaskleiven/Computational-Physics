#include<iostream>
#include<vector>
#include<armadillo>
extern "C" void dpttrf_(int* n, double* diag, double* e, int* info);
extern "C" void dpttrs_(int* n, int* nrhs, double* diag, double* subdiag, double* rhs, int* nx, int* info);

using namespace std;

//Declare variables
static int nx = 40;                       //Number of length steps
static int nt = 10000;                     //Number of time steps

static double dx = 1.0/nx;                //Size of length steps
static double dt = 0.1*dx*dx*1;           //Size of time steps

static arma::vec u(nx);                   //unknown u at new time level
static arma::vec u_l(nx);                 //u at previous time levelS

static double F = dt/(dx*dx);

static arma::vec diagonal(nx);             //Diagonal
static arma::vec subDiagonal(nx);          //SubDiagonal

static arma::vec diagonalB(nx);           //Diag and Subdiag for matrix B in Crank Nicoloson
static arma::vec subDiagonalB(nx);

static int counter = 0;


enum class Scheme_t{
  FORWARD_EULER, BACKWARDS_EULER, CRANK_NICOLSON
};

class ConstantDiffusion{
public:
  double operator()(double x){
    return 1.0;
  }
};

void generateA(Scheme_t scheme){
  switch (scheme) {
    case Scheme_t::FORWARD_EULER:
      u_l.fill(0);
      u_l(nx/2) = 1;                 //Initial conditions
      return;

    case Scheme_t::BACKWARDS_EULER:
      u.fill(0);
      u(nx/2) = 1;                   //Initial condition

      diagonal.fill(1 + 2*F);        //Set diagonal and subdiagonal
      subDiagonal.fill(-F);
      return;

    case Scheme_t::CRANK_NICOLSON:
      diagonal.fill(1+F);
      subDiagonal.fill(-F/2);
      diagonalB.fill(1-F);
      subDiagonalB.fill(F/2);
      u.fill(0);
      u(nx/2)=1;
      return;
      }
}







void explicitEuler(bool reflective){

  generateA(Scheme_t::FORWARD_EULER);

  arma::mat results(nx, 5);

  for (int i=0; i<nt;i++){
    for (int j=1; j<nx-1; j++){
        u(j) = u_l(j) + F*(u_l(j-1) - 2*u_l(j) + u_l(j+1));
    }
    if(reflective){
      u(0) = u(1);
      u(nx-1) = u(nx-2);
    }

    //Write to file every 5th iteration
    if((i*5)%nt == 0){
      results.col(counter++) = u;
    }

    //Update u_l before next step
    u_l = u;
  }
  results.save("explicitEuler.csv", arma::csv_ascii);
}


void implicitEuler(bool reflective){
  int info;                           //Needed for LAPACK

  generateA(Scheme_t::BACKWARDS_EULER);

  arma::mat results(nx, 5);

  dpttrf_(&nx, diagonal.memptr(), subDiagonal.memptr(), &info);

  //Solve the tridiag problem
  for(int i = 0; i<nt; i++){
    int nrhs = 1;
    dpttrs_(&nx, &nrhs, diagonal.memptr(), subDiagonal.memptr(), u.memptr(), &nx, &info);
    if(reflective){
      u(0) = u(1);
      u(nx-1) = u(nx-2);
    }


    if((i*5)%nt == 0){
      results.col(counter++) = u;
    }

  }

  results.save("implicitEuler.csv", arma::csv_ascii);
}



  void crankNicolson(bool reflective){
    int info;
    arma::mat results(nx, nt);

    generateA(Scheme_t::CRANK_NICOLSON);

    //Compute factorisation
    dpttrf_(&nx, diagonal.memptr(), subDiagonal.memptr(), &info);

    for(int i = 0; i < nt; ++i){
      int nrhs = 1;
      for(int i=0; i < nx-2; ++i){
        u(i+1) = subDiagonalB(i)*u(i) + diagonalB(i+1)*u(i+1) + subDiagonalB(i+1)*u(i+2);
      }
      if(reflective){
        u(0) = u(1);
        u(nx-1) = u(nx-2);
      }

      //Solve the system Ax = Bb
      dpttrs_(&nx, &nrhs, diagonal.memptr(), subDiagonal.memptr(), u.memptr(), &nx,  &info);

      //if((i*5)%nt == 0){
        results.col(counter++) = u;
      //}
    }

    results.save("crankNicolson.csv", arma::csv_ascii);



  }






int main(){
  //explicitEuler(true);
  //implicitEuler(false);
  crankNicolson(true);
  return 0;
}
;
