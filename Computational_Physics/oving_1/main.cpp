#include<iostream>
#include<vector>
#include<armadillo>
extern "C" void dpttrf_(int* n, double* diag, double* e, int* info);
extern "C" void dpttrs_(int* n, int* nrhs, double* diag, double* subdiag, double* rhs, int* nx, int* info);



using namespace std;


enum class Scheme_t{
  FORWARD_EULER, BACKWARDS_EULER, CRANK_NICOLSON
};

class ConstantDiffusion{
public:
  double operator()(double x){
    return 1.0;
  }
};








void explicitEuler(int nt, int nx, bool reflective){

  double dx = 1.0/nx;
  double dt = 0.1*dx*dx*1;

  arma::vec u(nx);          //unknown u at new time level
  arma::vec u_l(nx);        //u at previous time levelS

  int counter = 0;
  arma::mat results(nx, nt);

  u_l.fill(0);
  u_l(nx/2) = 1;                 //Initial condition

  double F = dt/(dx*dx);

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


void implicitEuler(int nt, int nx, bool reflective){
  double dx = 1.0/nx;
  double dt = 0.1*dx*dx*1;
  int info;                           //Needed for LAPACK

  arma::vec u(nx);                    //unknown u at new time level

  arma::vec diagonal(nx);             //Diagonal
  arma::vec subDiagonal(nx);          //SubDiagonal

  u.fill(0);
  u(nx/2) = 1;                        //Initial condition

  int counter = 0;
  arma::mat results(nx, 5);

  double F = dt/(dx*dx);

  //Set diagonal and subdiagonal
  diagonal.fill(1 + 2*F);
  subDiagonal.fill(-F);

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



  void crankNicolson(int nt, int nx){
    double dx = 1.0/nx;
    double dt = 0.1*dx*dx*1;
    double F = dt/(dx*dx);

    int info;

    int counter = 0;

    arma::mat results(nx, nt);

    arma::vec u(nx);

    arma::vec diagonalA(nx);
    arma::vec subDiagonalA(nx-1);

    arma::vec diagonalB(nx);
    arma::vec subDiagonalB(nx-1);

    diagonalA.fill(1+F);
    subDiagonalA.fill(-F/2);
    diagonalB.fill(1-F);
    subDiagonalB.fill(F/2);
    u.fill(0);
    u(nx/2)=1;

    //Compute factorisation
    dpttrf_(&nx, subDiagonalA.memptr(), subDiagonalA.memptr(), &info);

    for(int i = 0; i < nt; ++i){
      int nrhs = 1;
      for(int i=0; i < nx-2; ++i){
        u(i+1) = subDiagonalB(i)*u(i) + diagonalB(i+1)*u(i+1) + subDiagonalB(i+1)*u(i+2);
      }
      //u(0) = diagonalB(0)*u(0)+subDiagonalB(0)*u(1);
      //u(nx-1) = diagonalB(0)*u(nx-1)+subDiagonalB(0)*u(nx-1);

      //Solve the system Ax = b
      dpttrs_(&nx, &nrhs, diagonalA.memptr(), subDiagonalA.memptr(), u.memptr(), &nx,  &info);

      //if((i*5)%nt == 0){
        results.col(counter++) = u;
      //}
    }

    results.save("crankNicolson.csv", arma::csv_ascii);



  }






int main(){
  explicitEuler(1000, 20, true);
  implicitEuler(300, 20, true);
  //crankNicolson(500,15);
  return 0;
}
;
