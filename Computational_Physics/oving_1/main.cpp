#include<iostream>
#include<vector>
#include<armadillo>
extern "C" void dpttrf_(int* n, double* diag, double* e, int* info);
extern "C" void dpttrs_(int* n, int* nrhs, double* diag, double* subdiag, double* rhs, int* nx, int* info);

using namespace std;

void explicitEuler(int nt, int nx, bool reflective){

  double dx = 1.0/nx;
  double dt = 0.1*dx*dx*1;

  arma::vec u(nx);          //unknown u at new time level
  arma::vec u_l(nx);        //u at previous time levelS

  int counter = 0;
  arma::mat results(nx, 5);

  u_l.fill(0);
  u_l(nx/2) = 1;                 //Initial condition

  double F = dt/(dx*dx);
  for (int i=0; i<nt;i++){
    for (int j=1; j<nx-1; j++){
      if(!reflective){
        u(j) = u_l[j] + F*(u_l[j-1] - 2*u_l[j] + u_l[j+1]);}
    }

    //Boundry conditions
    u[0]  = 0;
    u[nx-1] = 0;

    if((i*5)%nt == 0){
      results.col(counter++) = u;
    }

    //Update u_l before next step
    u_l = u;
  }
  results.save("explicitEuler.csv", arma::csv_ascii);
}


void implicitEuler(int nt, int nx){
  double dx = 1.0/nx;
  double dt = 0.1*dx*dx*1;
  int info;


  arma::vec u(nx);          //unknown u at new time level

  arma::vec diagonal(nx);           //Diagonal
  arma::vec subDiagonal(nx);        //SubDiagonal

  u.fill(0);
  u(nx/2) = 1;                 //Initial condition


  int counter = 0;
  arma::mat results(nx, 5);


  double F = dt/(dx*dx);


  //Set diagonal and subdiagonal
  diagonal.fill(1 + 2*F);
  subDiagonal.fill(-F);

  dpttrf_(&nx, diagonal.memptr(), subDiagonal.memptr(), &info);


  //Compute b
  for(int i = 0; i<nt; i++){
    int nrhs = 1;
    dpttrs_(&nx, &nrhs, diagonal.memptr(), subDiagonal.memptr(), u.memptr(), &nx, &info);

    if((i*5)%nt == 0){
      results.col(counter++) = u;
    }

  }




  results.save("implicitEuler.csv", arma::csv_ascii);
  //Update u_l before next step
  //u_l = u;
}
  //results.save("implicitEuler.csv", arma::csv_ascii);



void crankNicolson(int nx, int nt){

  double dx = 1.0/nx;
  double dt = 0.1*dx*dx*1;

  double F = dt/(dx*dx);

  arma::vec Adiag(nx);               //Diagonal
  arma::vec AsubDiagonal(nx);        //SubDiagonal

  arma::vec Bdiag(nx);               //Diagonal
  arma::vec BsubDiagonal(nx);        //SubDiagonal

  arma::vec u(nx);                   //unknown u at time level

  u.fill(0);
  u(nx/2) = 1;                       //Initial condition

  Adiag.fill(1+F);
  AsubDiagonal.fill(-F/2);
  Bdiag.fill(1-F);
  BsubDiagonal.fill(F/2);



}





int main(){
  explicitEuler(100, 400000, false);
  //implicitEuler(100,40000);
  return 0;
}
