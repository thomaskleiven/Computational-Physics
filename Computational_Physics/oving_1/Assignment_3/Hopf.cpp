#include "Hopf.hpp"
#include<armadillo>
#include<iostream>
#include<sstream>
#include<string>
#include<cmath>

using namespace std;

int Hopf::nt = 10000;
double Hopf::dx = 1.0;
double Hopf::dt = 0.00025;



Hopf::Hopf(int nx, const char* name):name(name){
  dx = 10.0/nx;
  solution.set_size(nx);
  solution.fill(0);
}

Hopf::Hopf(int nx):Hopf(nx, "HopfE"){};

void Hopf::save(arma::vec &sol, string const &name){
  sol.save(name, arma::csv_ascii);
}

void Hopf::setInitialCondition(){
  for (int i = 0; i< solution.n_elem; i++){
    double x = i*dx;
    //solution(i) = (1-(4*(x-0.5)*(x-0.5))) > 0 ? (1-(4*(x-0.5)*(x-0.5))) : 0;
    solution(i) = exp(-0.5*pow((x-5),2));
  }
}



void Upwind::applyScheme(){
  arma::vec u_last(solution);
  solution(solution.n_elem/2) = 0.5;
  solution(0) = -getDxDt()*(u_last(0)) - u_last(0);
  for (int i = 1; i<solution.n_elem;i++){
    solution(i) = -getDxDt()*(u_last(i) - u_last(i-1)) + u_last(i);
  }

}

void Downwind::applyScheme(){
  arma::vec u_last(solution);
  for(int i = 0; i<(solution.n_elem-2);i++){
    solution(i) = -getDxDt()*(u_last(i+1) - u_last(i)) + u_last(i);
  }
  solution(solution.n_elem-1) = -getDxDt()*u_last(solution.n_elem-1) + u_last(solution.n_elem-1);
}

void ExplicitCentered::applyScheme(){
  arma::vec u_last(solution);
  for(int i = 1; i<(solution.n_elem-2);i++){
    solution(i) = (-getDxDt()/2)*(u_last(i+1) - u_last(i-1)) + u_last(i);
  }
  solution(0) = -(getDxDt()/2)*(u_last(1)) + u_last(0);
  solution(solution.n_elem-1) = -(getDxDt()/2)*(u_last(solution.n_elem-2)) + u_last(solution.n_elem-1);
}

void ImplicitCentered::applyScheme(){
  return;
}

void LaxFriedrich::applyScheme(){
  arma::vec u_last(solution);
  for(int i = 1; i<(solution.n_elem-2); i++){
    solution(i) = -(getDxDt()*(u_last(i+1) - u_last(i-1)) + u_last(i+1) - u_last(i-1))/2;
}

  solution(0) = -(getDxDt()*u_last(1) + u_last(1))/2;
  solution(solution.n_elem-1) = -(getDxDt()*u_last(solution.n_elem-2) - u_last(solution.n_elem-2))/2;
}

void Hopf::transport(){
  exportSolution.set_size(solution.n_elem, n_save);
  int counter = 0;
  for(int i = 0; i<nt; i++){
    applyScheme();
    if((i*n_save)%nt == 0){
      exportSolution.col(counter) = solution;
      counter++;
      cout << counter << endl;
    }
  }
  string filename = name+".csv";
  exportSolution.save(filename, arma::csv_ascii);
}

void LaxWendroff::applyScheme(){
  arma::vec u_last(solution);
  for(int i = 1; i<(solution.n_elem-2); i++){
    solution(i) = ((dt/2.0)*(1.0/(dx*dx))*(u_last(i-1)-2*u_last(i) + u_last(i+1)) - (1.0/(2.0*dx))*(u_last(i+1) - u_last(i-1)))*dt + u_last(i);
  }
  solution(0) = ((dt/2.0)*(1.0/(dx*dx))*(-2*u_last(0) + u_last(1)) - (1.0/(2.0*dx))*(u_last(1)))*dt + u_last(0);
  solution(solution.n_elem-1) = ((dt/2.0)*(1.0/(dx*dx))*(u_last(solution.n_elem-2)-2*u_last(solution.n_elem-1)) - (1.0/(2.0*dx))*(- u_last(solution.n_elem-2)))*dt + u_last(solution.n_elem-1);
}


void ConservativeLaxWendroff::applyScheme(){
  arma::vec u_last(solution);
  double dtdx  = dt/(4*dx);
  double dtdx2 = pow(dt,2)/(8*pow(dx,2));
  double first;
  double second;
  double third;
  for (int i = 1; i<(solution.n_elem-2); i++){
    first = pow(u_last(i+1),2)-pow(u_last(i-1),2);
    second = (u_last(i+1) + u_last(i)) * (pow(u_last(i+1),2) - pow(u_last(i),2));
    third = (u_last(i) + u_last(i-1))*(pow(u_last(i),2) - pow(u_last(i-1),2));

    solution(i) = u_last(i) - dtdx*first + dtdx2*(second - third);
  }
}
