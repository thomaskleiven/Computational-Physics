template<class Function>
double TimeEvolution::project(Function &func){
  double dx = 1.0/(matrix.n_rows+1);
  double dy = 1.0/(matrix.n_cols+1);
  arma::vec initialConditionVector(matrix.n_rows);
  arma::vec integralVector(matrix.n_cols);

  for(int i = 0; i < matrix.n_cols; i++){
    double x = (i+1)*dx;
    for(int j = 0; j < matrix.n_rows; j++){
      double y = (j+1)*dy;
        initialConditionVector(j) = matrix(j,i) * func(x,y);
      }
      integralVector(i) = trap(initialConditionVector);
    }

    return trap(integralVector);

}
