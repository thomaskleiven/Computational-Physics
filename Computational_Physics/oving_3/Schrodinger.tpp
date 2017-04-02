#ifndef SCHRODINGER_TPP
#define SCHRODINGER_TPP
template<class V>
void Schrodinger::initDiagonals(const V &potential){
  buildDiag(potential);
  buildSubDiag();
}


template<class V>
void Schrodinger::buildDiag(const V &potential){
  for (int i = 0; i < nx; i++){
      diagonal(i) = 2.0/(dx*dx) + potential(i*dx);
  }
}
#endif
