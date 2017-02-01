#include<iostream>
#include<vector>
#include<cstdlib>
#include<iostream>
#include<fstream>
#include<iterator>
#include<vector>

using namespace std;

vector<int> list;

double fRand(double fMin, double fMax);
int randomWalk();
void writeToFile();
int count = 0;


int randomWalk(){
  double position = 0.0;
  int sign = 1;
  int maxIter = 1E6;
  for (int i =0; i < maxIter; i++){
    sign = position > 0.0 ? 1:-1;
    position += fRand(-1,1);
    int newsign = position > 0.0 ? 1:-1;
    //if(position > 1000){count++;}
    if ( newsign != sign && i>0 ) {return i;}
  }
  count++;
  return maxIter;
}



double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

void writeToFile(){
  ofstream output_file("./data.txt");
  ostream_iterator<int> output_iterator(output_file, "\n");
  copy(list.begin(), list.end(), output_iterator);
}



int main() {
  srand((unsigned)time(0));
  for (int i = 0; i<1E6;i++){
    list.push_back(randomWalk());
  }



  writeToFile();

  //vector<int> v { 34,23 };

  //cout << v[0] << endl;


}
