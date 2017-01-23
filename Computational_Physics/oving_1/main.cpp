#include<iostream>
#include<vector>
#include<cstdlib>
#include<iostream>
#include<fstream>
#include<iterator>
#include<vector>

using namespace std;

vector<double> list;

double fRand(double fMin, double fMax);
void randomInteger();
void writeToFile();
int crossedZero;


void randomInteger(){
  srand((unsigned)time(0));
  for (int i =0; i < 100; i++){
          list.push_back(fRand(-1,1));
          if(i>0){
            list[i] = list[i-1] + list[i];
          }
          if(list[i] > 0 && list[i-1] < 0){
            crossedZero += 1;
          }else if(list[i] < 0 && list[i-1] > 0){
            crossedZero += 1;
          }

      }
}



double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

void writeToFile(){
  ofstream output_file("./data.txt");
  ostream_iterator<double> output_iterator(output_file, "\n");
  copy(list.begin(), list.end(), output_iterator);
}



int main() {

  randomInteger();
  writeToFile();

  //vector<int> v { 34,23 };

  //cout << v[0] << endl;


}
