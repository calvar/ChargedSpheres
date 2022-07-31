//Nummerically integrate a tabulated function
#include <iostream>
#include <fstream>
#include <deque>
#include <string> 
#include <cstdlib>
using namespace std;

int main(int argc, char* argv[]){
  int columnX = atoi(argv[2]);//argv[0] is the program namespace
  int columnY = atoi(argv[3]);
  deque<deque<double> > data;
  deque<double> pair;
  
  string dummy;
  ifstream In(argv[1]);
  if(! In){
    In.close();
    cout << "problem opening " << argv[1] << endl;
    return 1;
  }
  int count = 0;
  while(!In.eof()){
    if(argc == 5){
      if(count > atoi(argv[4])) break;
    }
    pair.clear();
    for(int i = 0; i < columnX; ++i)
      In >> dummy;
    pair.push_back(atof(dummy.c_str()));
    for(int i = 0; i < columnY-columnX; ++i)
      In >> dummy;
    pair.push_back(atof(dummy.c_str()));
    if(In.eof()) break;
    In.ignore(1023,'\n');
    data.push_back(pair);
    count++;
    //
    cout<<pair[0]<<","<<pair[1]<<endl;
    //
  }
  In.close();
  
  double integral = 0;
  int size = data.size();
  for(int i = 1; i < size; ++i){
    //cout << data[i][0] << " " << data[i][1] << endl;
    double x1 = data[i-1][0]; double y1 = data[i-1][1];
    double x2 = data[i][0]; double y2 = data[i][1];
    double dx = x2 - x1;
    double dy = y2 - y1;
    integral += dy * (x2*x2-x1*x1) / (2*dx) + y1*dx - x1*dy;
  }
  
  cout << integral << endl;
  
  return 0;
}
