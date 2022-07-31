//Nummerically integrate a tabulated function
#include <iostream>
#include <fstream>
#include <deque>
#include <string> 
#include <cstdlib>
#include <cmath>
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
    // //
    // cout<<pair[0]<<","<<pair[1]<<endl;
    // //
  }
  In.close();
  

  deque<deque<double> > pmf;
  int size = data.size();
  for(int i = 0; i < size; ++i){
    //cout << data[i][0] << " " << data[i][1] << endl; 
    double x = data[i][0]; 
    double y = 1.e6;
    if(data[i][1] > 0.) 
      y = -log(data[i][1]);
    double f = 0.;
    if(i > 1 && (x-pmf[i-1][0])>1.e-6)
      f = -(y-pmf[i-1][1])/(x-pmf[i-1][0]);
    
    cout << f << endl;
    pair.clear();
    pair.push_back(x);
    pair.push_back(y);
    pair.push_back(f);
    pmf.push_back(pair);
  }
  
  ofstream Out("pmf.dat");
  if(! Out){
    cout << "problem opening pmf.dat" << endl;
    Out.close();
    return false;
  }
  for(int i = 0; i < size; ++i){
    Out << pmf[i][0] << " " << pmf[i][1] << " " << pmf[i][2] << "\n";
  }
  Out.close();
  
  return 0;
}
