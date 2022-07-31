#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <string>
#include <cstdlib>
using namespace std;

int trajectory(int Nfiles, const vector<int>& np_tc, int N, int I){
  int NC = np_tc.size();
  string line;
  char name[30];
  double Dt;
  
  vector<vector<double> > r_dat;
  r_dat.resize(Nfiles);
  for(int t = 0; t < Nfiles; ++t){
    r_dat[t].resize(3);
  }

  for(int t = 0; t < Nfiles; ++t){
    sprintf(name, "pos_vel/pos_vel_%d.dat", t);
    ifstream Con(name);
    if(! Con){
      Con.close();
      // //
      // cout<<name<<endl;
      // //
      return 1;
    }
    for(int m = 0; m < NC; ++m){
      int MM = np_tc[m];
      for(int i = 0; i < MM; ++i){
  	if(m == N && i == I){
  	  for(int k = 0; k < 3; ++k)
  	    Con >> r_dat[t][k];
  	}else{
  	  for(int k = 0; k < 3; ++k)
  	    Con >> line;
  	}
  	for(int k = 0; k < 3; ++k)
  	  Con >> line;
      }
    }
    if(t == 0)
      Con >> Dt;    
    Con.close();
  }
  
  // //
  // cout<<Dt<<endl;
  // //

  sprintf(name, "trajectory.dat");
  ofstream Out(name);
  if(!Out){
    Out.close();
    return 2;
  }
  for(int t = 0; t < Nfiles; ++t){
    Out << t*Dt << " ";
    for(int k = 0; k < 3; ++k)
      Out << r_dat[t][k] << " ";
    Out << endl;
  }
  Out.close();
  return 0;
}


int main(int argc, char* argv[]){
  //argv[0] is the program namespace
  int N = atoi(argv[1]);
  int I = atoi(argv[2]);

  int Ncomp;
  vector<int> np_tc;
  int mx_tcf;
  int Nfiles;

  ifstream IN("tc_input.txt");
  if(! IN){
    cout<<"Could not read tc_input.txt\n";
    IN.close();
    return 1;
  }
  IN >> Ncomp;
  for(int n = 0; n < Ncomp; ++n){
    int sp;
    IN >> sp;
    np_tc.push_back(sp);
  }
  IN >> mx_tcf;
  IN >> Nfiles;
  IN.close();

  //Num. of particles taken into account for computing time correlations
  for(int i = 0; i < Ncomp; ++i){
    if(np_tc[i] > mx_tcf) np_tc[i] = mx_tcf;
  }

  trajectory(Nfiles, np_tc, N, I);
  
  return 0;
}
