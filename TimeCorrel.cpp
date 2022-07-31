#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <string>
using namespace std;

int diffusion(int Nfiles, const vector<int>& np_tc){ //v autocorrel is non-normalized!
  int NC = np_tc.size();
  string line;
  char name[30];
  
  vector<vector<vector<double> > > r0, v0;
  r0.resize(NC);
  v0.resize(NC);
  for(int i = 0; i < NC; ++i){
    int NN = np_tc[i];
    r0[i].resize(NN);
    v0[i].resize(NN);
    for(int j = 0; j < NN; ++j){
      r0[i][j].resize(3);
      v0[i][j].resize(3);
    }
  }
  double rt[3], vt[3], d[3], Dt;
  
  vector<vector<int> > norm(NC, vector<int>(Nfiles+1, 0));
  vector<vector<vector<double> > > msd(NC);
  vector<vector<vector<double> > > vacf(NC);
  for(int i = 0; i < NC; ++i){
    msd[i].resize(Nfiles+1);
    vacf[i].resize(Nfiles+1);
    for(int j = 0; j < Nfiles+1; ++j){
      msd[i][j].resize(3);
      vacf[i][j].resize(3);
      for(int k = 0; k < 3; ++k){
	msd[i][j][k] = 0.;
	vacf[i][j][k] = 0.;
      }
    }
  }
  
  for(int t0 = 0; t0 < Nfiles; ++t0){
    sprintf(name, "pos_vel/pos_vel_%d.dat", t0);
    std::ifstream Con(name);
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
	for(int k = 0; k < 3; ++k){
	  Con >> r0[m][i][k];          
	}                              
	for(int k = 0; k < 3; ++k){
	  Con >> v0[m][i][k];         
	}
	// for(int k = 0; k < 7; ++k){
	//   Con >> line;
	// }
      }
    }
    if(t0 == 0)
      Con >> Dt;    
    Con.close();
    
    //autocorrelation Dt=0
    for(int n = 0; n < NC; ++n){
      int NN = np_tc[n];
      for(int j = 0; j < NN; ++j){
	for(int k = 0; k < 3; ++k){
	  d[k] = 0.;
	  msd[n][0][k] += 0.;
	  vacf[n][0][k] += v0[n][j][k]*v0[n][j][k];
	}
	norm[n][0]++;
      }
    }
    //correlation with other times
    for(int t = t0+1; t < Nfiles; ++t){
      sprintf(name, "pos_vel/pos_vel_%d.dat", t);
      std::ifstream Con1(name);
      if(! Con1){
	Con1.close();
	// //
	// cout<<name<<endl;
	// //
	return 1;
      }
      
      for(int n = 0; n < NC; ++n){
	int NN = np_tc[n];
	for(int j = 0; j < NN; ++j){
	  for(int l = 0; l < 3; ++l)
	    Con1 >> rt[l];
	  for(int l = 0; l < 3; ++l)
	    Con1 >> vt[l];
	  // for(int k = 0; k < 7; ++k)
	  //   Con1 >> line;
	  for(int k = 0; k < 3; ++k){
	    d[k] = rt[k]-r0[n][j][k];
	    msd[n][t-t0][k] += d[k]*d[k];
	    vacf[n][t-t0][k] += v0[n][j][k]*vt[k];
	  }
	  norm[n][t-t0]++;
	}
      }
      Con1.close();
    }
  }
  
  
  sprintf(name, "difudata.dat");
  std::ofstream Out(name);
  if(! Out){
    Out.close();
    return 2;
  }
  Out << setiosflags(ios::fixed) << setprecision(6);
  for(int t = 0; t < Nfiles; ++t){
    Out << t*Dt << " |";
    for(int n = 0; n < NC; ++n){ 
      Out << "|r2 ";
      for(int k = 0; k < 3; ++k)
	Out << msd[n][t][k]/norm[n][t] << " ";
      Out << "|vv ";
      for(int k = 0; k < 3; ++k)
	Out << vacf[n][t][k]/norm[n][t] << " "; 
      Out << "|";
    }
    Out << endl;
  }
  Out.close();
  
  return 0;
}


int main(){
  int Ncomp;
  vector<int> np_tc;
  int mx_tcf;
  int Nfiles;
  // int Ndat;
  // double avP;
  // double V;

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
  // IN >> Ndat;
  // IN >> avP;
  // IN >> V;
  IN.close();

  
  //   cout << "Enter number of species: ";
  //   cin >> Ncomp;
  //   
  //   cout << "Enter number of particle per species.\n";
  //   for(int n = 0; n < Ncomp; ++n){
  //     int sp;
  //     cout << "Species " << n << ": ";
  //     cin >> sp;
  //     np_tc.push_back(sp);
  //   }
  //   
  //   cout << "Enter the maximum number of particles per species taken into account for the diffusion computation: ";
  //   cin >> mx_tcf;
  //   
  //   cout << "Enter the number of files to read for computing diffusion: ";
  //   cin >> Nfiles;
  //   
  //   cout << "Enter the number of rows in ptens.dat: ";
  //   cin >> Ndat;
  //   
  //   cout << "Enter the average pressure in the system: ";
  //   cin >> avP;
  //   
  //   cout << "Finally.. enter the volume of the system: ";
  //   cin >> V;

  //Num. of particles taken into account for computing time correlations
  for(int i = 0; i < Ncomp; ++i){
    if(np_tc[i] > mx_tcf) np_tc[i] = mx_tcf;
  }

  //error
  ofstream Err("ErrorTC.dat");
  Err.close();
  
  // if(! shear_visc(Ndat, V, avP)){
  //   Err.open("ErrorTC.dat", ios::app);
  //   Err << "problem opening viscodata.dat\n"; 
  //   Err.close();
  //   return 1;
  // }
  
  Err.open("ErrorTC.dat", ios::app);
  Err << "Computing diffusion graphs...\n";
  Err.close();
  int code = diffusion(Nfiles, np_tc);
  if(code != 0){
    Err.open("ErrorTC.dat", ios::app);
    if(code == 1)
      Err << "problem opening pos_vel.dat\n";
    else if(code == 2)
      Err << "problem opening difudata.dat\n"; 
    Err.close();
    return 1;
  }
  Err.open("ErrorTC.dat", ios::app);
  Err << "Diffusion graphs ready!\n";
  Err.close();
  
  return 0;
}
