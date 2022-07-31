
#include <cstring>
#include <complex>
#include <deque>
#include <unistd.h>
#include <sys/stat.h>
#include "functions.hpp"
using namespace std;

#define OMP_THR 6

//Read input
bool r_input(int& Nsteps, double& pre_time, int& iprint, int& signal, int& Nkinds, vector<int>& N, vector<double>& m, vector<double>& r, vector<double>& q, double& dens, double& Dt, double& ArLx, double& ArLy, double& temp, double& Q, double& gamma, double& sr_cut, int& K_max, double& alpha, double& bar_w, double& corrf, double& tcfd, int& mx_tcf){
  N.clear();
  m.clear();
  r.clear();
  q.clear();
  string line;
  ifstream In("input.txt");
  if(! In){
    ofstream Err("Error.dat");
    Err << "Couldn't open " << "input.txt" << endl;
    Err.close();
    In.close();
    return false;
  }
  In >> line >> line >> line >> Nsteps;
  In >> line >> line >> line >> pre_time;
  In >> line >> line >> iprint;
  In >> line >> signal;
  In >> line >> line >> line >> Nkinds;
  In >> line >> line >> line;
  for(int i = 0; i < Nkinds; ++i){
    int a;
    In >> a;
    N.push_back(a);
  }
  In >> line;
  for(int i = 0; i < Nkinds; ++i){
    double a;
    In >> a;
    m.push_back(a);
  }
  In >> line;
  for(int i = 0; i < Nkinds; ++i){
    double a;
    In >> a;
    r.push_back(a);
  }
  In >> line;
  for(int i = 0; i < Nkinds; ++i){
    double a;
    In >> a;
    q.push_back(a);
  }
  In >> line >> dens;
  In >> line >> line >> line >> Dt;
  In >> line >> line >> line >> line >> line;
  In >> ArLx >> ArLy; 
  In >> line >> temp;
  In >> line >> Q;
  In >> line >> gamma;
  In >> line >> line >> sr_cut;
  In >> line >> line >> line >> K_max;
  In >> line >> line >> alpha;
  In >> line >> line >> line >> bar_w;
  In >> line >> line >> line >> corrf;
  In >> line >> line >> tcfd;
  In >> line >> line >> line >> line >> mx_tcf;
  In.close();
  return true;
}


//Test for overlaps
bool overlap_test(const Particles& part, const double L[3], double frac){
  int Tpart = 0;
  for(int i = 0; i < part.get_Nkinds(); ++i)
    Tpart += part.get_N(i);
  double srt = pow(2., 1./6);
  double ri[3], rj[3];
  
  for(int i = 0; i < Tpart-1; ++i){
    int m, mi;
    part.get_kind(i, m, mi);
    for(int a = 0; a < 3; ++a)
      ri[a] = part.get_pos(m, mi, a);

    for(int j = i+1; j < Tpart; ++j){
      int n, nj;
      part.get_kind(j, n, nj);
      for(int a = 0; a < 3; ++a)
	rj[a] = part.get_pos(n, nj, a);
      
      double rij[3];
      for(int a = 0; a < 3; ++a){
	rij[a] = ri[a] - rj[a];
	rij[a] -= L[a] * floor(rij[a] / L[a] + 0.5);
      }
      double RIJSQ;
      dot(rij, rij, RIJSQ);
      double RIJ = sqrt(RIJSQ);

      double sigma = part.get_rad(m) + part.get_rad(n);
      double sig_sq = sigma * sigma;
      double rsh = srt * sigma;

      if(RIJ < frac * sigma){ //check for overlaps
	double p[3];
	ofstream Err("Error.dat", ios::app);
	Err << "sigma= "<<sigma<<" rij="<<RIJ<<"\n";
	for(int a = 0; a < 3; ++a)
	  p[a] = part.get_mom(m, mi, a);
	Err <<"r "<<m<<","<<mi<<" =("<<ri[0]<<","<<ri[1]<<","<<ri[2]<<") "
	    <<"p "<<m<<","<<mi<<" =("<<p[0]<<","<<p[1]<<","<<p[2]<<")\n";
	for(int a = 0; a < 3; ++a)
	  p[a] = part.get_mom(n, nj, a);
	Err <<"r "<<n<<","<<nj<<" =("<<rj[0]<<","<<rj[1]<<","<<rj[2]<<") "
	    <<"p "<<n<<","<<nj<<" =("<<p[0]<<","<<p[1]<<","<<p[2]<<")\n";
	Err.close();
	return true;
      }
    }
  }
  return false;
}


//Initial configuration
int ini_pos(Particles& part, const double L[3]){
  int NC = part.get_Nkinds();                                  
  vector<double> rad(NC);
  for(int n = 0; n < NC; ++n){
    rad[n] = part.get_rad(n);
  }

  double r[3] = {0., 0., 0.};

  //place the particles
  double shift_m = 1. + 1.e-3; 
  
  double sq = sqrt(3.) / 2.;
  double sqt = 1./ (2 * sqrt(3.));
  double sqtt = sqrt(6.) / 3;
  for(int n = 0; n < NC; ++n){
    int NN = part.get_N(n);
    int num1 = static_cast<int>(max(floor((L[0]-1.e-2)/(2*rad[n] * shift_m)), 1.));
    int num2 = static_cast<int>(max(floor((L[1]-1.e-2)/(2*rad[n]*sq * shift_m)), 1.));
    int num3 = static_cast<int>(max(floor((L[2]-1.e-2)/(2*rad[n]*sqtt * shift_m)), 1.));
    int cont = 0;
    while(cont < NN){
      bool escape = false;
      for(int nz = 0; nz < num3; ++nz){
	for(int ny = 0; ny < num2; ++ny){
	  for(int nx = 0; nx < num1; ++nx){
	    // //
	    // cout<<cont<<endl;
	    // //
	    if(cont < NN){
	      if(fmod(static_cast<double>(nz),2.) == 0){
		if(fmod(static_cast<double>(ny),2.) == 0){
		  r[0] = 2 * rad[n] * static_cast<double>(nx) /* shift_m*/;
		}else{
		  r[0] = 2 * rad[n] * (static_cast<double>(nx) /* shift_m*/ + 0.5);
		}
		r[1] = 2 * rad[n] * static_cast<double>(ny) * sq /* shift_m*/;
	      }else{
		if(fmod(static_cast<double>(ny),2.) != 0){
		  r[0] = 2 * rad[n] * static_cast<double>(nx) /* shift_m*/;
		}else{
		  r[0] = 2 * rad[n] * (static_cast<double>(nx) /* shift_m*/ + 0.5);
		}
		r[1] = 2 * rad[n] * (static_cast<double>(ny) * sq + sqt) /* shift_m*/;
	      }
	      r[2] = 2 * rad[n] * static_cast<double>(nz) * sqtt /* shift_m*/;
	      
	      for(int a = 0; a < 3; ++a)
		r[a] -= L[a] * floor(r[a] / L[a] + 0.5);
	
	      for(int a = 0; a < 3; ++a)
		part.set_pos(n, cont, a, r[a]);
	      // //
	      // cout<<"Done. Checking for overlaps...\n";
	      // //
	      
	      if(n > 0){
		int con = 0;
		double rj[3];
		for(int m = 0; m < n; ++m){
		  int MM = part.get_N(m);
		  double rr = rad[n] + rad[m];
		  double rsq = rr * rr;
		  bool test = false;
		  for(int j = 0; j < MM; ++j){
		    for(int a = 0; a < 3; ++a)
		      rj[a] = part.get_pos(m, j, a);
		    
		    double XIJ = r[0] - rj[0];
		    XIJ -= L[0] * floor(XIJ / L[0] + 0.5);
		    double YIJ = r[1] - rj[1];
		    YIJ -= L[1] * floor(YIJ / L[1] + 0.5);
		    double ZIJ = r[2] - rj[2];
		    ZIJ -= L[2] * floor(ZIJ / L[2] + 0.5);
		    double RIJSQ = XIJ * XIJ + YIJ * YIJ + ZIJ * ZIJ;
		    
		    if(RIJSQ < rsq){
		      test = true;
		      break;
		    }
		  }
		  if(test){
		    break;
		  }
		  if(m == n-1){
		    cont++;
		  }
		  con++;
		  // //
		  // cout<<"Ovrlp="<<test<<endl;
		  // //
		}
	      }else{
		cont++;
	      }
	      // //
	      // cout<<"Done checking.\n";
	      // //
	    }else{
	      escape = true;
	    }

	    // //
	    // cout<<"esc "<<escape<<endl;
	    // //
	    
	    if(escape)
	      break;
	  }
	  if(escape)
	    break;
	}
	if(escape)
	  break;
      }
    }
  }
 
  //print
  ofstream InPos("iniPos.dat");
  if(! InPos){
    ofstream Err("Error.dat", ios::app);
    Err << "Couldn't open " << "iniPos.dat" << endl;
    Err.close();
    InPos.close();
    return -1;
  }
  for(int n = 0; n < NC; ++n) //charges
    InPos << "0 ";
  InPos << endl;
  InPos << L[0] << " " << L[1] << " " << L[2] << endl;
  
  for(int n = 0; n < NC; ++n){
    int NN = part.get_N(n);
    for(int i = 0; i < NN; i++){
      for(int a = 0; a < 3; ++a)
	r[a] = part.get_pos(n, i, a);
      
      InPos << r[0] << " " << r[1] << " " << r[2];
      InPos << " 0 " << "0 " << "0 " << endl; 
    }
  }
  InPos.close();

  //test for overlaps
  if(overlap_test(part, L, 1.)){
    cout<<"overlap at ini_pos!"<<endl;
    return 1;
  }
  
  return 0;
}

void ini_mom(Particles& part, double temp, const gsl_rng* ran){
  int NC = part.get_Nkinds();
  int Ntot = 0;
  for(int i = 0; i < NC; ++i)
    Ntot += part.get_N(i);
  
  //for translation velocities************************************************
  vector<vector<vector<double> > > rnd(3,vector<vector<double> >(NC, vector<double>(1)));
  //rnd[0].resize(NC); rnd[1].resize(NC); rnd[2].resize(NC); 
  for(int n = 0; n < NC; ++n){
    for(int j = 0; j < 3; ++j)
      rnd[j][n].resize(part.get_N(n)+1);
  }
  for(int n = 0; n < NC; ++n){
    int NN = part.get_N(n);
    for(int i = 0; i < NN; ++i){
      for(int j = 0; j < 3; ++j)
	rnd[j][n][i] = gsl_rng_uniform(ran);
    }
  }
  for(int n = 0; n < NC; ++n){
    int NN = part.get_N(n);
    double mtemp = temp * part.get_mass(n) * 1.e-6;
    for(int i = 0; i < NN; i += 2){
      if((fmod(static_cast<float>(NN),2) != 0.) && (i == NN-1)){
	bool test = true;
	while(test){
	  for(int j = 0; j < 3; j++){
	    double rr = gsl_rng_uniform(ran);
	    gauss(mtemp, rr, rnd[j][n][i]);
	    if(rnd[j][n][i] > -1 && rnd[j][n][i] < 1){
	      test = false;
	    }else if(rr > -1 && rr < 1){
	      rnd[j][n][i] = rr;
	      test = false;
	    }
	  }
	}
      }else{
	bool test = true;
	while(test){
	  for(int j = 0; j < 3; j++){
	    gauss(mtemp, rnd[j][n][i], rnd[j][n][i+1]);
	    if(rnd[j][n][i] > -1 && rnd[j][n][i] < 1 && rnd[j][n][i+1] > -1 && rnd[j][n][i+1] < 1){
	      test = false;
	    }
	  }
	}
      }
    }
  }

  double P[3] = {0., 0., 0.};
  //double pp[3][Ntot+1];
  vector<vector<double> > pp(3, vector<double>(Ntot+1));
  
  int cnt = 0;
  for(int n = 0; n < NC; ++n){
    int NN = part.get_N(n);
    for(int i = 0; i < NN; ++i){
      for(int j = 0; j < 3; j++){
	pp[j][cnt] = rnd[j][n][i];
	P[j] += pp[j][cnt];
      }
      cnt++;
    }
  }

  double pr = 1.;
  
//   //
//   cout<<"begin linear X..\n";
//   //
  
 //  // X part/////////////////////////////////////////////////
  int II = 0;/*, nn=0, ii=0;*/
  long scont = 0;
  while(fabs(P[0]) > pr){ //NOTE: in fabsP > x, if the value of x is very 
                           // small the program will have an error. 
                           //(what is small depends on Ntot) 
                           // for N=256 x=0.1 is fine.
    while(P[0] > pr){
      if(pp[0][II] > 0){ 
	pp[0][II] = -pp[0][II];
// 	rev_II(NP, nn, ii, II);
// 	rnd[0][nn][ii] *= -1;
	P[0] += 2 * pp[0][II];
      }
      II++;
      if(II == Ntot - 1) II = 0;
    }
    while(P[0] < -pr){
      if(pp[0][II] < 0){ 
	pp[0][II] = -pp[0][II];
// 	rev_II(NP, nn, ii, II);
// 	rnd[0][nn][ii] *= -1;
	P[0] += 2 * pp[0][II];
      }
      II++;
      if(II == Ntot - 1) II = 0;
    }
    if(scont > 1e6) break;
    scont++;
  }
  P[0] -= pp[0][Ntot-1];
  pp[0][Ntot-1] = - P[0];
//   rev_II(NP, nn, ii, Ntot-1);
//   rnd[0][nn][ii] = pp[0][Ntot-1];
  P[0] = 0.;
  cnt = 0;
  for(int n = 0; n < NC; ++n){
    int NN = part.get_N(n);
    for(int i = 0; i < NN; ++i){
      rnd[0][n][i] = pp[0][cnt];
      P[0] += rnd[0][n][i];
      cnt++;
    }
  }
  if(fabs(P[0]) > 1e-12){
    ofstream Err("Error.dat", ios::app);
    Err << "Px!=0 !!\n";
    Err.close();
  }
  
//   //
//   cout<<"begin linear Y..\n";
//   //
  
  //  // Y part/////////////////////////////////////////////////
  II = 0;
  scont = 0;
  while(fabs(P[1]) > pr){ //NOTE: in fabsP > x, if the value of x is very 
                           // small the program will have an error. 
                           //(what is small depends on Ntot) 
                           // for N=256 x=0.1 is fine.
    while(P[1] > pr){
      if(pp[1][II] > 0){ 
	pp[1][II] = -pp[1][II];
// 	rev_II(NP, nn, ii, II);
// 	rnd[1][nn][ii] *= -1;
	P[1] += 2 * pp[1][II];
      }
      II++;
      if(II == Ntot - 1) II = 0;
    }
    while(P[1] < -pr){
      if(pp[1][II] < 0){ 
	pp[1][II] = -pp[1][II];
// 	rev_II(NP, nn, ii, II);
// 	rnd[1][nn][ii] *= -1;
	P[1] += 2 * pp[1][II];
      }
      II++;
      if(II == Ntot - 1) II = 0;
    }
    if(scont > 1e6) break;
    scont++;
  }
  P[1] -= pp[1][Ntot-1];
  pp[1][Ntot-1] = - P[1];
//   rev_II(NP, nn, ii, Ntot-1);
//   rnd[1][nn][ii] = pp[1][Ntot-1];
  P[1] = 0.;
  cnt = 0;
  for(int n = 0; n < NC; ++n){
    int NN = part.get_N(n);
    for(int i = 0; i < NN; ++i){
      rnd[1][n][i] = pp[1][cnt];
      P[1] += rnd[1][n][i];
      cnt++;
    }
  }
  if(fabs(P[1]) > 1e-12){
    ofstream Err("Error.dat", ios::app);
    Err << "Py!=0 !!\n";
    Err.close();
  }

//   //
//   cout<<"begin linear Z..\n";
//   //
  
  //  // Z part/////////////////////////////////////////////////
  II = 0;
  scont = 0;
  while(fabs(P[2]) > pr){ //NOTE: in fabsP > x, if the value of x is very 
                           // small the program will have an error. 
                           //(what is small depends on Ntot) 
                           // for N=256 x=0.1 is fine.
    while(P[2] > pr){
      if(pp[2][II] > 0){ 
	pp[2][II] = -pp[2][II];
// 	rev_II(NP, nn, ii, II);
// 	rnd[2][nn][ii] *= -1;
	P[2] += 2 * pp[2][II];
      }
      II++;
      if(II == Ntot - 1) II = 0;
    }
    while(P[2] < -pr){
      if(pp[2][II] < 0){ 
	pp[2][II] = -pp[2][II];
// 	rev_II(NP, nn, ii, II);
// 	rnd[2][nn][ii] *= -1;
	P[2] += 2 * pp[2][II];
      }
      II++;
      if(II == Ntot - 1) II = 0;
    }
    if(scont > 1e6) break;
    scont++;
  }
  P[2] -= pp[2][Ntot-1];
  pp[2][Ntot-1] = - P[2];
//   rev_II(NP, nn, ii, Ntot-1);
//   rnd[2][nn][ii] = pp[2][Ntot-1];
  P[2] = 0;
  cnt = 0;
  for(int n = 0; n < NC; ++n){
    int NN = part.get_N(n);
    for(int i = 0; i < NN; ++i){
      rnd[2][n][i] = pp[2][cnt];
      P[2] += rnd[2][n][i];
      cnt++;
    }
  }
  if(fabs(P[2]) > 1e-12){
    ofstream Err("Error.dat", ios::app);
    Err << "Pz!=0 !!\n";
    Err.close();
  }

//   ofstream InVel("iniVel.dat");
//   if(! InVel){
//     ofstream Err("Error.dat", ios::app);
//     Err << "Couldn't open " << "iniVel.dat" << endl;
//     Err.close();
//     InVel.close();
//     return;
//   }
  
  double ptot[3] = {0,0,0};
  for(int n = 0; n < NC; ++n){
    int NN = part.get_N(n);
    for(int i = 0; i < NN; ++i){
      for(int k = 0; k < 3; ++k){
	part.set_mom(n, i, k, rnd[k][n][i]);
  	ptot[k] += part.get_mom(n, i, k);
      }
    }
  }
  // if(ptot[0] > 1e-8 || ptot[2] > 1e-8 || ptot[2] > 1e-8){
  //   Err.open("Error.dat", ios::app);
  //   Err << "Total momentum is not zero! : "
  // 	<<ptot[0]<<" , "<<ptot[1]<<" , "<<ptot[2]<<endl;
  //   Err.close();
  //   gsl_rng_free(ran);
  //    return 1;
  // }
//   InVel.close();
}

//Box-Muller algorithm
void gauss(double variance, double& rnd1, double& rnd2){
  double rndA, rndB;
  rndA = sqrt(- 2 * log(rnd1)) * cos(2 * M_PI * rnd2);
  rndB = sqrt(- 2 * log(rnd1)) * sin(2 * M_PI * rnd2);
  double sqrtVar = sqrt(variance);
  rnd1 = rndA * sqrtVar;
  rnd2 = rndB * sqrtVar;
}

//Print configuration
bool print_conf(const Particles& part, const Thermostat& thermo, const double L[3], double dt, const int scal_cnt[10], bool fin){
  int NC = part.get_Nkinds();
  char name[30];
  sprintf(name, "conf.dat");
  std::ofstream Con(name);
  if(! Con){
    Con.close();
    return false;
  }
  Con << setiosflags(ios::fixed) << setprecision(12);
  Con << dt << endl;  //time step
  for(int i = 0; i <10; ++i)  //couters for dt scaling
    Con << scal_cnt[i] << " ";
  Con << endl;
  for(int n = 0; n < NC; ++n) //charges
    Con << part.get_charge(n) << " ";
  Con << endl;
  Con << L[0] << " " << L[1] << " " << L[2] << endl;
  for(int n = 0; n < NC; ++n){
    int NN = part.get_N(n);
    for(int i = 0; i < NN; i++){
      for(int a = 0; a < 3; ++a)
	Con << part.get_pos(n, i, a) << " ";
      for(int a = 0; a < 3; ++a)
	Con << part.get_mom(n, i, a) << " ";
      Con << endl; 
    }
  }
  Con << thermo.get_coo() << " " << thermo.get_strn() << endl;
  Con.close();

  //print state of random generators//
  if(! thermo.print_ran(0, fin)){
    ofstream Err("Error.dat", ios::app);
    Err << "problem printing ran_state" << 0 << ".dat\n"; 
    Err.close();
    return false;
  }
// 	//
// 	cout<<" t"<<thermo.get_coo()<<" ";
// 	//

  return true;
}

//Read configuration
bool chrg_conf(Particles& part, Thermostat& thermo, double L[3], double& dt, int scal_cnt[10]){
  std::string line;
  int NC = part.get_Nkinds();
  char name[30];
  sprintf(name, "conf.dat");
  std::ifstream Con(name);
  if(! Con){
    Con.close();
    return false;
  }

  Con >> dt;  //time step
  for(int i = 0; i < 10; ++i) //counters for dt scalling
    Con >> scal_cnt[i];

  for(int n = 0; n < NC; ++n){ //charges
    double q;
    Con >> q;
    part.set_charge(n, q);
  }
  Con >> L[0] >> L[1] >> L[2];
      
  for(int n = 0; n < NC; ++n){
    int NN = part.get_N(n);
    for(int i = 0; i < NN; i++){
      for(int a = 0; a < 3; ++a){
	double p;
	Con >> p;
	part.set_pos(n, i, a, p); 
      }
      for(int a = 0; a < 3; ++a){
	double m;
	Con >> m;
	part.set_mom(n, i, a, m);
      }
    }
  }
  double s, cn;
  Con >> s >> cn;
  thermo.set_coo(s); 
  thermo.set_strn(cn);
  Con.close();

  //charge thermo rng state
  if(! thermo.read_ran(0)){
    ofstream Err("Error.dat", ios::app);
    Err << "problem reading ran_state" << 0 << ".dat\n"; 
    Err.close();
    return false;
  }

  return true;
}


void tabul(double *Tab, double alpha){
  double extL = 200 * 1.001;
  double alsq = alpha * alpha;
  double spi = sqrt(M_PI);
  for(int i = 1; i < TAB_SIZE+1; i++){
    double r =  0.5 * extL * i / TAB_SIZE;
    double rsq = r * r;
    double sr1 = 1.0 / r;
    double ERFC = erfc(alpha * r);
    double aer = 2.0 * alpha * exp(- alsq * rsq) / spi;
    Tab[i-1] = ERFC * sr1;
    Tab[TAB_SIZE+(i-1)] = Tab[i-1] + aer;
  }
  // ofstream Table("TabBC.dat");
//   for(int i = 0; i < TAB_SIZE; i++){
//     Table << 0.5*extL*(i+1)/TAB_SIZE << " " << Tab[i] << " "
//      << Tab[TAB_SIZE+i] << endl;
//   }
  // Table.close();
}

bool force_R(Particles& part, const double L[3], const double *Tab, double cut_d, double& U_sr, double& W_sr, double& U_c, double& W_c, double Dt, bool *reduce){
  int Tpart = 0;
  for(int i = 0; i < part.get_Nkinds(); ++i)
    Tpart += part.get_N(i);

  //assume there is no need for time step reduction or maintaining current reduction
  //reduce[0] refers to wether one needs to scale dt by 1/10
  //reduce[1] indicates wether one needs to scale dt further by 1/10 (so 1/100)
  //in case no reduction is necessary, reduce[2] indicates if current scaling
  // should be maintained for a longer number of steps
  reduce[0] = false; reduce[1] = false; reduce[2] = false;
  
  //set forces to 0
  part.set_zero_F();
  //short range
  double Usr = 0.;
  double Wsr = 0.;
  double srt = pow(2., 1./6);
  //Coulomb
  double Uc = 0.;
  double Wc = 0.;
  double Lh = 0.5 * (200*1.001) / TAB_SIZE;

  //  #pragma omp_set_num_threads(OMP_THR);
  bool abort = false;
  #pragma omp parallel default(shared)
  {
    int tasks2 = 16 * omp_get_num_threads();
    int chunk = static_cast<int>(floor(static_cast<float>(Tpart)/tasks2));
    #pragma omp for schedule(dynamic, chunk) reduction(+: Usr, Wsr, Uc, Wc) nowait
    for(int i = 0; i < Tpart; ++i){
      #pragma omp flush(abort)
      if(!abort){
        double ri[3], rj[3], f[3];
	double vi[3], vj[3];
	double qi, qj;
	double massi, massj;
	
	int m, mi;
	part.get_kind(i, m, mi);
	massi = part.get_mass(m);
	for(int a = 0; a < 3; ++a){
	  ri[a] = part.get_pos(m, mi, a);
	  vi[a] = part.get_mom(m, mi, a);
	  vi[a] /= massi;
	}
	qi = part.get_charge(m);

	//even out the charge among threads
	int mx = static_cast<int>(ceil(static_cast<float>(Tpart-1)/2));
	if(fmod(static_cast<float>(Tpart),2) == 0. && i >= Tpart/2)
	  mx = static_cast<int>(floor(static_cast<float>(Tpart-1)/2));

	int j = i+1 - Tpart*static_cast<int>(floor((i+1)/Tpart + 0.5));
	int cnt = 0;
	while(cnt < mx){
	  // //
	  // cout<<cnt<<" "<<mx<<" "<<i<<","<<j<<"\n";
	  // //
	  int n, nj;
	  part.get_kind(j, n, nj);
	  massj = part.get_mass(n);
	  for(int a = 0; a < 3; ++a){
	    rj[a] = part.get_pos(n, nj, a);
	    vj[a] = part.get_mom(n, nj, a);
	    vj[a] /= massj;
	  }
	  qj = part.get_charge(n);
      
	  double rij[3];
	  for(int a = 0; a < 3; ++a){
	    rij[a] = ri[a] - rj[a];
	    rij[a] -= L[a] * floor(rij[a] / L[a] + 0.5);
	  }

	  //cutoff distance for the real part (cube)
	  if(fabs(rij[0])<cut_d && fabs(rij[1])<cut_d && fabs(rij[2])<cut_d){
	    double RIJSQ;
	    dot(rij, rij, RIJSQ);
	    double RIJ = sqrt(RIJSQ);
	    
	    double sigma = part.get_rad(m) + part.get_rad(n);
	    double sig_sq = sigma * sigma;
	    double rsh = srt * sigma;


	    //Check for posible future overlaps----
	    //relative velocity
	    double vij[3], vr=0.;
	    for(int a = 0; a < 3; ++a){
	      vij[a] = vi[a] - vj[a];
	      vr += vij[a] * rij[a]/RIJ;
	    }
	    if(vr < 0.){ //if particles are approaching...

              // ////
              // if(-vr > 100.){
              //   double ff1[3], ff2[3];
	      // 	double f1=0.; 
	      // 	double f2=0.;
	      // 	for(int a = 0; a < 3; ++a){
	      // 	  ff1[a] = part.get_F(m,mi,a);
	      // 	  ff2[a] = part.get_F(n,nj,a);
	      // 	  f1 += ff1[a]*ff1[a];//rij[a]/RIJ;
	      // 	  f2 += ff2[a]*ff2[a];//rij[a]/RIJ;
              //   }
	      // 	printf("(%d,%d) vs (%d,%d) vij=%f f1ij=%f f2ij=%f\n",m,mi,n,nj,vr,sqrt(f1),sqrt(f2));
              // }
	      // ////

	      double min_d = RIJ - MX_OVRLP*sigma;
	      if(!reduce[0]){ //if no reason to reduce has been found yet
	  	double travel_d = -vr*Dt;
	  	if(TSTEP_MULT * travel_d > min_d)
	  	  reduce[0] = true;
	      }
	      if(reduce[0] && !reduce[1]){
	  	double travel_d = -vr*Dt/TIME_MULT;
	  	if(TSTEP_MULT * travel_d > min_d) //further reduction is needed
	  	  reduce[1] = true;
	      }

	      if(!reduce[0] && !reduce[2]){ //check if the reduction must me maintained(re-check at the end of this function)
	  	double travel_d = -vr*Dt*TIME_MULT;
	  	if(TSTEP_MULT * travel_d > min_d)
	  	  reduce[2] = true;
	      }
	    }
	    // //
	    // if(m==1 && mi==175 && n==2 && nj==12){
	    //   cout<<"["<<reduce[0]<<","<<reduce[1]<<","<<reduce[2]<<"]\n";
	    //   cout<<"rij ("<<rij[0]<<","<<rij[1]<<","<<rij[2]<<") ";
	    //   cout<<"vij ("<<vij[0]<<","<<vij[1]<<","<<vij[2]<<") ";
	    //   cout<<"r.v= "<<vr<<" RIJ= "<<RIJ<<" min= "<<MX_OVRLP*sigma<<endl;  
	    // }
	    // //
	    //-----------------------

	    
	    if(RIJ < MX_OVRLP * sigma){ //check for overlaps
	      double p[3];
	      ofstream Err("Error.dat", ios::app);
	      Err << "sigma= "<<sigma<<" rij="<<RIJ<<"\n";
	      for(int a = 0; a < 3; ++a)
	        p[a] = part.get_mom(m, mi, a);
	      Err <<"r "<<m<<","<<mi<<" =("<<ri[0]<<","<<ri[1]<<","<<ri[2]<<") "
	      	<<"p "<<m<<","<<mi<<" =("<<p[0]<<","<<p[1]<<","<<p[2]<<")\n";
	      for(int a = 0; a < 3; ++a)
	        p[a] = part.get_mom(n, nj, a);
	      Err <<"r "<<n<<","<<nj<<" =("<<rj[0]<<","<<rj[1]<<","<<rj[2]<<") "
	      	<<"p "<<n<<","<<nj<<" =("<<p[0]<<","<<p[1]<<","<<p[2]<<")\n";
	      Err.close();
	      abort = true;
              #pragma omp flush(abort) //only safe with bool var
	    }
	      
	    //Lennard-Jones shifted potential, virial and force----------------------
	    double SR2s = sig_sq / RIJSQ;
	    double SR6s = SR2s * SR2s * SR2s;
	    double pot = 4 * SR6s * (SR6s - 1.);
	    double UIJ = 0.;
	    double WIJ = 0.;
	    if(RIJ < rsh){
	      UIJ = pot + 1.;//shifted potential (epsilon = 1)
	      WIJ = 6 * (pot + 4*SR6s*SR6s);
	    }
	    // ////
	    // cout<<" ("<<i<<","<<j<<") "<<UIJ<<","<<RIJ<<","<<rsh<<" ";
	    // ////
	    Usr += UIJ;
	    Wsr += WIJ;

	    //short range force
	    double FIJsr = WIJ / RIJSQ;
	    for(int a = 0; a < 3; ++a){
	      f[a] = FIJsr * rij[a];
	      part.add_F(m, mi, a, f[a]);
	    }
	    for(int a = 0; a < 3; ++a){
	      f[a] *= -1;
	      part.add_F(n, nj, a, f[a]);
	    }
	    // //
	    // cout<<"("<<i<<","<<j<<") "<<f[0]<<","<<f[1]<<","<<f[2];
	    // cout<<" F "<<part.get_F(m, mi, 0)<<","<<part.get_F(m, mi, 1)<<","<<part.get_F(m, mi, 2)<<" | ";
	    // cout<<part.get_F(n, nj, 0)<<","<<part.get_F(n, nj, 1)<<","<<part.get_F(n, nj, 2)<<"\n";
	    // //
	    //----------------------------------------------------------------------

	    double qq = qi * qj;
	    //Coulomb potential, virial and force-----------------------------------
	    if(fabs(qq) > 1.e-12){
	      int ri = static_cast<int>(floor(RIJ / Lh));
	      double R_L0 = RIJ - Lh * ri;
	      double R_L1 = RIJ - Lh * (ri+1);
	      double A = (Tab[ri] * R_L0 - Tab[ri-1] * R_L1) / Lh;
	      double B = (Tab[ri+TAB_SIZE] * R_L0 - Tab[(ri-1)+TAB_SIZE] * R_L1) / Lh;
	      UIJ = qq * A;
	      WIJ = qq * B;
	      Uc += UIJ;
	      Wc += WIJ;
	      // //
	      // cout<<i<<","<<j<<" "<<Tab[ri][0]<<" "<<Tab[ri-1][0]<<" "<<Tab[ri][1]<<" "<<Tab[ri-1][1]<<endl;
	      // //
	      
	      //Coulomb force
	      for(int a = 0; a < 3; ++a){
	      	f[a] = WIJ * rij[a] / RIJSQ;
	      	part.add_F(m, mi, a, f[a]);
	      }
	      for(int a = 0; a < 3; ++a){
	      	f[a] *= -1;
	      	part.add_F(n, nj, a, f[a]);
	      }
	    }
	    //----------------------------------------------------------------------
	    // //
	    // cout<<"("<<i<<","<<j<<") "<<f[0]<<","<<f[1]<<","<<f[2];
	    // cout<<" F "<<part.get_F(m, mi, 0)<<","<<part.get_F(m, mi, 1)<<","<<part.get_F(m, mi, 2)<<" | ";
	    // cout<<part.get_F(n, nj, 0)<<","<<part.get_F(n, nj, 1)<<","<<part.get_F(n, nj, 2)<<"\n";
	    // //
	    // ////
	    // double fc = 0.;
	    // for(int a = 0; a < 3; ++a)
	    //   fc += f[a]*f[a];
	    // if(sqrt(fc) > 50.){
	    //   printf("(%d,%d) vs (%d,%d) F=%f\n",m,mi,n,nj,sqrt(fc));
	    // }
	    // ////
	  }
	  j += 1 - Tpart*static_cast<int>(floor((j+1)/Tpart + 0.5));
	  cnt++;
        }
      }
    }
  }//end parallel
  if(abort)
    return false;

  U_sr = Usr;
  W_sr = Wsr;
  U_c = Uc;
  W_c = Wc;

  
  if(reduce[0]){
    if(reduce[2]) //maintain reduction?
      reduce[2] = false;
    //if(reduce[1])
    //  reduce[0] = false;
  }

  // //
  // cout<<reduce[0]<<" "<<reduce[1]<<" "<<reduce[2]<<endl;
  // // 
  
  ////
  //cout<<"no gpu: "<<Usr<<" "<<Wsr<<" "<<Uc<<" "<<Wc<<endl;
  // int NC = part.get_Nkinds();
  // for(int n = 0; n < NC; ++n){
  //   for(int i = 0; i < part.get_N(n); ++i){
  //     for(int a = 0; a < 3; ++a)
  // 	cout << part.get_F(n, i, a) << " ";
  //     cout << endl;
  //   }
  // }
  ////
  
  return true;
}

void wave_v(double *Kvec, const double L[3], double alpha, int K_max){
  double c = 4. * alpha * alpha;
  double P2[3] = {2.*M_PI/L[0], 2.*M_PI/L[1], 2.*M_PI/L[2]};
  int sz = K_max+1;
  int sz2 = sz*sz;
  for(int nz = 0; nz <= K_max; ++nz){
    double Kz = P2[2] * nz;
    for(int ny = 0; ny <= K_max; ++ny){
      double Ky = P2[1] * ny;
      for(int nx = 0; nx <= K_max; ++nx){
	double Kx = P2[0] * nx;
	double Ksq = Kx * Kx + Ky * Ky + Kz * Kz;
	if(Ksq != 0.){
	  Kvec[nz*sz2 + ny*sz + nx] = 4. * M_PI * exp(- Ksq / c) / Ksq; 
	}else{
	  Kvec[nz*sz2 + ny*sz + nx] = 0.;
	}
      }
    }
  }
  // ofstream Vec("Kvec.dat");
  // for(int z = 0; z <= K_max; ++z){
  //   for(int y = 0; y <= K_max; ++y){
  //     for(int x = 0; x <= K_max; ++x){
  // 	Vec << Kvec[z*sz2+y*sz+x] << " ";
  //     }
  //     Vec << "\n";
  //   }
  //   Vec << "\n";
  // }
  // Vec.close();
}


void force_K(Particles& part, const double L[3], int K_max, double alpha, const double *Kvec, double& U, double& W){
  int Tpart = 0;
  for(int i = 0; i < part.get_Nkinds(); ++i)
    Tpart += part.get_N(i);
  
//   double Lmax = max(L[2], max(L[1], L[0]));
  double V = L[0] * L[1] * L[2];
  double P2[3] = {2.*M_PI/L[0], 2.*M_PI/L[1], 2.*M_PI/L[2]};
//   double K_msq = 4*M_PI*M_PI*K_max*K_max/(Lmax*Lmax) * 1.001;
  K_max += 1;
  int K_max2 = K_max * K_max;
  int K_max3 = K_max2 * K_max;
  double Uc = 0.;
  double Wc = 0.;

  //  #pragma omp_set_num_threads(OMP_THR);
  #pragma omp parallel default(shared)
  {
    int tasks2 = 8 * omp_get_num_threads();
    int chunk = static_cast<int>(floor(static_cast<float>(K_max3)/tasks2));
    #pragma omp for schedule(dynamic, chunk) reduction(+: Uc, Wc) nowait
    for(int kn = 0; kn < K_max3; ++kn){
      int nz = static_cast<int>(floor(static_cast<float>(kn)/K_max2));
      float l = kn - K_max2 * nz;
      int ny = static_cast<int>(floor(l/K_max));
      int nx = static_cast<int>(l) - K_max * ny;
      double Kz = P2[2] * nz;
      double Ky = P2[1] * ny;
      double Kx = P2[0] * nx;
      double Ksq = Kx * Kx + Ky * Ky + Kz * Kz;
      
      if(Ksq != 0/* && Ksq < K_msq*/){
	//
	if(kn != (nz*K_max2+ny*K_max+nx)){
	  cout<<"kn != nz*K_max2+ny*K_max+nx ! ";
	  cout<<kn<<" != "<<nz*K_max2+ny*K_max+nx<<endl;
	}
	//
	double KK = Kvec[kn];
	double val = 2. * KK / V;  //mult by 2 for symmetry
	
	double K[4][3];
	K[0][0] = Kx; K[0][1] = Ky; K[0][2] = Kz;
	K[1][0] = -Kx; K[1][1] = Ky; K[1][2] = Kz;
	K[2][0] = Kx; K[2][1] = -Ky; K[2][2] = Kz;
	K[3][0] = -Kx; K[3][1] = -Ky; K[3][2] = Kz;
	
	//sum array
	complex<double> cc(0., 0.);
	complex<double> sum[8];
	for(int b = 0; b < 8; ++b)
	  sum[b] = cc;
	
	double self = 0.; //extract the self energy corresponding to each term (Kmax finite)
	
	//individual values array for each particle
	vector<vector<complex<double> > >  M(4, vector<complex<double> >(Tpart, cc));
	
	double ri[3], qi;
	
	//generate values
	for(int i = 0; i < Tpart; ++i){
	  int n, ni;
	  part.get_kind(i, n, ni);
	  for(int a = 0; a < 3; ++a)
	    ri[a] = part.get_pos(n, ni, a);
	  qi = part.get_charge(n);

	  if(qi != 0.){
	    double RIK[4];
	    complex<double> z[4];
	    for(int b = 0; b < 4; ++b){
	      dot(ri, K[b], RIK[b]);
	      z[b] = polar(1., RIK[b]);
	    }

	    for(int b = 0; b < 4; ++b){
	      M[b][i] = conj(z[b]);      
	      sum[b] += qi * M[b][i];	
	      sum[b+4] += qi * RIK[b] * z[b];
	    }
	    self += qi * qi;

	    // //
	    // for(int b = 0; b < 4; ++b)
	    //   cout<<"["<<b<<" "<<z[b]<<"] ";
	    // cout<<endl;
	    // //
	  }
	}


        // //
	// for(int l = 0; l < Tpart; ++l){
	//   for(int b = 0; b < 4; ++b)
	//     cout<<M[b][l]<<" ";
	//   cout<<endl;
	// }
	// cout<<endl;
	// //
	//cout<<"Kv="<<kn<<" "<<sum[2]<<endl;
	
	
	//compute energy - virial
	double Upar = 0, Wpar = 0;
	for(int b = 0; b < 4; ++b){
	  Upar += val * 0.5 * (norm(sum[b]) - self);
	  Wpar -= val * imag(sum[b+4] * sum[b]);
	}
	if((nx == 0 && ny == 0) || (nx == 0 && nz == 0) || (ny == 0 && nz == 0)){
	  Upar /= 4.;
	  Wpar /= 4.;
	}else if((nz == 0 && nx != 0 && ny != 0) || (ny == 0 && nz != 0 && nx != 0) || (nx == 0 && ny != 0 && nz != 0)){
	  Upar /= 2.;
	  Wpar /= 2.;
	}
	Uc += Upar;       
	Wc += Wpar;

	// //
	// cout<<"Kv="<<kn<<", "<<Upar<<" "<<Wpar<<endl;
	// //

	
	//force
	double f[3];
	for(int i = 0; i < Tpart; ++i){
	  int n, ni;
	  part.get_kind(i, n, ni);
	  for(int a = 0; a < 3; ++a)
	    ri[a] = part.get_pos(n, ni, a);
	  qi = part.get_charge(n);

	  if(qi != 0.){
	    double fact = 0.;
	    //as force is not scalar, different values of "a" are not always equivalent.
	    for(int j = 0; j < 3; ++j)
	      f[j] = 0.;
	    for(int b = 0; b < 4; ++b){
	      fact = -val * imag(qi * M[b][i] * conj(sum[b])); //self force term = 0
	      for(int j = 0; j < 3; ++j)
		f[j] += fact * K[b][j];
	    }
	    if((nx == 0 && ny == 0) || (nx == 0 && nz == 0) || (ny == 0 && nz == 0)){
	      for(int j = 0; j < 3; ++j)
		f[j] /= 4.;
	    }else if((nz == 0 && nx != 0 && ny != 0) || (ny == 0 && nz != 0 && nx != 0) || (nx == 0 && ny != 0 && nz != 0)){
	      for(int j = 0; j < 3; ++j)
		f[j] /= 2.;
	      
	    }
	    for(int a = 0; a < 3; ++a)
	      part.add_F(n, ni, a, f[a]);

	  }
	}
	
      }
    }
  }//end parallel
  U = Uc;
  W = Wc;

  // //
  // for(int i = 0; i < Tpart; ++i){
  //   int n, ni;
  //   part.get_kind(i, n, ni);
  //   cout<<part.get_F(n,ni,0)<<" "<<part.get_F(n,ni,1)<<" "<<part.get_F(n,ni,2)<<endl;
  // }
  // cout<<endl;
  // //
}


 
 
void moveX(Particles& part, double Dt){
  int Tpart = 0;
  for(int i = 0; i < part.get_Nkinds(); ++i)
    Tpart += part.get_N(i);

  //  #pragma omp_set_num_threads(OMP_THR);
  #pragma omp parallel default(shared)
  {
    int tasks2 = 8 * omp_get_num_threads();
    int chunk = static_cast<int>(floor(static_cast<float>(Tpart)/tasks2));
    #pragma omp for schedule(dynamic, chunk) nowait
    for(int i = 0; i < Tpart; ++i){
      int m, mi;
      part.get_kind(i, m, mi);
      
      double mass = part.get_mass(m);
      
      for(int a = 0; a < 3; a++){
	double x = part.get_pos(m, mi, a) + Dt * part.get_mom(m, mi, a) / mass; 
	part.set_pos(m, mi, a, x);
      }
    }
  }//end parallel
}

void moveP(Particles& part, const Thermostat& thermo, double Dt, double& K, bool last){
  int Tpart = 0;
  for(int i = 0; i < part.get_Nkinds(); ++i)
    Tpart += part.get_N(i);
  
  double Dt2 = Dt / 2;
  double Kc0 = 0.;

  //  #pragma omp_set_num_threads(OMP_THR);
  #pragma omp parallel default(shared)
  {
    int tasks2 = 8 * omp_get_num_threads();
    int chunk = static_cast<int>(floor(static_cast<float>(Tpart)/tasks2));
    #pragma omp for schedule(dynamic, chunk) reduction(+: Kc0) nowait
    for(int i = 0; i < Tpart; ++i){ 
      int m, mi;
      part.get_kind(i, m, mi);
      double p;
      double psq = 0.;
      double s = thermo.get_coo();
      double mass = part.get_mass(m);
      
      if(!last){
	double sfac = 1. + Dt2 * s;
	for(int a = 0; a < 3; a++){
	  p = (part.get_mom(m, mi, a) + Dt2 * part.get_F(m, mi, a)) / sfac; 
	  part.set_mom(m, mi, a, p);
	  psq += p * p;
	}
      }else{
	double sfac = 1. - Dt2 * s;
	for(int a = 0; a < 3; a++){
	  p = part.get_mom(m, mi, a) * sfac + Dt2 * part.get_F(m, mi, a);
	  part.set_mom(m, mi, a, p);
	  psq += p * p;
	}
      }
      Kc0 += psq / mass;
      
    }
  }//end parallel
  K = Kc0 / 2;
}

void moveThermo(Thermostat& thermo, double temp, double Dt, double K, double Tpart, int mode){
  double rnd0, rnd1, g, s, ex1;
  double Q = thermo.get_mass();
  
  switch(mode){
   case 0:
     rnd0 = thermo.get_rand();
     rnd1 = thermo.get_rand();
     gauss(1, rnd0, rnd1);
     thermo.set_strn(rnd1);
     g = thermo.get_shr();
     ex1 = exp(-g*Dt/2);
     s = ex1*thermo.get_coo() + sqrt(2.*temp*g*Dt/Q)*rnd0;
     thermo.set_coo(s);
     break;
   case 1:
     s = thermo.get_coo() + Dt*(2.*K - 3.*Tpart*temp)/Q;
     thermo.set_coo(s);
     break;
   case 2:
     rnd0 = thermo.get_strn();
     g = thermo.get_shr();
     ex1 = exp(-g*Dt/2);
     s = ex1*thermo.get_coo() + sqrt(2.*temp*g*Dt/Q)*rnd0;
     thermo.set_coo(s);
     break;
   default:
     break;
  }
}



bool save_cont(double Time, long cont, int Ntcf, const accumulator& AcT, const accumulator& AcTsq){
  char name[30];
  sprintf(name, "counters.dat");
  ofstream Con(name);
  if(! Con){
    Con.close();
    return false;
  }
  Con << setiosflags(ios::fixed) << setprecision(12);
  Con << Time << " " << cont << " " << Ntcf << " ";
  Con << AcT.get_K() << " " << AcTsq.get_K() << " "
      << AcT.get_U() << " " << AcTsq.get_U() << " "
      << AcT.get_P() << " " << AcTsq.get_P() << " "
      << AcT.get_Uss() << " " << AcTsq.get_Uss() << " "
      << AcT.get_Pss() << " " << AcTsq.get_Pss() << " "
      << AcT.get_UchR() << " " << AcTsq.get_UchR() << " "
      << AcT.get_PchR() << " " << AcTsq.get_PchR() << " "
      << AcT.get_UchK() << " " << AcTsq.get_UchK() << " "
      << AcT.get_PchK() << " " << AcTsq.get_PchK() << " ";
  Con << endl;
  Con.close();
  return true;
}

bool chrg_cont(double& Time, long& cont, int& Ntcf, accumulator& AcT, accumulator& AcTsq){
  double var;
  char name[30];
  sprintf(name, "counters.dat");
  ifstream Con(name);
  if(! Con){
    Con.close();
    return false;
  }
  Con >> Time >> cont >> Ntcf;
  Con >> var; AcT.set_K(var);
  Con >> var; AcTsq.set_K(var);
  Con >> var; AcT.set_U(var);
  Con >> var; AcTsq.set_U(var);
  Con >> var; AcT.set_P(var);
  Con >> var; AcTsq.set_P(var);
  Con >> var; AcT.set_Uss(var);
  Con >> var; AcTsq.set_Uss(var);
  Con >> var; AcT.set_Pss(var);
  Con >> var; AcTsq.set_Pss(var);
  Con >> var; AcT.set_UchR(var);
  Con >> var; AcTsq.set_UchR(var);
  Con >> var; AcT.set_PchR(var);
  Con >> var; AcTsq.set_PchR(var);
  Con >> var; AcT.set_UchK(var);
  Con >> var; AcTsq.set_UchK(var);
  Con >> var; AcT.set_PchK(var);
  Con >> var; AcTsq.set_PchK(var);
  Con.close();
  return true;
}



void acc_correl(const Particles& part, Hist& histograms, const double L[3], double bar_w){
  int NC = part.get_Nkinds();
  //int NS = (NC-1)*(NC-1) - ((NC-1)*(NC-2))/2;
  int Tpart = 0;
  for(int i = 0; i < NC; ++i)
    Tpart += part.get_N(i);
  
  //int cst;
  int n_bars = histograms.get_nbars();	


  #pragma omp parallel default(shared) //private(P)
  {
    int tasks2 = 16 * omp_get_num_threads();
    int chunk = static_cast<int>(static_cast<float>(Tpart)/tasks2);
    #pragma omp for schedule(dynamic, chunk) nowait
    for(int f = 0; f < Tpart; ++f){
      double ri[3], rj[3];
      // m and i
      int m=0, i=0;
      part.get_kind(f, m, i);
      for(int a = 0; a < 3; ++a)
	ri[a] = part.get_pos(m, i, a);
      //
      double radi = part.get_rad(m);
      //
            
      //even out the charge among threads
      int mx = static_cast<int>(ceil(static_cast<float>(Tpart-1)/2));
      if(fmod(static_cast<float>(Tpart),2) == 0. && f >= Tpart/2)
	mx = static_cast<int>(floor(static_cast<float>(Tpart-1)/2));
      
      int g = f+1 - Tpart*static_cast<int>(floor((f+1)/Tpart + 0.5));
      int cnt = 0;
      while(cnt < mx){
	// n and j
	int n=0, j=0;
	part.get_kind(g, n, j);
	for(int a = 0; a < 3; ++a)
	  rj[a] = part.get_pos(n, j, a);
	//
	double radj = part.get_rad(n);
	//

	double rij[3];
	for(int a = 0; a < 3; ++a){
	  rij[a] = ri[a] - rj[a];
	  rij[a] -= L[a] * floor(rij[a] / L[a] + 0.5);
	}

	double RIJSQ;
	dot(rij, rij, RIJSQ);
	double RIJ = sqrt(RIJSQ);    
	
	int bin = static_cast<int>(floor(RIJ / bar_w));
	if(bin < n_bars){
	  int mx_k = max(m,n);
	  int mn_k = min(m,n);
	  int pos_hist = 0;
	  for(int l = 0; l < mn_k; ++l)
	    pos_hist += NC - l;
	  pos_hist += mx_k - mn_k;
	  pos_hist *= n_bars;
	  pos_hist += bin;

	  // //
	  // if((bin<50) && (m == 0 || n == 0) && histograms.get_hist(pos_hist)!=0.){
	  //   printf("bin %d RIJ=%f <? %f !\n", bin, RIJ, (radi+radj)*MX_OVRLP);
	  // }
	  // //

	  double val = 1.;
	  if(n == m) val = 2.;
	  histograms.add(pos_hist, val); //000
	}
	g += 1 - Tpart*static_cast<int>(floor((g+1)/Tpart + 0.5));
	cnt++;
      }

    }
  }//end parallel
}


bool nor_correl(const Hist& histograms, const Particles& part, double Lmin, double bar_w, double cor_cont){
  int NC = part.get_Nkinds();
	  
  int n_bars = histograms.get_nbars();

  char name[30];
  sprintf(name, "correlation.dat");
  ofstream HIST(name);
  if(! HIST){
    HIST.close();
    return false;
  }
  
  double pd = 4. * M_PI;
  double bar3 = bar_w * bar_w * bar_w / 12.;
  double sq2 = sqrt(2.);
  double hnrm;
  
  for(int bin = 0; bin < n_bars; ++bin){
    double pos = (static_cast<double>(bin) + 0.5) * bar_w;
    double norm;
    if(pos >= Lmin / 2. && pos <= Lmin / sq2){
      norm = pd * pos * pos * ( 3 * Lmin / (2. * pos) - 2.) * bar_w;
    }// else if(pos > Lmin / sq2){
    //}
    else
      norm = pd * (pos * pos * bar_w + bar3);
    double nnorm = norm * cor_cont;

    HIST << cor_cont << " " 
	 << setiosflags(ios::fixed) << setprecision(8) << norm << " "
	 << setiosflags(ios::fixed) << setprecision(4) << pos << " & ";
    HIST << setiosflags(ios::fixed) << setprecision(8);
    
    for(int m = 0; m < NC; ++m){
      int MM = part.get_N(m);
      for(int n = m; n < NC; ++n){
	int NN = part.get_N(n);
	if(m == n)
	  NN -= 1;

	int pos_hist = 0;
	for(int l = 0; l < m; ++l)
	  pos_hist += NC - l;
	pos_hist += n - m;
	pos_hist *= n_bars;
	pos_hist += bin;

	// //
	// if((bin<50) && (m == 0 || n == 0) && histograms.get_hist(pos_hist)!=0.){
	//   cout<<"bin="<<bin<<" pos="<<pos<<" pos_hist="<<pos_hist<<endl;
	// }
	// //
	
	if(NN > 0){
	  hnrm = histograms.get_hist(pos_hist) / (nnorm * MM * NN);
	  HIST << " " << hnrm;
	}else
	  HIST << " " << 0;
	HIST << " |";
      }
      HIST << "|";
    }
    HIST << endl;
  }
  HIST.close();
  return true;
}

bool chrg_correl(double& cor_cont, Hist& histograms, const Particles& part, double bar_w){
  int NC = part.get_Nkinds();
 
  int n_bars = histograms.get_nbars();

  char name[30];
  sprintf(name, "correlation.dat");
  ifstream HIST(name);
  if(! HIST){
    HIST.close();
    return false;
  }
  string line;
  double pos, norm;
  double hnrm;
    
  for(int bin = 0; bin < n_bars; ++bin){
    HIST >> cor_cont >> norm >> pos >> line;
    double nnorm = norm * cor_cont;

    for(int m = 0; m < NC; ++m){
      int MM = part.get_N(m);
      for(int n = m; n < NC; ++n){
	int NN = part.get_N(n);
	if(m == n)
	  NN -= 1;
	
	int pos_hist = 0;
	for(int l = 0; l < m; ++l)
	  pos_hist += NC - l;
	pos_hist += n - m;
	pos_hist *= n_bars;
	pos_hist += bin;

	HIST >> hnrm;
	histograms.set_hist(pos_hist, hnrm*(nnorm * MM * NN)); 
	HIST >> line;
      }
    }
  }
  HIST.close();
  return true;
}


//for computing mean square disp. and velocity autocorrelation
bool pos_vel(const Particles& part, vector<int>& np_tc, double Dt, int num){
  int NC = np_tc.size();
  
  char name[130];
  getcwd(name,130);
  sprintf(name+strlen(name),"/pos_vel");
  mkdir(name, 0777);
  sprintf(name+strlen(name), "/pos_vel_%d.dat", num);
  std::ofstream Con(name);
  if(! Con){
    Con.close();
    return false;
  }
  Con << setiosflags(ios::fixed) << setprecision(12);
  double r[3], p[3];
  for(int n = 0; n < NC; ++n){
    int NN = np_tc[n];
    double m = part.get_mass(n);
    for(int i = 0; i < NN; ++i){
      for(int a = 0; a < 3; ++a){
	r[a] = part.get_pos(n, i, a);
	p[a] = part.get_mom(n, i, a);
      }
      Con << r[0] << " " << r[1] << " " << r[2] << " "
	  << p[0]/m << " " << p[1]/m << " " << p[2]/m << endl;
    }
  }
  Con << Dt << endl;
  Con.close();
  return true;
}

 
void dot(const double a[3], const double b[3], double& ans){
  ans = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

bool print_out(double Time, const accumulator& var, long telap, long tcpu, int iprint, double tpu){
  //time difference evaluation
  double t_elap = (telap / iprint) * 1.e-6; //time / step
  double t_cpu = (static_cast<double>(tcpu) / iprint) * 1.e-6; //cpu time / step
  double t_unit = (telap / tpu) * 1.e-6; //time / simulation time unit
  //print
  char name[30];
  sprintf(name, "output.dat");
  ofstream Out(name, ios::app);
  if(! Out){
    Out.close();
    return false;
  }
  Out << Time << " " << t_elap << " " << t_cpu << " " << t_unit <<" | ";
  Out << setiosflags(ios::fixed) << setprecision(12);
  Out << var.get_K() << " " << var.get_U() << " " << var.get_P() << " " 
      << var.get_Uss() << " " << var.get_Pss() << " " 
      << var.get_UchR() << " " << var.get_PchR() << " "
      << var.get_UchK() << " " << var.get_PchK() << endl;
  Out.close();
  return true;
}

bool print_final(double Time, double avK, double flK, double avU, double flU, double avP, double flP, double avUss, double flUss, double avPss, double flPss, double avUchR, double flUchR, double avPchR, double flPchR, double avUchK, double flUchK, double avPchK, double flPchK, const double L[3], const Particles& part){
  int NC = part.get_Nkinds();
  double *mass = part.get_M();
  double *rad = part.get_R();
  double *charge = part.get_Q();
  //print
  char name[30];
  sprintf(name, "final.dat");
  ofstream Out(name, ios::app);
  if(! Out){
    Out.close();
    return false;
  }
  Out << "NVT ensemble\n";
  Out << "Side of the box lenght (Lx,Ly,Lz): " 
      << L[0] << "," << L[1] << "," << L[2] << endl;
  Out << "Elapsed time: " << Time << endl;
  Out << "Component | No. of particles | mass | charge | rad | Vol. frac.\n"; 
  double V = L[0] * L[1] * L[2];
  for(int n = 0; n < NC; ++n){
    int NN =part.get_N(n);
    Out << n << " | " << NN << " | " << mass[n] << " | " << charge[n] << " | " << rad[n] << " | " << NN * 4. * M_PI * rad[n] * rad[n] * rad[n] / (3 * V) << endl;
  }
  Out << "Av. kinetic energy: " << avK << "+/-" << flK << endl;
  Out << "Av. potential energy: " << avU << "+/-" << flU << endl;
  Out << "Av. soft sphere en.: " << avUss << "+/-" << flUss 
      << " Av. Coul.R en.: " << avUchR << "+/-" << flUchR 
      << " Av. Coul.K en.: " << avUchK << "+/-" << flUchK << endl;
  Out << "Av. pressure: " << avP << "+/-" << flP << endl;
  Out << "Av. soft sphere pres.: " << avPss << "+/-" << flPss 
      << " Av. Coul.R pres.: " << avPchR << "+/-" << flPchR
      << " Av. Coul.K pres.: " << avPchK << "+/-" << flPchK << endl;
  Out.close();
  return true;
}








//mix initial configuration (Geometric Clustering Algorithm)
bool GCAmix(Particles& part, const double L[3], int id, const gsl_rng* ran){
  int NC = part.get_Nkinds();
  int Tpart = 0;
  for(int i = 0; i < part.get_Nkinds(); ++i)
    Tpart += part.get_N(i);
  
  //added particles
  bool added[Tpart];
  for(int i = 0; i < Tpart; ++i)
    added[i] = false;
  
  //choose pivot
  double pivot[3];
  double rnd = gsl_rng_uniform(ran);
  pivot[0] = (2 * rnd - 1) * L[0] / 2;
  rnd = gsl_rng_uniform(ran);
  pivot[1] = (2 * rnd - 1) * L[1] / 2;
  rnd = gsl_rng_uniform(ran);
  pivot[2] = (2 * rnd - 1) * L[2] / 2;
  
  // //choose seed particle at random
  // rnd = gsl_rng_uniform(ran);
  // int id = static_cast<int>(floor(rnd*Tpart));
  int m, mi;
  part.get_kind(id, m, mi);

  //add first particle
  added[id] = true;
  deque<int> queue;
  queue.push_back(id);

  
  bool test = false;

  //construct cluster
  double ri[3], rj[3];
  while(queue.size() > 0){
    id = queue[0];
    part.get_kind(id, m, mi);
    //reflect
    for(int a = 0; a < 3; ++a)
      ri[a] = part.get_pos(m, mi, a);
    double XIJ = pivot[0] - ri[0];
    XIJ -= L[0] * floor(XIJ / L[0] + 0.5);
    double YIJ = pivot[1] - ri[1];
    YIJ -= L[1] * floor(YIJ / L[1] + 0.5);
    double ZIJ = pivot[2] - ri[2];
    ZIJ -= L[2] * floor(ZIJ / L[2] + 0.5);
    ri[0] += 2*XIJ; ri[1] += 2*YIJ; ri[2] += 2*ZIJ; 
    ri[0] -= L[0] * floor(ri[0] / L[0] + 0.5);
    ri[1] -= L[1] * floor(ri[1] / L[1] + 0.5);
    ri[2] -= L[2] * floor(ri[2] / L[2] + 0.5);
    for(int a = 0; a < 3; ++a)
      part.set_pos(m, mi, a, ri[a]);
    
    //if particle is a test particle add the other one
    if(!test){
      if(id == 0){
	added[1] = true;
	queue.push_back(1);
	test = true;
      }else if(id == 1){
	added[0] = true;
	queue.push_back(0);
	test = true;
      }
    }

    for(int jd = 0; jd < Tpart; ++jd){
      if(!added[jd]){
	int n, nj;
	part.get_kind(jd, n, nj);
	double sigma = part.get_rad(m) + part.get_rad(n);
	double sig_sq = sigma * sigma;
	for(int a = 0; a < 3; ++a)
	  rj[a] = part.get_pos(n, nj, a);

	double XIJ = ri[0] - rj[0];
	XIJ -= L[0] * floor(XIJ / L[0] + 0.5);
	double YIJ = ri[1] - rj[1];
	YIJ -= L[1] * floor(YIJ / L[1] + 0.5);
	double ZIJ = ri[2] - rj[2];
	ZIJ -= L[2] * floor(ZIJ / L[2] + 0.5);
	double RIJSQ = XIJ * XIJ + YIJ * YIJ + ZIJ * ZIJ;
	
	if(RIJSQ < sig_sq){
	  queue.push_back(jd);
	  added[jd] = true;
	  
	  //if particle is a test particle add the other one
	  if(!test){
            if(id == 0){
              added[1] = true;
              queue.push_back(1);
	      test = true;
            }else if(jd == 1){
              added[0] = true;
	      queue.push_back(0);
	      test = true;
            }
          }

	}
      }
    }
    queue.pop_front();
  }

  // //test for overlaps
  // if(overlap_test(part, L, 1.)){
  //   cout<<"overlap at GCAmix!"<<endl;
  //   return false;
  // }
  return true;
}
