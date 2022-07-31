#include "functions.hpp"
#include "dev_functions.hpp"

#define DEV false

int main(){
  int omp_get_thread_num();
  # pragma omp parallel
  {
    printf("Thread rank: %d\n", omp_get_thread_num());
  }

  int Nsteps;
  double pre_time; //time lapse before averaging 
  int iprint; //prints/saves frequency in steps
  int signal;
  int Nkinds; //how many kinds of particles
  vector<int> N; //No. of particles per kind
  vector<double> m; //masses
  vector<double> r; //radii
  vector<double> q; //charges
  int K_max;  //number of K vectors per dimension
  double dens; //density
  double Dt; //time interval
  double bar_w; //histogram bar width
  double ArLx, ArLy; //aspect ratio of the box (Lx/Lz and Ly/Lz)
  double temp; //temperature
  double Q; //thermostat mass
  double gamma; //thermostat shear rate
  double cut; //real part cutoff multipier (1. is minimum image convention)
  double alpha; //Ewald convergencve parameter
  double corrf; //correlation averaging frequency
  double tcfd; // time correlation sampling frequency (diffusion)
  int mx_tcf; //maximum No. of part. per species to compute time correl.


  r_input(Nsteps, pre_time, iprint, signal, Nkinds, N, m, r, q, dens, Dt, ArLx, ArLy, temp, Q, gamma, cut, K_max, alpha, bar_w, corrf, tcfd, mx_tcf);

  double block_time = 1000*Dt;

  // // *******
  // cout << Nsteps << "\n";
  // cout << pre_time << "\n";
  // cout << iprint << "\n";
  // cout << signal << "\n";
  // cout << Nkinds << "\n";
  // for(int i = 0; i < Nkinds; ++i)
  //   cout << N[i] << " ";
  // cout << "\n";
  // for(int i = 0; i < Nkinds; ++i)
  //   cout << m[i] << " ";
  // cout << "\n";
  // for(int i = 0; i < Nkinds; ++i)
  //   cout << r[i] << " ";
  // cout << "\n";
  // for(int i = 0; i < Nkinds; ++i)
  //   cout << q[i] << " ";
  // cout << "\n";
  // cout << dens << "\n";
  // cout << Dt << "\n";
  // cout << ArLx << " " << ArLy << "\n"; 
  // cout << temp << "\n";
  // cout << Q << "\n";
  // cout << gamma << "\n";
  // cout << cut << "\n";
  // cout << K_max << "\n";
  // cout << alpha << "\n";
  // cout << bar_w << "\n";
  // cout << corrf << "\n";
  // cout << tcfd << "\n";
  // cout << mx_tcf << "\n";
  // // *******


  //number of part to take time correlations
  vector<int> np_tc(N.size());
  for(int i = 0; i < N.size(); ++i)
    N[i] <= mx_tcf ? np_tc[i] = N[i] : np_tc[i] = mx_tcf;

  //variable time step
  double dt = Dt;
  
  
  //initialize particles
  Particles part(Nkinds, N, m, r, q);
  int NC = part.get_Nkinds();
  //int NS = (NC-1)*(NC-1) - ((NC-1)*(NC-2))/2;
  
  int Tpart = 0;
  for(int i = 0; i < part.get_Nkinds(); ++i)
    Tpart += part.get_N(i);

  //error
  ofstream Err("Error.dat");
  Err.close();

  //check number of particles
  for(int n = 0; n < NC; ++n){
    if(part.get_N(n) < 1){
      Err.open("Error.dat", ios::app);
      Err << "No. of particles per species must be >= 1" << endl;
      Err.close();
      return 1;
    }
  }

  //box dimensions
  double L[3], Lmin, cut_d;
  
  //initialize global average counters
  long cont = 0;
  double Time = 0.;
  accumulator AcT;
  accumulator AcTsq;
  int Ntcf = 0; //number of time correlation files
  // int Ntcv = 0; //size of pres. t. correlation file
  

  //initialize random num. gen.
  gsl_rng* ran = gsl_rng_alloc(gsl_rng_mt19937); //mersenne twister rng
  long seed = 5.;//time(NULL);
  gsl_rng_set(ran, seed);

  //initialize thermostat
  Thermostat thermo; //for translation 
  thermo.set_mass(Q);
  thermo.set_shr(gamma);
  thermo.set_seed(333/*time(NULL)*/);


  //correlation histograms
  double cor_cont = 0;
  int n_bars = 0;
  Hist histograms;
  
  //counter to reset Dt (admits up to 10 time scales)
  int scal_cnt[10] = {0};
  
  //begin or continue run
  if(signal == 0){ //start new run
    L[2] = pow(static_cast<double>(Tpart) / (ArLx*ArLy*dens), 1./3);
    L[1] = ArLy * L[2];
    L[0] = ArLx * L[2];

    // ////
    // L[0]=110.; L[1]=110.; L[2]=110.;
    // ////

    //min L
    Lmin = min(L[2], min(L[1], L[0]));
    //real part cutoff distance
    cut_d = Lmin * cut / 2;
    
    //
    //This is done because if the initial box is too large some particles seem
    // to adquire a lot of momentum. WHY?? THIS NEEDS TO BE SOLVED!
    double LL[3]={min(50.,L[0]), min(50.,L[1]), min(50.,L[2])};
    //

    //initial position
    Err.open("Error.dat", ios::app);
    Err << "generating positions..." << endl;
    Err.close();
    ini_pos(part, LL/*L*/);
    Err.open("Error.dat", ios::app);
    Err << "positions generated." << endl;
    Err.close();
    for(int s = 0; s < Tpart; ++s){
      if(!GCAmix(part, LL/*L*/, s, ran)){
    	Err.open("Error.dat", ios::app);
    	Err<<"ovrlp. at "<<s<<" mix\n";
    	Err.close();
    	//break;
      }
    }
    //initial momenta
    Err.open("Error.dat", ios::app);
    Err << "generating velocities..." << endl;
    Err.close();
    ini_mom(part, temp, ran);
    Err.open("Error.dat", ios::app);
    Err << "velocities generated." << endl;
    Err.close();
    //print configuration
    if(! print_conf(part, thermo, L, dt, scal_cnt, true)){
      Err.open("Error.dat", ios::app);
      Err << "problem opening conf.dat" << endl;
      Err.close();
      gsl_rng_free(ran);
      return 1;
    }
    //allocate histograms
    n_bars = static_cast<int>(floor(L[2] / bar_w));
    histograms.set_size(NC, n_bars);
  }else if(signal == 1){  //charge config. Do not continue averaging.
    //charge configuration
    if(! chrg_conf(part, thermo, L, dt, scal_cnt)){
      Err.open("Error.dat", ios::app);
      Err << "problem opening conf.dat for reading\n";
      Err.close();
      gsl_rng_free(ran);
      return 1;
    }
    //min L
    Lmin = min(L[2], min(L[1], L[0]));
    //real part cutoff distance
    cut_d = Lmin * cut / 2;
    //charge
    for(int n = 0; n < NC; ++n)
      part.set_charge(n, q[n]);
    //allocate histograms
    n_bars = static_cast<int>(floor(L[2] / bar_w));
    histograms.set_size(NC, n_bars);
  }else if(signal == 2){
    //charge counters
    if(! chrg_cont(Time, cont, Ntcf, AcT, AcTsq)){
      Err.open("Error.dat", ios::app);
      Err << "problem opening counters.dat for reading\n";
      Err.close();
      gsl_rng_free(ran);
      return 1;
    }
    //charge configuration
    if(! chrg_conf(part, thermo, L, dt, scal_cnt)){
      Err.open("Error.dat", ios::app);
      Err << "problem opening conf.dat for reading\n";
      Err.close();
      gsl_rng_free(ran);
      return 1;
    }
    //min L
    Lmin = min(L[2], min(L[1], L[0]));   
    //real part cutoff distance
    cut_d = Lmin * cut / 2;
    //allocate histograms
    n_bars = static_cast<int>(floor(L[2] / bar_w));
    histograms.set_size(NC, n_bars);
    //charge correlations
    if(! chrg_correl(cor_cont, histograms, part, bar_w)){
      Err.open("Error.dat", ios::app);
      Err << "problem opening correlation.dat for reading\n";
      Err.close();
      gsl_rng_free(ran);
      return 1;
    }
  }else{
    Err.open("Error.dat", ios::app);
    Err << "Incorrect signal given.\n";
    Err.close();
    gsl_rng_free(ran);
    return 1;
  }

  
  //check total linear momentum 
  double ptot[3] = {0,0,0};
  for(int n = 0; n < NC; ++n){
    int NN = part.get_N(n);
    for(int i = 0; i < NN; ++i){
      for(int k = 0; k < 3; ++k)
  	ptot[k] += part.get_mom(n, i, k);
    }
  }
  if(ptot[0]+ptot[1]+ptot[2] > 1.e-6){
    Err.open("Error.dat", ios::app);
    Err << "total momentum is not 0! ";
    Err << ptot[0]<< "," << ptot[1] << "," << ptot[2] <<endl;
    Err.close();
    gsl_rng_free(ran);
    return 1;
  }
  
  //
  //sum squared momenta on each direction
  ofstream MOM("mom.dat", ios::app);
  double sqp[3] = {0,0,0};
  MOM.close();
  //


  //check charge neutrality
  double QQ = 0.;
  for(int n = 0; n < NC; ++n)
    QQ += part.get_N(n) * part.get_charge(n);
  if(fabs(QQ) > 1.e-12){
    Err.open("Error.dat", ios::app);
    Err << "total charge is not 0! " << QQ <<endl;
    Err.close();
    gsl_rng_free(ran);
    return 1;
  }

  double V = L[0] * L[1] * L[2];

  //partial counters
  long scont = 0;
  accumulator Ac;

  //TD quantities
  accumulator var;
  double K=0, U=0, P=0; //kinetic energy, potential energy and Pressure
  double Uss=0, Wss=0, Pss=0; //soft sphere contributions
  double UchR=0, WchR=0, PchR=0; //real sum charge contributions
  double UchK=0, WchK=0, PchK=0; //reciprocal sum charge contributions
  
  
  //tabulation for the real part of energy
  if(Lmin > 200){
    Err.open("Error.dat", ios::app);
    Err << "Tab is not large enough. Lmin=" << Lmin << " > 200\n";
    Err.close();
    gsl_rng_free(ran);
    return 1;
  }  
  double *Tab;
  Tab = new double[2*TAB_SIZE];
  for(int i = 0; i < 2*TAB_SIZE; ++i)
    Tab[i] = 0.;
  tabul(Tab, alpha);
  //tabulation for the reciprocal part of energy
  int K_size = (K_max+1)*(K_max+1)*(K_max+1);
  double *Kvec;
  Kvec = new double[K_size];
  wave_v(Kvec, L, alpha, K_max);

  //
  cout<<"start\n";
  //

  bool *reduce = new bool[3];
  reduce[0] = false; reduce[1] = false; reduce[2] = false; 

  if(!DEV){
    //
    cout<<"ini r\n";
    //
    //initial force R
    if(! force_R(part, L, Tab, cut_d, Uss, Wss, UchR, WchR, dt, reduce)){
      Err.open("Error.dat", ios::app);
      Err << "initial overlap!\n";
      Err.close();
      gsl_rng_free(ran);
      delete[] Tab;
      delete[] Kvec;
      return 1;
    }
    cout<<"ini k\n";
    //initial force K
    force_K(part, L, K_max, alpha, Kvec, UchK, WchK);
  }else{
    //Initial force R (with device)
    if(! dev_force_R(part, L, Tab, cut_d, Uss, Wss, UchR, WchR, dt, reduce)){
      Err.open("Error.dat", ios::app);
      Err << "initial overlap!\n";
      Err.close();
      gsl_rng_free(ran);
      delete[] Tab;
      delete[] Kvec;
      return 1;
    }
    //initial force K (with device)
    dev_force_K(part, L, K_max, alpha, Kvec, UchK, WchK);
  }

  delete[] reduce;
  
  // //
  // cout<<Uss<<" "<<UchR<<" "<<UchK<<endl;
  // //

  //initialize calculation time
  Timer stepT;
  double telap=0.; long tcpu=0;
  Timer force_chaRT;
  double force_Rel=0.; long force_Rpu=0;
  Timer force_chaKT;
  double force_Kel=0.; long force_Kpu=0;
  Timer movXT;
  double movX_el=0.; long movX_pu=0;
  Timer movPT;
  double movP_el=0.; long movP_pu=0;
  Timer corT;
  double cor_el=0.; long cor_pu=0;
  
  
  //time per iprint steps
  double tpu = 1.e-6;
  //to check weather or not to accumulate
  double prev_tc = 0.;
  double prev_tacc = 0.;
  double prev_tcf = 0.;
  
  //Begin step loop************************************************************
  for(int step = 0; step < Nsteps; ++step){
    // //
    if(fmod(static_cast<float>(step), iprint) == 0.){
      // cout<<step<<endl;
      sqp[0] = 0.; sqp[1] = 0.; sqp[2] = 0.;
      MOM.open("mom.dat", ios::app);
      for(int n = 0; n < NC; ++n){
	int NN = part.get_N(n);
	for(int i = 0; i < NN; ++i){
	  for(int k = 0; k < 3; ++k){
	    double p = part.get_mom(n, i, k);
	    sqp[k] += p*p;
	  }
	}
      }
      MOM<<Time<<" "<<sqp[0]<<" "<<sqp[1]<<" "<<sqp[2]<<endl;
      MOM.close();
    }
    // //

    //LEAPFROG VERLET--------------------------------------------------------
    double Ktr = 0.;
    //thermostat first step
    moveThermo(thermo, temp, dt, Ktr, Tpart, 0);
    //find intermediate velocities
    if(!DEV)
      moveP(part, thermo, dt, Ktr, false);
    else
      dev_moveP(part, thermo, dt, Ktr, false);
    //thermostat second step
    moveThermo(thermo, temp, dt, Ktr, Tpart, 1);

    //find new positions
    movXT.start();
    if(!DEV)
      moveX(part, dt);
    else
      dev_moveX(part, dt);
    movXT.stop(telap, tcpu);
    movX_el += telap; movX_pu += tcpu;
    

    reduce = new bool[3];
    reduce[0] = false; reduce[1] = false; reduce[2] = false; 
    
    //calculate new forces R
    force_chaRT.start();
    if(!DEV){
      if(! force_R(part, L, Tab, cut_d, Uss, Wss, UchR, WchR, dt, reduce)){
  	Err.open("Error.dat", ios::app);
  	Err << "overlap at step "<<step<<"!\n";
  	Err.close();
  	gsl_rng_free(ran);
  	delete[] Tab;
  	delete[] Kvec;
	delete[] reduce;
  	return 1;
      }
    }else{
      if(! dev_force_R(part, L, Tab, cut_d, Uss, Wss, UchR, WchR, dt, reduce)){
    	Err.open("Error.dat", ios::app);
    	Err << "overlap at step "<<step<<"!\n";
    	Err.close();
    	gsl_rng_free(ran);
    	delete[] Tab;
    	delete[] Kvec;
	delete[] reduce;
    	return 1;
      }
    }
    force_chaRT.stop(telap, tcpu);
    force_Rel += telap; force_Rpu += tcpu;

    //calculate new forces K
    force_chaKT.start();
    if(!DEV)
      force_K(part, L, K_max, alpha, Kvec, UchK, WchK);
    else
      dev_force_K(part, L, K_max, alpha, Kvec, UchK, WchK);
    force_chaKT.stop(telap, tcpu);
    force_Kel += telap; force_Kpu += tcpu;

    //find new velocity
    movPT.start();
    if(!DEV)
      moveP(part, thermo, dt, Ktr, true);
    else
      dev_moveP(part, thermo, dt, Ktr, true);
    movPT.stop(telap, tcpu);
    movP_el += telap; movP_pu += tcpu;
    
    //thermostat third step
    moveThermo(thermo, temp, dt, Ktr, Tpart, 2);
    //-----------------------------------------------------------------------
    
    //check total linear momentum conservation
    for(int a = 0; a < 3; ++a)
      ptot[a] = 0.;
    for(int n = 0; n < NC; ++n){
      int NN = part.get_N(n);
      for(int i = 0; i < NN; ++i){
    	for(int k = 0; k < 3; ++k)
    	  ptot[k] += part.get_mom(n, i, k);
      }
    }
    if(ptot[0]+ptot[1]+ptot[2] > 1.e-6){
      Err.open("Error.dat", ios::app);
      Err << "total momentum is not 0! ";
      Err << ptot[0]<< "," << ptot[1] << "," << ptot[2] <<endl;
      Err.close();
      gsl_rng_free(ran);
      delete[] Tab;
      delete[] Kvec;
      delete[] reduce;
      return 1;
    }
    
    

    //time step resizing--------------
    // //
    // cout<<reduce[0]<<" "<<reduce[1]<<" "<<reduce[2]<<endl;
    // //
    //check for rescaling level
    int tc = 0;
    while(scal_cnt[tc] > 0)
      tc++;
    if(tc < 10){
      if(reduce[0]){ //if reduction of time step is needed...
  	dt /= TIME_MULT;
  	scal_cnt[tc] += TSTEP_MULT * 10;
/*}else*/if(reduce[1] && tc < 9){ //if further reduction is needed...
	  dt /= TIME_MULT;//100;
	  //scal_cnt[tc] += TSTEP_MULT * 10;
  	  scal_cnt[tc+1] += TSTEP_MULT * 10;
	}
      }else if(tc > 0 && reduce[2] && scal_cnt[tc-1]==1){ //if a pair of particles needs more time
  	scal_cnt[tc-1] += TSTEP_MULT * 10; //mantain the reduction
      }else if(tc > 0 && fabs(dt-Dt) > 1.e-10){ //if dt different from Dt...
  	scal_cnt[tc-1]--;
  	if(scal_cnt[tc-1] == 0) //increase the time step
  	  dt *= TIME_MULT;
      }
    }
    
    // //
    // //if(scal_cnt[0] > 0){
    // cout<<"step "<<step<<" ( ";
    //   for(int tc = 0; tc < 10; ++tc)
    // 	cout<<scal_cnt[tc]<<" ";
    //   cout<<") dt= "<<dt<<endl;
    //   cout<<endl;
    //   //}
    // //

    delete[] reduce;

    
    K = Ktr / Tpart;
    Uss /= Tpart; UchR /= Tpart; UchK /= Tpart;
    U = Uss + UchR + UchK;
    Pss = Wss / (3 * V); PchR = WchR / (3 * V); PchK = WchK / (3 * V);
    double ctemp = 2 * K / 3;
    P = ctemp * dens + Pss + PchR + PchK;

    var.set_K(K); var.set_U(U); var.set_P(P);
    var.set_Uss(Uss); var.set_Pss(Pss); 
    var.set_UchR(UchR); var.set_PchR(PchR);
    var.set_UchK(UchK); var.set_PchK(PchK); 

    Time += dt;

    //accumulate data
    if(Time >= pre_time){
      
      if(prev_tc > corrf){
  	//accumulate histograms
  	corT.start();
  	if(!DEV)
  	  acc_correl(part, histograms, L, bar_w);
  	else
  	  dev_acc_correl(part, histograms, L, bar_w);
  	cor_cont += 1. / V;
  	corT.stop(telap, tcpu);
  	cor_el += telap;
  	cor_pu += tcpu;
  	prev_tc = 0.;
      }
      prev_tc += dt;

      //print files for time correlation functions
      if(prev_tcf > tcfd){
  	if(! pos_vel(part, np_tc, tcfd, Ntcf)){
  	  Err.open("Error.dat", ios::app);
  	  Err << "problem opening pos_vel_"<< Ntcf <<".dat." << endl;
  	  Err.close();
  	}
  	prev_tcf = 0.;
  	Ntcf++;
      }
      prev_tcf += dt;
      
      //accumulate counters
      Ac += var;
      scont++;

      //acumulate global counters
      if(prev_tacc > block_time){
      	accumulator Acv = Ac / scont;
      	AcT += Acv;
      	AcTsq += (Acv * Acv);
      	cont++;
  	prev_tacc = 0.;
	Ac.set_zero();
	scont = 0;
      }
      prev_tacc += dt;
    }
    
    //save data
    if(fmod(static_cast<float>(step), iprint) == 0.){
      //write configuration//
      if(! print_conf(part, thermo, L, dt, scal_cnt, true)){
  	Err.open("Error.dat", ios::app);
  	Err << "problem opening conf.dat\n"; 
  	Err.close();
  	delete[] Tab;
  	delete[] Kvec;
  	gsl_rng_free(ran);
  	return 1;
      }
      //print counters//
      if(! save_cont(Time, cont, Ntcf, AcT, AcTsq)){
  	Err.open("Error.dat", ios::app);
  	Err << "problem opening counters.dat\n";
  	Err.close();
      }
      //normalize and print correlations//
      if(Time >= pre_time){
  	if(! nor_correl(histograms, part, L[2], bar_w, cor_cont)){
  	  Err.open("Error.dat", ios::app);
  	  Err << "problem opening correlation.dat" << endl;
  	  Err.close();
	  delete[] Tab;
	  delete[] Kvec;
	  gsl_rng_free(ran);
  	  return 1;
  	}
      }
      //print output//
      stepT.stop(telap, tcpu);
      // //
      // cout<<telap<<" "<<tpu<<endl;
      // //
      if(! print_out(Time, var, telap, tcpu, iprint, tpu)){
  	Err.open("Error.dat", ios::app);
  	Err << "problem opening output.dat\n"; 
  	Err.close();
  	delete[] Tab;
  	delete[] Kvec;
  	gsl_rng_free(ran);
  	return 1;
      }
      stepT.start();
      tpu = 0;
    }


    tpu += dt;
  }
  //End step loop*************************************************************

  //free Tables
  delete[] Tab;
  delete[] Kvec;
  //free rng
  gsl_rng_free(ran);

  //
  ofstream Tfile("time.dat");
  Tfile << "Force_chaR time. elap: " << (force_Rel/(Nsteps))*1e-6 
  	<< " cpu: " << (static_cast<double>(force_Rpu)/(Nsteps))*1e-6 << endl;
  Tfile << "Force_chaK time. elap: " << (force_Kel/(Nsteps))*1e-6 
  	<< " cpu: " << (static_cast<double>(force_Kpu)/(Nsteps))*1e-6 << endl;
  Tfile << "X move time. elap: " << (movX_el/(Nsteps))*1e-6 
  	<< " cpu: " << (static_cast<double>(movX_pu)/(Nsteps))*1e-6  << endl;
  Tfile << "P move time. elap: " << (movP_el/(Nsteps))*1e-6 
  	<< " cpu: " << (static_cast<double>(movP_pu)/(Nsteps))*1e-6  << endl;
  Tfile << "correl. time. elap: " << (cor_el/(V*cor_cont))*1e-6 
  	<< " cpu: " << (static_cast<double>(cor_pu)/(V*cor_cont))*1e-6 << endl;
  Tfile.close();
  //

  //calculate total averages and fluctuations

  double avK = AcT.get_K() / cont;
  double flK = sqrt(AcTsq.get_K() / cont - avK * avK);
  double avU = AcT.get_U() / cont;
  double flU = sqrt(AcTsq.get_U() / cont - avU * avU);
  double avP = AcT.get_P() / cont;
  double flP = sqrt(AcTsq.get_P() / cont - avP * avP);
  double avUss = AcT.get_Uss() / cont;
  double flUss = sqrt(AcTsq.get_Uss() / cont - avUss * avUss);
  double avPss = AcT.get_Pss() / cont;
  double flPss = sqrt(AcTsq.get_Pss() / cont - avPss * avPss);
  double avUchR = AcT.get_UchR() / cont;
  double flUchR = sqrt(AcTsq.get_UchR() / cont - avUchR * avUchR);
  double avPchR = AcT.get_PchR() / cont;
  double flPchR = sqrt(AcTsq.get_PchR() / cont - avPchR * avPchR);
  double avUchK = AcT.get_UchK() / cont;
  double flUchK = sqrt(AcTsq.get_UchK() / cont - avUchK * avUchK);
  double avPchK = AcT.get_PchK() / cont;
  double flPchK = sqrt(AcTsq.get_PchK() / cont - avPchK * avPchK);
  
  print_final(Time, avK, flK, avU, flU, avP, flP, avUss, flUss, avPss, flPss, avUchR, flUchR, avPchR, flPchR, avUchK, flUchK, avPchK, flPchK, L, part);
  
  
  return 0;
}
