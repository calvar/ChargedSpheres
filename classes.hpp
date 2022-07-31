#ifndef __CLASSES_HPP
#define __CLASSES_HPP

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <omp.h>
#include <sys/time.h>
#include <gsl/gsl_rng.h>
using namespace std;

#define MX_OVRLP 0.65
#define TSTEP_MULT 2
#define TIME_MULT 3
#define TAB_SIZE 32768 //2 to the power 15

//Class that stores particle positions and properties in 
// 1D dynamicaly alocated arays, so that copying them
// to the device is easy.************************************************
class Particles {
private:
  int Nkinds; //How many kinds of particles are there
  int* Npart; //How many particles of each kind are there
  double* M; //Mass of each kind of particle
  double* R; //Radius of each kind of particle
  double* Q; //Charge of each kind of particle

  double* pos; //Position of each particle (linearized 3*N array)
  double* mom; //Momentum of each particle (linearized 3*N array)
  double* force; //force on each particle  (linearized 3*N array)
  
public:
  Particles();
  Particles(int, const vector<int>&, const vector<double>&, const vector<double>&, const vector<double>&);
  ~Particles();

  int get_Nkinds();
  int get_N(int);
  int* get_Npart();
  double get_mass(int);
  double* get_M();
  double get_rad(int);
  double* get_R();
  double get_charge(int);
  double* get_Q();
  int get_Nkinds() const;
  int get_N(int) const;
  int* get_Npart() const;
  double get_mass(int) const;
  double* get_M() const;
  double get_rad(int) const;
  double* get_R() const;
  double get_charge(int) const;
  double* get_Q() const;
  void set_charge(int, double); 

  double get_pos(int, int, int);
  double get_pos(int, int, int) const;
  void set_pos(int, int, int, double);
  double* get_X();
  double* get_X() const;
  // double get_pos_id(int, int);
  // double get_pos_id(int, int) const;
  // void set_pos_id(int, int, double);
  double get_mom(int, int, int);
  double get_mom(int, int, int) const;
  void set_mom(int, int, int, double);
  double* get_P();
  double* get_P() const;
  double get_F(int, int, int);
  double get_F(int, int, int) const;
  double* get_F();
  double* get_F() const;
  void set_zero_F();
  void add_F(int, int, int, double);
  void set_F(double*, int);
  void get_kind(int, int&, int&);
  void get_kind(int, int&, int&) const;
};


//Class to store the thermostat properties********************************
class Thermostat {
private:
  double m_Q; //"mass" (mass * area)
  double m_s; //generalized coord. (1. / time)
  double m_g; //"shear rate" on the thermostat (1. / time)
  gsl_rng* m_ran; //random number generator
  double m_strn; //stored normal random number 
  
public:
  Thermostat();
  Thermostat(double Q, double g, long seed);
  Thermostat(const Thermostat& T);
  ~Thermostat();
  
  void operator=(const Thermostat&);
  double get_mass();
  double get_mass() const;
  void set_mass(double Q);
  double get_coo();
  double get_coo() const;
  void set_coo(const double s);
  double get_shr();
  double get_shr() const;
  void set_shr(const double g);
  void set_seed(const long seed);
  double get_rand();
  void set_strn(const double r);
  double get_strn();
  double get_strn() const;
  bool print_ran(unsigned n, bool fin);
  bool print_ran(unsigned n, bool fin) const;
  bool read_ran(unsigned n);
};

//Class to accumulate all the thermodynamic variables*********************
class accumulator {
private:
  double m_K, m_U, m_P;
  double m_Uss, m_Pss;
  double m_UchR, m_PchR;
  double m_UchK, m_PchK;
  
public:
  accumulator();
  ~accumulator() {}
  
  void set_K(double K);
  double get_K();
  double get_K() const;
  void set_U(double U);
  double get_U();
  double get_U() const;
  void set_P(double P);
  double get_P();
  double get_P() const;
  void set_Uss(double Uss);
  double get_Uss();
  double get_Uss() const;
  void set_Pss(double Pss);
  double get_Pss();
  double get_Pss() const;
  void set_UchR(double UchR);
  double get_UchR();
  double get_UchR() const;
  void set_PchR(double PchR);
  double get_PchR();
  double get_PchR() const;
  void set_UchK(double UchK);
  double get_UchK();
  double get_UchK() const;
  void set_PchK(double PchK);
  double get_PchK();
  double get_PchK() const;
  void set_zero();
  void operator=(const accumulator&);
  void operator+=(const accumulator&);
  const accumulator operator*(const accumulator&);
  const accumulator operator/(long);
};


//Histograms**********************************************************
class Hist {
 private:
  int m_NC, m_nbars, m_size;
  double *m_hist;

 public:
  Hist();
  ~Hist();

  void set_size(int NC, int n_bars);
  int get_size();
  int get_size() const;
  int get_NC();
  int get_NC() const;
  int get_nbars();
  int get_nbars() const;
  double get_hist(int i);
  double get_hist(int i) const;
  double* get_hist();
  double* get_hist() const;
  void set_hist(int i, double val);
  void set_hist(double* hist, int size);
  void add(int i, double val);
};


//Timer***************************************************************
class Timer {
 private:
  timeval m_tim;
  long m_tcpu;
  
 public:
  Timer();
  ~Timer() {}
  
   void start();
   void stop(double&, long&);
};


// //Tabulation of predetermined functions to speed up calculations*********
// class Tab {
// private:
//   int size; //size of the table (No. of bins)
//   double extL; //maximum distance defined for the table
//   double* table;

// public:
//   Tab(double alpha);
//   ~Tab();

//   double get(int n, int i);
// };


#endif
