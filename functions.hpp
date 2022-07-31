#ifndef __FUNCTIONS_HPP
#define __FUNCTIONS_HPP


#include "classes.hpp"

bool r_input(int& Nsteps, double& pre_time, int& iprint, int& signal, int& Nkinds, vector<int>& N, vector<double>& m, vector<double>& r, vector<double>& q, double& dens, double& Dt, double& ArLx, double& ArLy, double& temp, double& Q, double& gamma, double& sr_cut, int& K_max, double& alpha, double& bar_w, double& corrf, double& tcfd, int& mx_tcf);

//Test for overlaps
bool overlap_test(const Particles& part, const double L[3], double frac);

int ini_pos(Particles& part, const double L[3]);

void ini_mom(Particles& part, double temp, const gsl_rng* ran);

void gauss(double variance, double& rnd1, double& rnd2);

bool print_conf(const Particles& part, const Thermostat& thermo, const double L[3], double dt, const int scal_cnt[10], bool fin);

bool chrg_conf(Particles& part, Thermostat& thermo, double L[3], double& dt, int scal_cnt[10]);

void tabul(double *Tab, double alpha);

bool force_R(Particles& part, const double L[3], const double *Tab, double cut_d, double& Usr, double& Wsr, double& Uc, double& Wc, double Dt, bool *reduce);

void wave_v(double *Kvec, const double L[3], double alpha, int K_max);

void force_K(Particles& part, const double L[3], int K_max, double alpha, const double *Kvec, double& U, double& W);

void moveX(Particles& part, double Dt);

void moveP(Particles& part, const Thermostat& thermo, double Dt, double& K, bool last);

void moveThermo(Thermostat& thermo, double temp, double Dt, double K, double Tpart, int mode);

bool save_cont(double Time, long cont, int Ntcf, const accumulator& AcT, const accumulator& AcTsq);

bool chrg_cont(double& Time, long& cont, int& Ntcf, accumulator& AcT, accumulator& AcTsq);

void acc_correl(const Particles& part, Hist& histograms, const double L[3], double bar_w);

bool nor_correl(const Hist& histograms, const Particles& part, double Lmin, double bar_w, double cor_cont);

bool chrg_correl(double& cor_cont, Hist& histograms, const Particles& part, double bar_w);

bool pos_vel(const Particles& part, vector<int>& np_tc, double Dt, int num);

void dot(const double a[3], const double b[3], double& ans);

bool print_out(double Time, const accumulator& var, long telap, long tcpu, int iprint, double tpu);

bool print_final(double Time, double avK, double flK, double avU, double flU, double avP, double flP, double avUss, double flUss, double avPss, double flPss, double avUchR, double flUchR, double avPchR, double flPchR, double avUchK, double flUchK, double avPchK, double flPchK, const double L[3], const Particles& part);

bool GCAmix(Particles& part, const double L[3], int id, const gsl_rng* ran);

#endif
