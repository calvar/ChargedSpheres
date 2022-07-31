#ifndef DEV_FUNC_HPP__
#define DEV_FUNC_HPP__

bool dev_force_R(Particles& part, const double L[3], const double *Tab, double cut_d, double& tot_Usr, double& tot_Wsr, double& tot_Uch, double& tot_Wch, double Dt, bool *reduce);

void dev_force_K(Particles& part, const double L[3], int K_max, double alpha, const double *Kvec, double& tot_Uch, double& tot_Wch);

void dev_moveX(Particles& part, double Dt);

void dev_moveP(Particles& part, const Thermostat& thermo, double Dt, double& K, bool last);

void dev_acc_correl(const Particles& part, Hist& histograms, const double L[3], double bar_w);

#endif
