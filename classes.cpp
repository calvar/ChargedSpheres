#include "classes.hpp"

//Particles******************************************************
Particles::Particles(): Nkinds(1) {
  Npart = new int[1];
  M = new double[1];
  R = new double[1];
  Q = new double[1];

  pos = new double[3];
  mom = new double[3];
  force = new double[3];
}

Particles::Particles(int kinds, const vector<int>& N, const vector<double>& m, const vector<double>& r, const vector<double>& q): Nkinds(kinds) {
  Npart = new int[Nkinds];
  M = new double[Nkinds];
  R = new double[Nkinds];  
  Q = new double[Nkinds];
  
  int tot = 0;
  for(int i = 0; i < Nkinds; ++i){
    Npart[i] = N[i];
    M[i] = m[i];  
    R[i] = r[i];
    Q[i] = q[i];
    tot += N[i];
  }

  pos = new double[3*tot];
  mom = new double[3*tot];
  force = new double[3*tot];
}

Particles::~Particles() {
  delete[] Npart;
  delete[] M;
  delete[] R;
  delete[] Q;
  
  delete[] pos;
  delete[] mom;
  delete[] force;
}

int Particles::get_Nkinds() {
  return Nkinds;
}

int Particles::get_N(int i) {
  return Npart[i];
}

int* Particles::get_Npart() {
  return Npart;
}

double Particles::get_mass(int i) {
  return M[i];
}

double* Particles::get_M() {
  return M;
}

double Particles::get_rad(int i) {
  return R[i];
}

double* Particles::get_R() {
  return R;
}

double Particles::get_charge(int i) {
  return Q[i];
}

double* Particles::get_Q() {
  return Q;
}

int Particles::get_Nkinds() const {
  return Nkinds;
}

int Particles::get_N(int i) const {
  return Npart[i];
}

int* Particles::get_Npart() const {
  return Npart;
}

double Particles::get_mass(int i) const {
  return M[i];
}

double* Particles::get_M() const{
  return M;
}

double Particles::get_rad(int i) const {
  return R[i];
}

double* Particles::get_R() const {
  return R;
}

double Particles::get_charge(int i) const {
  return Q[i];
}

double* Particles::get_Q() const {
  return Q;
}

void Particles::set_charge(int i, double q) {
  Q[i] = q;
}


//Get coordinate "a" from particle "i" of species "n"; 
double Particles::get_pos(int n, int i, int a) {
  //check boundaries
  if(i > Npart[n]-1) cout << "Particle index out of bounds.\n";
  if(a > 2) cout << "Component out of bounds.\n";

  //get the index of the particle
  int part = 0;
  for(int l = 0; l < n; ++l)
    part += Npart[l];
  part += i;

  //get component "a" of the position
  return pos[part*3+a];
}

//Get coordinate "a" from particle "i" of species "n"; 
double Particles::get_pos(int n, int i, int a) const {
  //check boundaries
  if(i > Npart[n]-1) cout << "Particle index out of bounds.\n";
  if(a > 2) cout << "Component out of bounds.\n";

  //get the index of the particle
  int part = 0;
  for(int l = 0; l < n; ++l)
    part += Npart[l];
  part += i;

  //get component "a" of the position
  return pos[part*3+a];
}

//Set coordinate "a" from particle "i" of species "n"
void Particles::set_pos(int n, int i, int a, double val){
  //check boundaries
  if(i > Npart[n]-1) cout << "Particle index out of bounds.\n";
  if(a > 2) cout << "Component out of bounds.\n";
  
  //get the index of the particle
  int part = 0;
  for(int l = 0; l < n; ++l)
    part += Npart[l];
  part += i;

  pos[part*3+a] = val;
}

double* Particles::get_X(){
  return pos;
}

double* Particles::get_X() const {
  return pos;
}

// //Get coordinate "a" from particle "id"
// double Particles::get_pos_id(int id, int a) {
//   //check boundaries
//   if(i > Nkinds-1) cout << "Particle index out of bounds.\n";
//   if(a > 2) cout << "Component out of bounds.\n";

//   //get component "a" of the position
//   return pos[id*3+a];
// }

// //Get coordinate "a" from particle "id"
// double Particles::get_pos_id(int id, int a) const{
//   //check boundaries
//   if(i > Nkinds-1) cout << "Particle index out of bounds.\n";
//   if(a > 2) cout << "Component out of bounds.\n";

//   //get component "a" of the position
//   return pos[id*3+a];
// }

// //Set coordinate "a" from particle "id"
// void Particles::set_pos_id(int id, int a, double val){
//   //check boundaries
//   if(i > Nkinds-1) cout << "Particle index out of bounds.\n";
//   if(a > 2) cout << "Component out of bounds.\n";

//   pos[id*3+a] = val;
// }

//Get momentum component "a" from particle "i" of species "n"; 
double Particles::get_mom(int n, int i, int a) {
  //check boundaries
  if(i > Npart[n]-1) cout << "Particle index out of bounds.\n";
  if(a > 2) cout << "Component out of bounds.\n";

  //get the index of the particle
  int part = 0;
  for(int l = 0; l < n; ++l)
    part += Npart[l];
  part += i;

  //get component "a" of the position
  return mom[part*3+a];
}

//Get momentum component "a" from particle "i" of species "n"; 
double Particles::get_mom(int n, int i, int a) const {
  //check boundaries
  if(i > Npart[n]-1) cout << "Particle index out of bounds.\n";
  if(a > 2) cout << "Component out of bounds.\n";

  //get the index of the particle
  int part = 0;
  for(int l = 0; l < n; ++l)
    part += Npart[l];
  part += i;

  //get component "a" of the position
  return mom[part*3+a];
}

//Set momentum component "a" from particle "i" of species "n"
void Particles::set_mom(int n, int i, int a, double val){
  //check boundaries
  if(i > Npart[n]-1) cout << "Particle index out of bounds.\n";
  if(a > 2) cout << "Component out of bounds.\n";
  
  //get the index of the particle
  int part = 0;
  for(int l = 0; l < n; ++l)
    part += Npart[l];
  part += i;

  mom[part*3+a] = val;
}

double* Particles::get_P(){
  return mom;
}

double* Particles::get_P() const {
  return mom;
}

//Get force component "a" on particle "i" of species "n"; 
double Particles::get_F(int n, int i, int a) {
  //check boundaries
  if(i > Npart[n]-1) cout << "Particle index out of bounds.\n";
  if(a > 2) cout << "Component out of bounds.\n";

  //get the index of the particle
  int part = 0;
  for(int l = 0; l < n; ++l)
    part += Npart[l];
  part += i;

  //get component "a" of the position
  return force[part*3+a];
}

//Get force component "a" on particle "i" of species "n"; 
double Particles::get_F(int n, int i, int a) const{
  //check boundaries
  if(i > Npart[n]-1) cout << "Particle index out of bounds.\n";
  if(a > 2) cout << "Component out of bounds.\n";

  //get the index of the particle
  int part = 0;
  for(int l = 0; l < n; ++l)
    part += Npart[l];
  part += i;

  //get component "a" of the position
  return force[part*3+a];
}

//Get force vector 
double* Particles::get_F() {
  return force;
}

//Get force vector 
double* Particles::get_F() const{
  return force;
}

//Set forces to zero
void Particles::set_zero_F(){
  int tot = 0;
  for(int i = 0; i < Nkinds; ++i)
    tot += Npart[i];
  
  for(int i = 0; i < 3*tot; ++i)
    force[i] = 0.;
}

//Add to force component "a" on particle "i" of species "n"
void Particles::add_F(int n, int i, int a, double val){
  //check boundaries
  if(i > Npart[n]-1) cout << "Particle index out of bounds.\n";
  if(a > 2) cout << "Component out of bounds.\n";
  
  //get the index of the particle
  int part = 0;
  for(int l = 0; l < n; ++l)
    part += Npart[l];
  part += i;

  #pragma omp atomic
  force[part*3+a] += val;
}

//Set the force vector
void Particles::set_F(double *Fin, int size){
  int Tpart = 0;
  for(int n =0; n < Nkinds; ++n)
    Tpart += Npart[n];
  if(size != 3*Tpart) cout << "Wrong array size in set_F!\n";

  for(int i = 0; i < size; ++i)
    force[i] = Fin[i];
}

//Get the kind and number of the particle
void Particles::get_kind(int id, int& n, int& i){
  int Nto=0, Ntn=0;
  for(int nn = 0; nn < Nkinds; ++nn){
    Nto = Ntn;
    Ntn += Npart[nn];
    if(id < Ntn){
      n = nn;
      i = id - Nto;
      break;
    }
  }
}

//Get the kind and number of the particle
void Particles::get_kind(int id, int& n, int& i) const{
  int Nto=0, Ntn=0;
  for(int nn = 0; nn < Nkinds; ++nn){
    Nto = Ntn;
    Ntn += Npart[nn];
    if(id < Ntn){
      n = nn;
      i = id - Nto;
      break;
    }
  }
}


//Thermostat********************************************************
Thermostat::Thermostat() : m_Q(1.), m_s(0.) , m_g(0.) { 
  m_ran = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(m_ran, 1);
}

Thermostat::Thermostat(double Q, double g, long seed) : m_Q(Q), m_s(0.), m_g(g) { 
  m_ran = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(m_ran, seed);
}

Thermostat::Thermostat(const Thermostat& T) {
  m_Q = T.m_Q;
  m_s = T.m_s;
  m_g = T.m_g;
  m_ran = T.m_ran;
  m_strn = T.m_strn;
}

Thermostat::~Thermostat() { gsl_rng_free(m_ran); }

void Thermostat::operator=(const Thermostat& T) {
  m_Q = T.m_Q;
  m_s = T.m_s;
  m_g = T.m_g;
  m_ran = T.m_ran;
  m_strn = T.m_strn;
}

double Thermostat::get_mass() { return m_Q; }

double Thermostat::get_mass() const{ return m_Q; }

void Thermostat::set_mass(double Q) { m_Q = Q; }

double Thermostat::get_coo() { return m_s; }

double Thermostat::get_coo() const{ return m_s; }

void Thermostat::set_coo(const double s) { m_s = s; }

double Thermostat::get_shr() { return m_g; }

double Thermostat::get_shr() const{ return m_g; }

void Thermostat::set_shr(const double g) { m_g = g; }

void Thermostat::set_seed(const long seed) { gsl_rng_set(m_ran, seed); }

double Thermostat::get_rand() { return gsl_rng_uniform(m_ran); }

void Thermostat::set_strn(const double r) { m_strn = r; }

double Thermostat::get_strn() { return m_strn; }

double Thermostat::get_strn() const{ return m_strn; }

bool Thermostat::print_ran(unsigned n, bool fin) {
  char name[30];
  char c;
  if(! fin) c = 'I';
  else c = 'F';
  sprintf(name, "ran_state%d%c.dat", n, c);
  FILE* Out = fopen(name, "w");
  if(Out == NULL)
    return false;
  gsl_rng_fwrite(Out, m_ran);
  fclose(Out);
  return true;
}

bool Thermostat::print_ran(unsigned n, bool fin) const{
  char name[30];
  char c;
  if(! fin) c = 'I';
  else c = 'F';
  sprintf(name, "ran_state%d%c.dat", n, c);
  FILE* Out = fopen(name, "w");
  if(Out == NULL)
    return false;
  gsl_rng_fwrite(Out, m_ran);
  fclose(Out);
  return true;
}

bool Thermostat::read_ran(unsigned n) {
  char name[30];
  sprintf(name, "ran_state%dF.dat", n);
  FILE* In = fopen(name, "r");
  if(In == NULL)
    return false;
  gsl_rng_fread(In, m_ran);
  fclose(In);
  return true;
}


//acumulator******************************************************
accumulator::accumulator() : m_K(0.), m_U(0.), m_P(0.), m_Uss(0.), m_Pss(0.),
		  m_UchR(0.), m_PchR(0.), m_UchK(0.), m_PchK(0.) {}

void accumulator::set_K(double K) { m_K = K; }

double accumulator::get_K() { return m_K; }

double accumulator::get_K() const{ return m_K; }

void accumulator::set_U(double U) { m_U = U; }

double accumulator::get_U() { return m_U; }

double accumulator::get_U() const{ return m_U; }

void accumulator::set_P(double P) { m_P = P; }

double accumulator::get_P() { return m_P; }

double accumulator::get_P() const{ return m_P; }

void accumulator::set_Uss(double Uss) { m_Uss = Uss; }

double accumulator::get_Uss() { return m_Uss; }

double accumulator::get_Uss() const{ return m_Uss; }

void accumulator::set_Pss(double Pss) { m_Pss = Pss; }

double accumulator::get_Pss() { return m_Pss; }

double accumulator::get_Pss() const{ return m_Pss; }

void accumulator::set_UchR(double UchR) { m_UchR = UchR; }

double accumulator::get_UchR() { return m_UchR; }

double accumulator::get_UchR() const{ return m_UchR; }

void accumulator::set_PchR(double PchR) { m_PchR = PchR; }

double accumulator::get_PchR() { return m_PchR; }
 
double accumulator::get_PchR() const{ return m_PchR; }

void accumulator::set_UchK(double UchK) { m_UchK = UchK; }

double accumulator::get_UchK() { return m_UchK; }

double accumulator::get_UchK() const{ return m_UchK; }

void accumulator::set_PchK(double PchK) { m_PchK = PchK; }

double accumulator::get_PchK() { return m_PchK; }

double accumulator::get_PchK() const{ return m_PchK; }

void accumulator::set_zero() {
  m_K = 0.; m_U = 0.; m_P = 0.;
  m_Uss = 0.; m_Pss = 0.; m_UchR = 0.; m_PchR = 0.; 
  m_UchK = 0.; m_PchK = 0.;
}

void accumulator::operator=(const accumulator& acc) {
  m_K = acc.m_K;
  m_U = acc.m_U;
  m_P = acc.m_P;
  m_Uss = acc.m_Uss;
  m_Pss = acc.m_Pss;
  m_UchR = acc.m_UchR;
  m_PchR = acc.m_PchR;
  m_UchK = acc.m_UchK;
  m_PchK = acc.m_PchK;
}

void accumulator::operator+=(const accumulator& acc) {
  m_K += acc.m_K;
  m_U += acc.m_U;
  m_P += acc.m_P;
  m_Uss += acc.m_Uss;
  m_Pss += acc.m_Pss;
  m_UchR += acc.m_UchR;
  m_PchR += acc.m_PchR;
  m_UchK += acc.m_UchK;
  m_PchK += acc.m_PchK;
}
  
const accumulator accumulator::operator*(const accumulator& acc) {
  accumulator result;
  result.m_K = m_K * acc.m_K;
  result.m_U = m_U * acc.m_U;
  result.m_P = m_P * acc.m_P;
  result.m_Uss = m_Uss * acc.m_Uss;
  result.m_Pss = m_Pss * acc.m_Pss;
  result.m_UchR = m_UchR * acc.m_UchR;
  result.m_PchR = m_PchR * acc.m_PchR;
  result.m_UchK = m_UchK * acc.m_UchK;
  result.m_PchK = m_PchK * acc.m_PchK;
  return result;
}

const accumulator accumulator::operator/(long a) {
  accumulator result;
  result.m_K = m_K / a;
  result.m_U = m_U / a;
  result.m_P = m_P / a;
  result.m_Uss = m_Uss / a;
  result.m_Pss = m_Pss / a;
  result.m_UchR = m_UchR / a;
  result.m_PchR = m_PchR / a;
  result.m_UchK = m_UchK / a;
  result.m_PchK = m_PchK / a;
  return result;
}

//Histograms*****************************************************
Hist::Hist() {
  m_NC = 0;
  m_nbars = 0;
  m_hist = new double[1];
}

Hist::~Hist() {
  delete[] m_hist;
}

void Hist::set_size(int NC, int n_bars) {
  m_NC = NC;
  m_nbars = n_bars;
  int totn = 0;
  for(int i = 0; i < m_NC; ++i)
    totn += m_NC - i;
  m_size = totn * m_nbars;
  delete[] m_hist;
  m_hist = new double[m_size];
  for(int i = 0; i < m_size; ++i)
    m_hist[i] = 0.;
}

int Hist::get_size() {
  return m_size;
}

int Hist::get_size() const{
  return m_size;
}

int Hist::get_NC() {
  return m_NC;
}

int Hist::get_NC() const{
  return m_NC;
}

int Hist::get_nbars() {
  return m_nbars;
}

int Hist::get_nbars() const{
  return m_nbars;
}

double Hist::get_hist(int i) {
  return m_hist[i];
}

double Hist::get_hist(int i) const{
  return m_hist[i];
}

double* Hist::get_hist() {
  return m_hist;
}

double* Hist::get_hist() const{
  return m_hist;
}

void Hist::set_hist(int i, double val) {
  if(i < m_size)
    m_hist[i] = val;
}

void Hist::set_hist(double *hist, int size) {
  if(size != m_size) cout << "Wrong array size in set_hist!\n";

  for(int i = 0; i < m_size; ++i)
    m_hist[i] = hist[i];
}

void Hist::add(int i, double val) {
  if(i < m_size){ 
    #pragma omp atomic
    m_hist[i] += val;
  }
}


//Timer**********************************************************
Timer::Timer() {
  gettimeofday(&m_tim, NULL); 
  m_tcpu = clock();
}

void Timer::start() {
  gettimeofday(&m_tim, NULL); 
  m_tcpu = clock();
}

void Timer::stop(double& telap, long& twork) { 
  double t1 = m_tim.tv_sec * 1e+6 + m_tim.tv_usec;
  long start = m_tcpu;
  gettimeofday(&m_tim, NULL);
  double t2 = m_tim.tv_sec * 1e+6 + m_tim.tv_usec;
  m_tcpu = clock();
  long end = m_tcpu;
  telap = t2 - t1;
  twork = end - start;
}



// //Tab********************************************************
// Tab::Tab(double alpha) {
//   size = 32768;
//   extL = 200 * 1.001;
//   table = new double[size*2];
//   double alsq = alpha * alpha;
//   double alft = alsq * alsq;
//   double spi = sqrt(M_PI);
//   for(int i = 1; i < size; i++){
//     double r =  0.5 * extL * i / size;
//     double rsq = r * r;
//     double sr1 = 1.0 / r;
//     double ERFC = erfc(alpha * r);
//     double aer = 2.0 * alpha * exp(- alsq * rsq) / spi;
//     table[(i-1)] = ERFC * sr1;
//     table[(i-1)+1] = (ERFC * sr1) + aer;
//   }
// }

// Tab::~Tab() {
//   delete[] table;
// }

// double Tab::get(int n, int i) {
//   //check boundaries
//   if(i > size-1) cout << "Index out of bounds in Tab.\n";
//   if(i > 1) cout << "Component out of bounds in Tab.\n";
  
//   return table[n*2+i];
// }
