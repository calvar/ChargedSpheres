
#include <cuda_runtime.h>
#include <cuComplex.h>
#include <cuda_profiler_api.h>
#include "classes.hpp"
#include "dev_functions.hpp"
using namespace std;

#define BLOCK_WIDTH 32
#define K_SIM_SIZE 4

//Implementation from the CUDA page: https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#atomic-functions
__device__ double atomicAdd(double* address, double val){
  unsigned long long int* address_as_ull =
    (unsigned long long int*)address;
  unsigned long long int old = *address_as_ull, assumed;
  
  do {
    assumed = old;
    old = atomicCAS(address_as_ull, assumed,
		    __double_as_longlong(val +
					 __longlong_as_double(assumed)));
    
    // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
  } while (assumed != old);
  
  return __longlong_as_double(old);
}

//Atomic add for cuDoubleComplex
__device__ void atAddComplex(cuDoubleComplex* a, cuDoubleComplex b){
  //transform the addresses of real and imag. parts to double pointers
  double *x = (double*)a;
  double *y = x+1; //add 1*sizeof(double) bits

  //use atomicAdd for double variables
  atomicAdd(x, cuCreal(b));
  atomicAdd(y, cuCimag(b));
}
  

//REAL PART--------------------------------------------------------------------

__global__ void force_kernel(const double *pos, const double *mom, const double *rad,
			     const double *chrg, const double *mass, const int *Npart,
			     int NC, const double* L, double cut_d, double Lh,
			     const double *Tab, double srt, double *U_sr,
			     double *W_sr, double *U_ch, double *W_ch,
			     double *F, double Dt, bool *reduce){
  int n = threadIdx.x + blockDim.x * blockIdx.x;
  // //initialize force
  // for(int a = 0; a < 3; ++a)
  //   F[3*n+a] = 0.;
  //shared variables for accumulated values
  __shared__ double sum_usr[BLOCK_WIDTH];
  __shared__ double sum_wsr[BLOCK_WIDTH];
  __shared__ double sum_uch[BLOCK_WIDTH];
  __shared__ double sum_wch[BLOCK_WIDTH];
  __syncthreads();
  
  double Usr = 0.;
  double Wsr = 0.;
  double Uch = 0.;
  double Wch = 0.;
  double rn, ri;
  double qn, qi;
  double mn, mi;
  
  int Tpart = 0;
  for(int l = 0; l < NC; ++l)
    Tpart += Npart[l];
  
  //to check for ovelaps
  double warn = 0.;
  
  if(n < Tpart){
    int k1 = 0;
    int Tt = 0;
    for(int l = 0; l < NC; ++l){
      Tt += Npart[l];
      if(n < Tt)
	break;
      k1++;
    }
    rn = rad[k1];
    qn = chrg[k1];
    mn = mass[k1];
    
    //even out the charge among threads
    int mx = static_cast<int>(ceil(static_cast<float>(Tpart-1)/2));
    if(fmod(static_cast<float>(Tpart),static_cast<float>(2.)) == 0. && n >= Tpart/2)
      mx = static_cast<int>(floor(static_cast<float>(Tpart-1)/2));
    
    int i = n+1 - Tpart*(static_cast<int>(floor((n+1)/Tpart + 0.5)));
    int cnt = 0;
    while(cnt < mx){
      int k2 = 0;
      Tt = 0;
      for(int l = 0; l < NC; ++l){
	Tt += Npart[l];
	if(i < Tt)
	  break;
	k2++;
      }
      ri = rad[k2];
      qi = chrg[k2];
      mi = mass[k2];
      
      double rij[3];
      for(int a = 0; a < 3; ++a){
	// //
	// printf("%f-%f\n", pos[n*3+a], pos[i*3+a]);
	// //
	rij[a] = pos[n*3+a] - pos[i*3+a];
	rij[a] -= L[a] * floor(rij[a] / L[a] + 0.5); 
      }
      
      //cutoff distance (cube)
      if(fabs(rij[0])<cut_d && fabs(rij[1])<cut_d && fabs(rij[2])<cut_d){
	double RIJsq = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
	double RIJ = sqrt(RIJsq);
	double sigma = rn + ri;
	double sig_sq = sigma * sigma;
	double rsh = srt * sigma;
	
	// ////
	//printf("(%d,%d) [%f] %f ",n,i,sigma,RIJ);
	// printf("(%d,%d) [%f,%f] %f ",n,i,qn,qi,RIJ);
	// ////


	//Check for posible future overlaps----
	//relative velocity
	double vij[3], vr=0.;
	for(int a = 0; a < 3; ++a){
	  vij[a] = mom[n*3+a]/mn - mom[i*3+a]/mi;
	  vr += vij[a] * rij[a]/RIJ;
	}
	if(vr < 0.){ //if particles are approaching...
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
	//----------------------------------
	
	
	//Check for overlaps
	warn = 0.;
	if(RIJ < MX_OVRLP * sigma)
	  warn = 1.e12;
	
	//Lennard-Jones shifted potential, virial and force---------------
	double SR2s = sig_sq / RIJsq;
	double SR6s = SR2s * SR2s * SR2s;
	double pot = 4 * SR6s * (SR6s - 1.);
	double UIJ = 0.;
	double WIJ = 0.;
	if(RIJ < rsh){
	  UIJ = pot + 1.;//shifted potential (epsilon = 1)
	  WIJ = 6 * (pot + 4*SR6s*SR6s);
	}
	// ////
	// printf("(%d,%d) %f,%f,%f ",n,i,UIJ,RIJ,rsh);
	// ////
	Usr += UIJ + warn;
	Wsr += WIJ;
	
	//short range force
	double f[3] = {0., 0., 0.};
	double FIJsr = WIJ / RIJsq;
	for(int a = 0; a < 3; ++a){
	  f[a] = FIJsr * rij[a];
	  atomicAdd(&F[n*3+a], f[a]);
	}
	for(int a = 0; a < 3; ++a){
	  atomicAdd(&F[i*3+a], -f[a]);
	}
	// //
	// printf("(%d,%d) %f,%f,%f F %f,%f,%f | %f,%f,%f\n",n,i,f[0],f[1],f[2],F[n*3+0],F[n*3+1],F[n*3+2],F[i*3+0],F[i*3+1],F[i*3+2]);
	// //
	//----------------------------------------------------------------
	
	double qq = qn * qi;
	//Coulomb potential, virial and force---------------------------------
	if(fabs(qq) > 1.e-12){
	  int ri = static_cast<int>(floor(RIJ / Lh));
	  double R_L0 = RIJ - Lh * ri;
	  double R_L1 = RIJ - Lh * (ri+1);
	  double A = (Tab[ri] * R_L0 - Tab[ri-1] * R_L1) / Lh;
	  double B = (Tab[ri+TAB_SIZE] * R_L0 - Tab[ri-1+TAB_SIZE] * R_L1) / Lh;
	  UIJ = qq * A;
	  WIJ = qq * B;
	  Uch += UIJ;
	  Wch += WIJ;
	  // //
	  // printf("%d,%d %f %f %f %f\n",n,i,Tab[ri],Tab[ri-1],Tab[2*ri],Tab[2*(ri-1)]);
	  // //

	  //Coulomb force
	  for(int a = 0; a < 3; ++a){
	    f[a] = WIJ * rij[a] / RIJsq;
	    atomicAdd(&F[n*3+a], f[a]);
	  }
	  for(int a = 0; a < 3; ++a){
	    atomicAdd(&F[i*3+a], -f[a]);
	  }
	}
	//-----------------------------------------------------------------
	// //
	// printf("(%d,%d) %f,%f,%f F %f,%f,%f | %f,%f,%f\n",n,i,f[0],f[1],f[2],F[n*3+0],F[n*3+1],F[n*3+2],F[i*3+0],F[i*3+1],F[i*3+2]);
	// //
      }
      i += 1 - Tpart*static_cast<int>(floor((i+1)/Tpart + 0.5));
      cnt++;
    }
  }

  if(reduce[0]){
    if(reduce[2]) //maintain reduction?
      reduce[2] = false;
    //if(reduce[1])
    //  reduce[0] = false;
  }
  
  sum_usr[threadIdx.x] = Usr;
  sum_wsr[threadIdx.x] = Wsr;
  sum_uch[threadIdx.x] = Uch;
  sum_wch[threadIdx.x] = Wch;
  __syncthreads();
  if(threadIdx.x == 0){ 
    U_sr[blockIdx.x] = 0.;
    for(int i = 0; i < BLOCK_WIDTH; ++i) U_sr[blockIdx.x] += sum_usr[i];
    W_sr[blockIdx.x] = 0.;
    for(int i = 0; i < BLOCK_WIDTH; ++i) W_sr[blockIdx.x] += sum_wsr[i];
    U_ch[blockIdx.x] = 0.;
    for(int i = 0; i < BLOCK_WIDTH; ++i) U_ch[blockIdx.x] += sum_uch[i];
    W_ch[blockIdx.x] = 0.;
    for(int i = 0; i < BLOCK_WIDTH; ++i) W_ch[blockIdx.x] += sum_wch[i];
    // //
    // printf("  %f\n", U_sr[blockIdx.x]);
    // //
  }
}


bool dev_force_R(Particles& part, const double L[3], const double *Tab, double cut_d, double& tot_Usr, double& tot_Wsr, double& tot_Uch, double& tot_Wch, double Dt, bool *reduce){
  //cudaProfilerStart(); //only needed when using the cuda profiler
  
  int NC = part.get_Nkinds();
  int Tpart = 0;
  for(int i = 0; i < part.get_Nkinds(); ++i)
    Tpart += part.get_N(i);
  
  double srt = pow(2., 1./6);
  
  double *d_pos, *d_chrg, *d_rad, *d_mass, *d_L;
  int *d_Nc;
  double *d_Usr, *d_Wsr, *d_F;
  double *d_Uch, *d_Wch;
  double *d_Tab;
  double *d_mom;
  bool *d_red;

  double* X = part.get_X();
  double* P = part.get_P();
  double* Q = part.get_Q();
  double* R = part.get_R();
  double* M = part.get_M();
  int* Nprt = part.get_Npart();

  double Lh = 0.5 * (200*1.001) / TAB_SIZE;
  
  //Set device
  cudaSetDevice(0);

  //Allocate and copy memory to device
  cudaMalloc((void**)&d_pos, sizeof(double)*(3*Tpart));
  cudaMemcpy(d_pos, X, sizeof(double)*(3*Tpart), cudaMemcpyHostToDevice);
  cudaMalloc((void**)&d_mom, sizeof(double)*(3*Tpart));
  cudaMemcpy(d_mom, P, sizeof(double)*(3*Tpart), cudaMemcpyHostToDevice);
  cudaMalloc((void**)&d_chrg, sizeof(double)*NC);
  cudaMemcpy(d_chrg, Q, sizeof(double)*NC, cudaMemcpyHostToDevice);
  cudaMalloc((void**)&d_rad, sizeof(double)*NC);
  cudaMemcpy(d_rad, R, sizeof(double)*NC, cudaMemcpyHostToDevice);
  cudaMalloc((void**)&d_mass, sizeof(double)*NC);
  cudaMemcpy(d_mass, M, sizeof(double)*NC, cudaMemcpyHostToDevice);
  cudaMalloc((void**)&d_Nc, sizeof(double)*NC);
  cudaMemcpy(d_Nc, Nprt, sizeof(double)*NC, cudaMemcpyHostToDevice);
  cudaMalloc((void**)&d_L, sizeof(double)*3);
  cudaMemcpy(d_L, L, sizeof(double)*3, cudaMemcpyHostToDevice);
  cudaMalloc((void**)&d_Tab, sizeof(double)*(2*TAB_SIZE));
  cudaMemcpy(d_Tab, Tab, sizeof(double)*(2*TAB_SIZE), cudaMemcpyHostToDevice);
  cudaMalloc((void**)&d_red, sizeof(bool)*3);
  cudaMemcpy(d_red, reduce, sizeof(bool)*3, cudaMemcpyHostToDevice);
  

  //Define grid and block size
  int blocks = Tpart / BLOCK_WIDTH;
  if(Tpart % BLOCK_WIDTH) blocks++;
  
  //Allocate memory for outputx
  cudaMalloc((void**)&d_Usr, sizeof(double)*blocks);
  cudaMalloc((void**)&d_Wsr, sizeof(double)*blocks);
  cudaMalloc((void**)&d_Uch, sizeof(double)*blocks);
  cudaMalloc((void**)&d_Wch, sizeof(double)*blocks);
  cudaMalloc((void**)&d_F, sizeof(double)*(3*Tpart));
  cudaMemset(d_F, 0, sizeof(double)*(3*Tpart));
  
  // Kernel invocation code
  dim3 dimGrid(blocks, 1, 1);
  dim3 dimBlock(BLOCK_WIDTH, 1, 1);
  force_kernel<<<dimGrid,dimBlock>>>(d_pos, d_mom, d_rad, d_chrg, d_mass, d_Nc, NC, d_L, cut_d, Lh, d_Tab, srt, d_Usr, d_Wsr, d_Uch, d_Wch, d_F, Dt, d_red);

  //copy output to host
  double *Usr, *Wsr, *Uch, *Wch, *F;
  Usr = new double[blocks];
  Wsr = new double[blocks];
  Uch = new double[blocks];
  Wch = new double[blocks];
  F = new double[3*Tpart];
  cudaMemcpy(Usr, d_Usr, sizeof(double)*blocks, cudaMemcpyDeviceToHost);
  cudaMemcpy(Wsr, d_Wsr, sizeof(double)*blocks, cudaMemcpyDeviceToHost);
  cudaMemcpy(Uch, d_Uch, sizeof(double)*blocks, cudaMemcpyDeviceToHost);
  cudaMemcpy(Wch, d_Wch, sizeof(double)*blocks, cudaMemcpyDeviceToHost);
  cudaMemcpy(F, d_F, sizeof(double)*(3*Tpart), cudaMemcpyDeviceToHost);
  part.set_F(F, 3*Tpart);
  cudaMemcpy(reduce, d_red, sizeof(bool)*3, cudaMemcpyDeviceToHost);
  if(reduce[0]){
    if(reduce[2]) //maintain reduction?
      reduce[2] = false;
    if(reduce[1])
      reduce[0] = false;
  }

  // //
  // cout<<reduce[0]<<" "<<reduce[1]<<" "<<reduce[2]<<endl;
  // // 
  
  tot_Usr = 0.;
  for(int i = 0; i < blocks; ++i)
    tot_Usr += Usr[i];
  tot_Wsr = 0.;
  for(int i = 0; i < blocks; ++i)
    tot_Wsr += Wsr[i];
  tot_Uch = 0.;
  for(int i = 0; i < blocks; ++i)
    tot_Uch += Uch[i];
  tot_Wch = 0.;
  for(int i = 0; i < blocks; ++i)
    tot_Wch += Wch[i];
  
  cudaFree(d_pos);
  cudaFree(d_rad);
  cudaFree(d_chrg);
  cudaFree(d_mass);
  cudaFree(d_L);
  cudaFree(d_Nc);
  cudaFree(d_Tab);
  cudaFree(d_mom);
  cudaFree(d_red);
  cudaFree(d_Usr);
  cudaFree(d_Wsr);
  cudaFree(d_Uch);
  cudaFree(d_Wch);
  cudaFree(d_F);
  delete[] Usr;
  delete[] Wsr;
  delete[] Uch;
  delete[] Wch;
  delete[] F;

  //cudaProfilerStop(); //only needed when using the cuda profiler
  
  //overlap
  if(tot_Usr > 1.e11){
    return false;
  }
  
  //cout << "gpu: "<<tot_Usr <<" "<<tot_Wsr<<" "<<tot_Uch<<" "<<tot_Wch<<endl;
  // for(int n = 0; n < NC; ++n){
  //   for(int i = 0; i < part.get_N(n); ++i){
  //     for(int a = 0; a < 3; ++a)
  // 	cout << part.get_F(n, i, a) << " ";
  //     cout << endl;
  //   }
  // }
  //cout<<endl;

  return true;
}
//-----------------------------------------------------------------------------




//RECIPROCAL PART -------------------------------------------------------------


__global__ void kvec_add(double *pos, double *chrg, int *Npart, int NC, double *K, cuDoubleComplex *sum, cuDoubleComplex *M, double *self){
  int i = threadIdx.x + blockDim.x * blockIdx.x;

  //shared variables for accumulated values
  __shared__ double sum_self[BLOCK_WIDTH];
  extern __shared__ cuDoubleComplex sum_sum[]; //size=2*K_SIM_SIZE*BLOCK_WIDTH
  __syncthreads();

  int Tpart = 0;
  for(int l = 0; l < NC; ++l)
    Tpart += Npart[l];

  double self_par = 0.;
  cuDoubleComplex sum_par[2*K_SIM_SIZE] = {0.};
  
  if(i < Tpart){
    //Particle charge
    int k1 = 0;
    int Tt = 0;
    for(int l = 0; l < NC; ++l){
      Tt += Npart[l];
      if(i < Tt)
        break;
      k1++;
    }
    double qi = chrg[k1];
    cuDoubleComplex qic = make_cuDoubleComplex(qi, 0.);

    //printf("p %d q %f\n",i, qi);

    if(qi != 0){
      //Particle position
      double ri[3];
      for(int a = 0; a < 3; ++a)
        ri[a] = pos[i*3+a];

      double RIK[K_SIM_SIZE];
      cuDoubleComplex z[K_SIM_SIZE];
      for(int b = 0; b < K_SIM_SIZE; ++b){
        RIK[b] = 0.;
        for(int a = 0; a < 3; ++a)
  	  RIK[b] += ri[a] * K[b*3 + a];
        z[b] = make_cuDoubleComplex(cos(RIK[b]), sin(RIK[b]));
      }

      for(int b = 0; b < K_SIM_SIZE; ++b){
        M[i*K_SIM_SIZE + b] = cuConj(z[b]);
        //atAddComplex(&sum[b], cuCmul(qic, M[i*K_SIM_SIZE + b]));
	sum_par[b] = cuCadd(sum_par[b], cuCmul(qic, M[i*K_SIM_SIZE + b]));
        cuDoubleComplex qiRIK = make_cuDoubleComplex(qi*RIK[b], 0.);
        //atAddComplex(&sum[b+K_SIM_SIZE], cuCmul(qiRIK, z[b]));
	sum_par[b+K_SIM_SIZE] = cuCadd(sum_par[b+K_SIM_SIZE], cuCmul(qiRIK, z[b]));
      }
      //atomicAdd(&self[0], qi*qi);
      self_par += qi*qi;
  
      // //
      // for(int b = 0; b < K_SIM_SIZE; ++b)
      //   printf("[%d (%f,%f)] ",b,cuCreal(z[b]),cuCimag(z[b]));
      // printf("\n");
      // //
    }
  }

  sum_self[threadIdx.x] = self_par;
  for(int b = 0; b < 2*K_SIM_SIZE; ++b)
    sum_sum[2*K_SIM_SIZE*threadIdx.x + b] = make_cuDoubleComplex(cuCreal(sum_par[b]), cuCimag(sum_par[b]));

  __syncthreads();
  if(threadIdx.x == 0){
    for(int i = 0; i < BLOCK_WIDTH; ++i) self[blockIdx.x] += sum_self[i];
    for(int i = 0; i < BLOCK_WIDTH; ++i){
      for(int b = 0; b < 2*K_SIM_SIZE; ++b)
	sum[2*K_SIM_SIZE*blockIdx.x + b] = cuCadd(sum[2*K_SIM_SIZE*blockIdx.x + b],sum_sum[2*K_SIM_SIZE*i + b]);
    }
  }
}

__global__ void k_force_kernel(double *chrg, int *Npart, int NC, cuDoubleComplex *sum, cuDoubleComplex *M, double *K, double *F, double val, int nx, int ny, int nz){
  int i = threadIdx.x + blockDim.x * blockIdx.x;

  int Tpart = 0;
  for(int l = 0; l < NC; ++l)
    Tpart += Npart[l];

  if(i < Tpart){
    //Particle charge
    int k1 = 0;
    int Tt = 0;
    for(int l = 0; l < NC; ++l){
      Tt += Npart[l];
      if(i < Tt)
        break;
      k1++;
    }
    cuDoubleComplex qic = make_cuDoubleComplex(chrg[k1], 0.);

    if(cuCreal(qic) != 0.){
      double fact = 0.;
      double f[3] = {0., 0., 0.};
      for(int b = 0; b < 4; ++b){
        cuDoubleComplex qiM = cuCmul(qic, M[i*K_SIM_SIZE + b]);
        fact = -val * cuCimag(cuCmul(qiM, cuConj(sum[b]))); //self force term = 0
        for(int a = 0; a < 3; ++a)
	  f[a] += fact * K[b*3 + a];
      }
      if((nx == 0 && ny == 0) || (nx == 0 && nz == 0) || (ny == 0 && nz == 0)){
        for(int a = 0; a < 3; ++a)
	  f[a] /= 4.;
      }else if((nz == 0 && nx != 0 && ny != 0) || (ny == 0 && nz != 0 && nx != 0) || (nx == 0 && ny != 0 && nz != 0)){
        for(int a = 0; a < 3; ++a)
	  f[a] /= 2.;
      }

      for(int a = 0; a < 3; ++a)
	F[i*3+a] += f[a];
        //atomicAdd(&F[i*3+a], f[a]);
    }
  }
}



void dev_force_K(Particles& part, const double L[3], int K_max, double alpha, const double *Kvec, double& tot_Uch, double& tot_Wch){
  int NC = part.get_Nkinds();
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
  tot_Uch = 0.;
  tot_Wch = 0.;

  //Set device
  cudaSetDevice(0);

  double* X = part.get_X();
  double* Q = part.get_Q();
  double* F = part.get_F();
  int* Nprt = part.get_Npart();

  //Allocate and copy memory to device
  double *d_pos, *d_chrg, *d_K;
  int *d_Nc;
  cudaMalloc((void**)&d_pos, sizeof(double)*(3*Tpart));
  cudaMemcpy(d_pos, X, sizeof(double)*(3*Tpart), cudaMemcpyHostToDevice);
  cudaMalloc((void**)&d_chrg, sizeof(double)*NC);
  cudaMemcpy(d_chrg, Q, sizeof(double)*NC, cudaMemcpyHostToDevice);
  cudaMalloc((void**)&d_Nc, sizeof(double)*NC);
  cudaMemcpy(d_Nc, Nprt, sizeof(double)*NC, cudaMemcpyHostToDevice);   
  cudaMalloc((void**)&d_K, sizeof(double)*(K_SIM_SIZE*3));
  
  //Define grid and block size
  int blocks = Tpart / BLOCK_WIDTH;
  if(Tpart % BLOCK_WIDTH) blocks++;
  
  //Allocate memory for accumulators
  cuDoubleComplex *d_M, *d_sum;
  double *d_self;
  cudaMalloc((void**)&d_M, sizeof(cuDoubleComplex)*(Tpart*K_SIM_SIZE));
  cudaMalloc((void**)&d_sum, sizeof(cuDoubleComplex)*(2*K_SIM_SIZE*blocks));
  cudaMalloc((void**)&d_self, sizeof(double)*blocks);
  

  cuDoubleComplex *d_sumc;
  double *d_force;
  cudaMalloc((void**)&d_sumc, sizeof(cuDoubleComplex)*(2*K_SIM_SIZE));
  cudaMalloc((void**)&d_force, sizeof(double)*(3*Tpart));
  
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
      double KK = Kvec[kn];
      double val = 2. * KK / V;  //mult by 2 for symmetry

      double *K;
      K = new double[K_SIM_SIZE*3];
      K[0] = Kx; K[1] = Ky; K[2] = Kz;
      K[3] = -Kx; K[4] = Ky; K[5] = Kz;
      K[6] = Kx; K[7] = -Ky; K[8] = Kz;
      K[9] = -Kx; K[10] = -Ky; K[11] = Kz;

      //copy memory to device 
      cudaMemcpy(d_K, K, sizeof(double)*(K_SIM_SIZE*3), cudaMemcpyHostToDevice);
    
      //initialize accumulators to 0 (d_M does not accumulate)
      cudaMemset(d_sum, 0, sizeof(cuDoubleComplex)*(2*K_SIM_SIZE*blocks));
      cudaMemset(d_self, 0, sizeof(double)*blocks);

      // add kernel invocation
      dim3 dimGrid(blocks, 1, 1);
      dim3 dimBlock(BLOCK_WIDTH, 1, 1);
      kvec_add<<<dimGrid,dimBlock,sizeof(cuDoubleComplex)*2*K_SIM_SIZE*BLOCK_WIDTH>>>(d_pos, d_chrg, d_Nc, NC, d_K, d_sum, d_M, d_self);

      //copy output to host
      double slf_en[blocks];
      cuDoubleComplex sum_p[2*K_SIM_SIZE*blocks];
      cudaMemcpy(slf_en, d_self, sizeof(double)*blocks, cudaMemcpyDeviceToHost);
      cudaMemcpy(sum_p, d_sum, sizeof(cuDoubleComplex)*(2*K_SIM_SIZE*blocks), cudaMemcpyDeviceToHost);
      double self = 0.;
      for(int i = 0; i < blocks; ++i)
      	self += slf_en[i];
      cuDoubleComplex sum[2*K_SIM_SIZE];
      for(int b = 0; b < 2*K_SIM_SIZE; ++b){
	sum[b] = make_cuDoubleComplex(0., 0.);
      	for(int i = 0; i < blocks; ++i)
      	  sum[b] = cuCadd(sum[b], sum_p[2*K_SIM_SIZE*i + b]);
      }
      
      //compute energy - virial
      double Upar = 0, Wpar = 0;
      for(int b = 0; b < 4; ++b){
        double norm = cuCreal(sum[b])*cuCreal(sum[b]) + cuCimag(sum[b])*cuCimag(sum[b]);
        Upar += val * 0.5 * (norm - self);
        Wpar -= val * cuCimag(cuCmul(sum[b+K_SIM_SIZE], sum[b]));
      }
      if((nx == 0 && ny == 0) || (nx == 0 && nz == 0) || (ny == 0 && nz == 0)){
        Upar /= 4.;
        Wpar /= 4.;
      }else if((nz == 0 && nx != 0 && ny != 0) || (ny == 0 && nz != 0 && nx != 0) || (nx == 0 && ny != 0 && nz != 0)){
        Upar /= 2.;
        Wpar /= 2.;
      }
      tot_Uch += Upar;       
      tot_Wch += Wpar;
      
      // //
      // printf("%.12e %.12e\n",Upar,Wpar);
      // //

      // //copy memory to device
      cudaMemcpy(d_sumc, sum, sizeof(cuDoubleComplex)*(2*K_SIM_SIZE), cudaMemcpyHostToDevice);
      cudaMemcpy(d_force, F, sizeof(double)*(3*Tpart), cudaMemcpyHostToDevice);
      
      // force kernel invocation 
      k_force_kernel<<<dimGrid,dimBlock>>>(d_chrg, d_Nc, NC, d_sumc, d_M, d_K, d_force, val, nx, ny, nz);
      
      //copy output to host
      cudaMemcpy(F, d_force, sizeof(double)*(3*Tpart), cudaMemcpyDeviceToHost);


      // //
      // for(int l = 0; l < Tpart; ++l){
      //   for(int b = 0; b < K_SIM_SIZE; ++b)
      //     printf("(%f,%f) ",cuCreal(M[l*K_SIM_SIZE+b]),cuCimag(M[l*K_SIM_SIZE+b]));
      //   printf("\n");
      // }
      // printf("\n");
      // // //
      // for(int b = 0; b < 2*K_SIM_SIZE; ++b)
      // 	printf("(%f,%f) ",cuCreal(sum[b]),cuCimag(sum[b]));
      // printf("\n");
      // // //
      // printf("%f\n",self[0]);
      // //

      delete[] K;      
    }

  }
  cudaFree(d_K);
  cudaFree(d_pos);
  cudaFree(d_chrg);
  cudaFree(d_Nc);
  cudaFree(d_M);
  cudaFree(d_self);
  cudaFree(d_sum);
  cudaFree(d_sumc);
  cudaFree(d_force);
      
  
  //part.set_F(F, 3*Tpart);

  // //
  // for(int i = 0; i < Tpart; ++i){
  //   for(int a = 0; a < 3; ++a)
  //     cout<<F[3*i+a]<<" ";
  //   cout<<endl;
  // }
  // //
  
  // //
  // for(int i = 0; i < Tpart; ++i){
  //   int n, ni;
  //   part.get_kind(i, n, ni);
  //   printf("%f %f %f\n",part.get_F(n,ni,0),part.get_F(n,ni,1),part.get_F(n,ni,2));
  // }
  // printf("\n");
  // //
}



//------------------------------------------------------------------------------------



//MOVE X------------------------------------------------------------------------------

__global__ void moveX_kernel(double *pos, const double *mom, const double *mass, const int *Npart, int NC, double Dt){
  int n = threadIdx.x + blockDim.x * blockIdx.x;
  
  int Tpart = 0;
  for(int l = 0; l < NC; ++l)
    Tpart += Npart[l];

  if(n < Tpart){ 
    int k = 0;
    int Tt = 0;
    for(int m = 0; m < NC; ++m){
      Tt += Npart[m];
      if(n < Tt)
	break;
      k++;
    }
    double mn = mass[k];
    
    for(int a = 0; a < 3; ++a){
      pos[n*3+a] += Dt * mom[n*3+a] / mn;
    }
  }
}

void dev_moveX(Particles& part, double Dt){
  int NC = part.get_Nkinds();
  int Tpart = 0;
  for(int i = 0; i < part.get_Nkinds(); ++i)
    Tpart += part.get_N(i);

  //Set device
  cudaSetDevice(0);
  
  double *X = part.get_X();
  double *P = part.get_P();
  double *M = part.get_M();
  int *Nprt = part.get_Npart();
  
  //Allocate and copy memory to device
  double *d_pos, *d_mom, *d_mass;
  int *d_Nc;
  cudaMalloc((void**)&d_pos, sizeof(double)*(3*Tpart));
  cudaMemcpy(d_pos, X, sizeof(double)*(3*Tpart), cudaMemcpyHostToDevice);
  cudaMalloc((void**)&d_mom, sizeof(double)*(3*Tpart));
  cudaMemcpy(d_mom, P, sizeof(double)*(3*Tpart), cudaMemcpyHostToDevice);
  cudaMalloc((void**)&d_mass, sizeof(double)*NC);
  cudaMemcpy(d_mass, M, sizeof(double)*NC, cudaMemcpyHostToDevice);
  cudaMalloc((void**)&d_Nc, sizeof(double)*NC);
  cudaMemcpy(d_Nc, Nprt, sizeof(double)*NC, cudaMemcpyHostToDevice);

  //Define grid and block size
  int blocks = Tpart / BLOCK_WIDTH;
  if(Tpart % BLOCK_WIDTH) blocks++;

  // Kernel invocation code
  dim3 dimGrid(blocks, 1, 1);
  dim3 dimBlock(BLOCK_WIDTH, 1, 1);
  moveX_kernel<<<dimGrid,dimBlock>>>(d_pos, d_mom, d_mass, d_Nc, NC, Dt);

  //copy output to host
  cudaMemcpy(X, d_pos, sizeof(double)*(3*Tpart), cudaMemcpyDeviceToHost);
  
  cudaFree(d_pos);
  cudaFree(d_mom);
  cudaFree(d_mass);
  cudaFree(d_Nc);
}

//------------------------------------------------------------------------------------


//MOVE P------------------------------------------------------------------------------

__global__ void moveP_kernel(double *mom, const double *force, const double *mass, const int *Npart, int NC, double *K, double thermo_s, double Dt, bool last){
  int n = threadIdx.x + blockDim.x * blockIdx.x;

  __shared__ double sum_k[BLOCK_WIDTH];
  __syncthreads();
  
  double Kc0 = 0.;

  int Tpart = 0;
  for(int l = 0; l < NC; ++l)
    Tpart += Npart[l];
  
  if(n < Tpart){
    int k = 0;
    int Tt = 0;
    for(int m = 0; m < NC; ++m){
      Tt += Npart[m];
      if(n < Tt)
	break;
      k++;
    }
    double mn = mass[k];
    double Dt2 = Dt / 2;
    double psq = 0.;
    
    if(!last){
      double sfac = 1. + Dt2 * thermo_s;
      for(int a = 0; a < 3; ++a){
	mom[n*3+a] = (mom[n*3+a] + Dt2 * force[n*3+a]) / sfac;
	psq += mom[n*3+a] * mom[n*3+a];
      }
    }else{
      double sfac = 1. - Dt2 * thermo_s;
      for(int a = 0; a < 3; ++a){
	mom[n*3+a] = mom[n*3+a] * sfac + Dt2 * force[n*3+a];
	psq += mom[n*3+a] * mom[n*3+a];
      }
    }
      
    Kc0 += psq / mn;
  }

  sum_k[threadIdx.x] = Kc0;
  __syncthreads();
  if(threadIdx.x == 0){ 
    K[blockIdx.x] = 0.;
    for(int i = 0; i < BLOCK_WIDTH; ++i) K[blockIdx.x] += sum_k[i];
  }
}


void dev_moveP(Particles& part, const Thermostat& thermo, double Dt, double& K, bool last){
  int NC = part.get_Nkinds();
  int Tpart = 0;
  for(int i = 0; i < part.get_Nkinds(); ++i)
    Tpart += part.get_N(i);

  double thermo_s = thermo.get_coo();
  
  //Set device
  cudaSetDevice(0);
  
  double *P = part.get_P();
  double *M = part.get_M();
  double *F = part.get_F();
  int *Nprt = part.get_Npart();
  
  //Allocate and copy memory to device
  double *d_mom, *d_mass, *d_force, *d_K;
  int *d_Nc;
  cudaMalloc((void**)&d_mom, sizeof(double)*(3*Tpart));
  cudaMemcpy(d_mom, P, sizeof(double)*(3*Tpart), cudaMemcpyHostToDevice);
  cudaMalloc((void**)&d_mass, sizeof(double)*NC);
  cudaMemcpy(d_mass, M, sizeof(double)*NC, cudaMemcpyHostToDevice);
  cudaMalloc((void**)&d_force, sizeof(double)*(3*Tpart));
  cudaMemcpy(d_force, F, sizeof(double)*(3*Tpart), cudaMemcpyHostToDevice);
  cudaMalloc((void**)&d_Nc, sizeof(double)*NC);
  cudaMemcpy(d_Nc, Nprt, sizeof(double)*NC, cudaMemcpyHostToDevice);

  //Define grid and block size
  int blocks = Tpart / BLOCK_WIDTH;
  if(Tpart % BLOCK_WIDTH) blocks++;

  //Allocate memory for outputx
  cudaMalloc((void**)&d_K, sizeof(double)*blocks);
  
  // Kernel invocation code
  dim3 dimGrid(blocks, 1, 1);
  dim3 dimBlock(BLOCK_WIDTH, 1, 1);
  moveP_kernel<<<dimGrid,dimBlock>>>(d_mom, d_force, d_mass, d_Nc, NC, d_K, thermo_s, Dt, last);

  //copy output to host
  double *Ke = new double[blocks];
  cudaMemcpy(P, d_mom, sizeof(double)*(3*Tpart), cudaMemcpyDeviceToHost);
  cudaMemcpy(Ke, d_K, sizeof(double)*blocks, cudaMemcpyDeviceToHost);
  K = 0.;
  for(int i = 0; i < blocks; ++i)
    K += Ke[i];
  K /= 2;
  
  delete[] Ke;
  cudaFree(d_force);
  cudaFree(d_mom);
  cudaFree(d_mass);
  cudaFree(d_Nc);
  cudaFree(d_K);
}
  
//------------------------------------------------------------------------------------




//CORELATIONS-------------------------------------------------------------------------

__global__ void correl_kernel(const double *pos, double *hist, const int *Npart, int NC, const double* L, int n_bars, double bar_w){
  int n = threadIdx.x + blockDim.x * blockIdx.x;
  
  int Tpart = 0;
  for(int l = 0; l < NC; ++l)
    Tpart += Npart[l];

  if(n < Tpart){
    int k1 = 0;
    int Tt = 0;
    for(int l = 0; l < NC; ++l){
      Tt += Npart[l];
      if(n < Tt)
	break;
      k1++;
    }
    
    //even out the charge among threads
    int mx = static_cast<int>(ceil(static_cast<float>(Tpart-1)/2));
    if(fmod(static_cast<float>(Tpart),static_cast<float>(2.)) == 0. && n >= Tpart/2)
      mx = static_cast<int>(floor(static_cast<float>(Tpart-1)/2));
    
    int i = n+1 - Tpart*(static_cast<int>(floor((n+1)/Tpart + 0.5)));
    int cnt = 0;
    while(cnt < mx){
      int k2 = 0;
      Tt = 0;
      for(int l = 0; l < NC; ++l){
	Tt += Npart[l];
	if(i < Tt)
	  break;
	k2++;
      }
      
      double rij[3];
      for(int a = 0; a < 3; ++a){
	// //
	// printf("%f-%f\n", pos[n*3+a], pos[i*3+a]);
	// //
	rij[a] = pos[n*3+a] - pos[i*3+a];
	rij[a] -= L[a] * floor(rij[a] / L[a] + 0.5); 
      }    

      double RIJsq = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
      double RIJ = sqrt(RIJsq);
      
      int bin = static_cast<int>(floor(RIJ / bar_w));
      if(bin < n_bars){
	int mx_k = max(k1,k2);
	int mn_k = min(k1,k2);
	int pos_hist = 0;
	for(int l = 0; l < mn_k; ++l)
	  pos_hist += NC - l;
	pos_hist += mx_k - mn_k;
	pos_hist *= n_bars;
	pos_hist += bin;

	double val = 1.;
	if(k1 == k2) val = 2.; 
	atomicAdd(&hist[pos_hist], val); //000
      }
      i += 1 - Tpart*static_cast<int>(floor((i+1)/Tpart + 0.5));
      cnt++;
    }
    
  }
}


void dev_acc_correl(const Particles& part, Hist& histograms, const double L[3], double bar_w){
  int NC = part.get_Nkinds();
  //int NS = (NC-1)*(NC-1) - ((NC-1)*(NC-2))/2;
  int Tpart = 0;
  for(int i = 0; i < part.get_Nkinds(); ++i)
    Tpart += part.get_N(i);

  int n_bars = histograms.get_nbars();

  //Set device
  cudaSetDevice(0);
  
  double *X = part.get_X();
  int *Nprt = part.get_Npart();
  
  //Allocate and copy memory to device
  double *d_pos, *d_hist, *d_L;
  int *d_Nc;
  cudaMalloc((void**)&d_pos, sizeof(double)*(3*Tpart));
  cudaMemcpy(d_pos, X, sizeof(double)*(3*Tpart), cudaMemcpyHostToDevice);
  cudaMalloc((void**)&d_Nc, sizeof(double)*NC);
  cudaMemcpy(d_Nc, Nprt, sizeof(double)*NC, cudaMemcpyHostToDevice);
  cudaMalloc((void**)&d_L, sizeof(double)*3);
  cudaMemcpy(d_L, L, sizeof(double)*3, cudaMemcpyHostToDevice);

  //Define grid and block size
  int blocks = Tpart / BLOCK_WIDTH;
  if(Tpart % BLOCK_WIDTH) blocks++;

  //Allocate memory for output
  int hsize = histograms.get_size();
  double *hist = histograms.get_hist();
  cudaMalloc((void**)&d_hist, sizeof(double)*hsize);
  cudaMemcpy(d_hist, hist, sizeof(double)*hsize, cudaMemcpyHostToDevice);


  // Kernel invocation code
  dim3 dimGrid(blocks, 1, 1);
  dim3 dimBlock(BLOCK_WIDTH, 1, 1);
  correl_kernel<<<dimGrid,dimBlock>>>(d_pos, d_hist, d_Nc, NC, d_L, n_bars, bar_w);

  //copy output to host
  cudaMemcpy(hist, d_hist, sizeof(double)*hsize, cudaMemcpyDeviceToHost);
  histograms.set_hist(hist, hsize);
  
  cudaFree(d_hist);
  cudaFree(d_pos);
  cudaFree(d_L);
  cudaFree(d_Nc);
}


//------------------------------------------------------------------------------------
