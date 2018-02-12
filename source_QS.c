/*
  General temporal cycle (Externally injected bunch + laser + plasma (fluid or PIC))

  - depositBeamCharge()
  - swipeParticles()
  - dump data (beam, laser, plasma)
  - pushParticle()
  - evolveLaserPulse() 

  * Define FLUID for fluid calculations
  * Define _VERBOSE_ for verbose output
  * If FLUID is not defined the calculation is PIC, then either _LINEAR_SHAPE_FUNCTION_ or _QUADRATIC_SHAPE_FUNCTION_ must be defined
  * If evolveLaserPulse is called do not forget to properly set the value self_consistent_laser (default is NO)

  NB: the externally injected beam(s) DOES(DO) NOT interact with the laser 

  Restart: all the functions are at the end
  - loadLaser (also set the propagation variable s)
  - loadLaserFromINFERNO (also set the propagation variable s)
  - loadBeam
  - loadBeamFromINFERNO
  
*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<complex.h>

#define MAX(x, y) ((x)>(y)?(x):(y))
#define MIN(x, y) ((x)<(y)?(x):(y))
#define SGN(x) ((x)>0?(+1):(-1))

#define _CHARGE_ELECTRON_ (-1)
#define _CHARGE_POSITRON_ (+1)
#define Q0 1.6021892e-19 // electron charge, SI [C]
#define Q0_cgs 4.8e-10 // electron charge, cgs [statC]
#define R0_cgs (2.8179403267e-13) // classical electron radius, cgs [cm]
#define ORIGINAL_LOCATION -9999999999.
#define ALL_PARTICLES -9999
#define ASCII 1
#define BINARY 2

#define VOL0_QUADRATIC (13./64.)
#define VOL1_QUADRATIC (193./192.)
#define VOL0_GEO (1./8.)

#define YES 1
#define NO 0

#ifdef PARAXIAL_APPROXIMATION
#define _NONPARAXIAL_TERMS_ 0
#else
#define _NONPARAXIAL_TERMS_ 1
#endif

// grid definition
int Nz, Nr;
int beamDriver=NO;
double zmin, zmax, dz;
double rmax, dr;
double smin=0, smax=0, s=0, ds=0; // time integration
char *DIRECTORY_OUTPUT;

// numerical parameters
#ifdef FLUID
#ifdef LINEAR
int computeAllTheFields=NO;
#endif
#endif

#ifdef PIC
int filterDensityForOutput=YES;
#else
int filterDensityForOutput=NO;
#endif

int filterProperDensity=YES, Npass_filterProperDensity=4;
int filterLaserIntensity=YES, Npass_filterLaserIntensity=4;
int filterTransverseCurrent=NO, Npass_filterTransverseCurrent=4;

double alpha_B=0.98; 
double alpha_A=1; 
double psi_min_0=-0.990;
double d_psi_min= 0.003;
double EPSI=.1;


// fields
#define NZ_GHOST 2
double **density, **Jz, **Jr, **Ez, **Bphi, **Er_Bphi, **rho, **B2, **rhob, **Jz_ionization;
double **ur, **uz, **psi, **gamma_; // only for fluid

//laser
int RESET_evolveLaserPulse=NO;
double **a2_half, **dr_a2_half;
complex **a, **a_old, **tmp_a;

//beam
int Npart=0;
int *self_consistent=NULL, *particle_beam_active=NULL;
double *xpb=NULL, *ypb=NULL, *zpb=NULL, *uxb=NULL, *uyb=NULL, *uzb=NULL, *q=NULL, *me_mb=NULL; 

// plasma
#define AB_order 5 // order of the Adam-Bashfort integrator
#if AB_order==5
double AB_coefficients[]={1901./720., -1387./360., 109./30., -637./360., 251./720.}; // coefficients 5th-order Adams-Bashfort
double AM_coefficients[]={251./720. , 646./720.  ,-264./720., 106./720., -19./720.}; // coefficients 5th-order Adams-Moulton
#endif

int Nptot, Nppc;
int *active, *active_p;
double *rp, *urp, *psip, *q0, *rp0, *rp_p, *urp_p, *psip_p;
double *d_rp[AB_order], *d_urp[AB_order], *d_psip[AB_order];
double *tmp_d_urp;

// laser driver
int self_consistent_laser=NO;
int n_supergaussian;
double kpz0, kpzf, k0_kp, a0, kpL, kpW;
double delta_g;

double (*backgroundDensityProfile) (double, double);
complex (*laserEnvelope) (double, double, double)=NULL;

void numericalParameters()
{
#ifdef FLUID
#ifndef LINEAR
  printf("FLUID plasma [NONLINEAR]\n");
#else  
#ifdef ACCURATE
  printf("FLUID plasma [LINEAR, ACCURATE]\n");
#else
  printf("FLUID plasma [LINEAR, NOT ACCURATE]\n");
#endif
#endif
#endif

#ifdef PIC
  printf("PIC plasma\n");

#ifdef _LINEAR_SHAPE_FUNCTION_
  printf("grid<->particle: linear shape function\n");
#endif

#ifdef _QUADRATIC_SHAPE_FUNCTION_
  printf("grid<->particle: quadratic shape function\n");
#endif

#endif

  printf("alpha_A: %f\n",alpha_A);
  printf("alpha_B: %f\n",alpha_B);
#ifdef PIC
  printf("psi_min_0: %f (%f)\n", psi_min_0, d_psi_min);
#endif
  printf("EPSI:%f\n", EPSI);

  printf("Filters:\n");
  printf("- filterProperDensity=%i [Npass_filterProperDensity=%i]\n", filterProperDensity, Npass_filterProperDensity);
  printf("- filterLaserIntensity=%i [Npass_filterLaserIntensity=%i]\n", filterLaserIntensity, Npass_filterLaserIntensity);
  printf("- filterTransverseCurrent=%i [Npass_filterTransverseCurrent=%i]\n", filterTransverseCurrent, Npass_filterTransverseCurrent);
  printf("- filterDensityForOutput=%i\n", filterDensityForOutput);
}

int compute_Nz()
{
  Nz=(int)((zmax-zmin)/dz+0.5);
  dz=(zmax-zmin)/(Nz-1);
}

int compute_Nr()
{
  Nr=(int)(rmax/dr+0.5);
  dr=rmax/(Nr-1);
}

void binomialFilter(double *data, int ndata, double alpha)
{
  int k;
  static int ndata_old=0;
  static double *unfiltered_data=0;

  if (ndata_old!=ndata)
    {
      unfiltered_data=(double*)realloc((void*)unfiltered_data, ndata*sizeof(double));
      ndata_old=ndata;
    }
  
  memcpy((void*)unfiltered_data, data, ndata*sizeof(double));

  for (k=1; k<ndata-1; k++)
    data[k]=alpha*unfiltered_data[k]+0.5*(1-alpha)*(unfiltered_data[k-1]+unfiltered_data[k+1]);
  
}

void filter4xB_C(double *data, int ndata, int Npass)
{
  /*
    Radial filter: 4x(1,2,1) + compensator
    N.B. Assumes data are EVEN for r->-r (e.g., density, proper density, Jz, etc..)
  */

  int k;
  static int ndata_old=0;
  static double *data_buffer=0;

  if (ndata_old!=ndata)
    {
      data_buffer=(double*)realloc((void*)data_buffer, (2*ndata-1)*sizeof(double));
      ndata_old=ndata;
    }

  for (k=0; k<ndata-1; k++)
    data_buffer[k]=data[ndata-1-k];
  memcpy(&data_buffer[ndata-1], data, ndata*sizeof(double));
  
  for (k=0; k<Npass; k++)
    binomialFilter(data_buffer, 2*ndata-1, 0.5);
  binomialFilter(data_buffer, 2*ndata-1, 0.5*Npass+1); // compensator

  memcpy(data, &data_buffer[ndata-1], ndata*sizeof(double));
}

void filterNxB_C_odd(double *data, int ndata, int Npass)
{
  /*
    Radial filter: Npass x (1,2,1) + compensator
    N.B. Assumes data are ODD for r->-r (e.g., Jr, Er, etc..)
  */

  int k;
  static int ndata_old=0;
  static double *data_buffer=0;

  if (ndata_old!=ndata)
    {
      data_buffer=(double*)realloc((void*)data_buffer, (2*ndata-1)*sizeof(double));
      ndata_old=ndata;
    }

  for (k=0; k<ndata-1; k++)
    data_buffer[k]=-data[ndata-1-k];
  memcpy(&data_buffer[ndata-1], data, ndata*sizeof(double));
  
  for (k=0; k<Npass; k++)
    binomialFilter(data_buffer, 2*ndata-1, 0.5);
  binomialFilter(data_buffer, 2*ndata-1, 0.5*Npass+1); // compensator

  memcpy(data, &data_buffer[ndata-1], ndata*sizeof(double));
}

double **allocatePointer(int Nz_ghost)
{
  int i;
  double **pointer;

  pointer=(double**)malloc((Nz+Nz_ghost)*sizeof(double*));
  for (i=0; i<Nz+Nz_ghost; i++)
    {
      pointer[i]=(double*)malloc(Nr*sizeof(double));
      memset((void*)pointer[i], 0, Nr*sizeof(double));
    }
  
  return pointer;
}

complex **allocatePointer_complex(int Nz_ghost)
{
  int i;
  complex **pointer;

  pointer=(complex**)malloc((Nz+Nz_ghost)*sizeof(complex*));
  for (i=0; i<Nz+Nz_ghost; i++)
    {
      pointer[i]=(complex*)malloc(Nr*sizeof(complex));
      memset((void*)pointer[i], 0, Nr*sizeof(complex));
    }
  
  return pointer;
}

complex zeroLaserField(double z, double r, double s)
{
  return 0.;
}

void allocateFields()
{
  int i, j;

  density=allocatePointer(NZ_GHOST);
  Jz=allocatePointer(NZ_GHOST);
  Jr=allocatePointer(NZ_GHOST);
  rho=allocatePointer(NZ_GHOST);
  rhob=allocatePointer(NZ_GHOST);
  Jz_ionization=allocatePointer(NZ_GHOST);
  B2=allocatePointer(NZ_GHOST);

  Er_Bphi=allocatePointer(NZ_GHOST);
  Ez=allocatePointer(NZ_GHOST);
  Bphi=allocatePointer(NZ_GHOST);

#ifdef FLUID
  ur=allocatePointer(NZ_GHOST);
  uz=allocatePointer(NZ_GHOST);
  psi=allocatePointer(NZ_GHOST);
  gamma_=allocatePointer(NZ_GHOST);
#endif

  a=allocatePointer_complex(NZ_GHOST);
  a_old=allocatePointer_complex(NZ_GHOST);
  tmp_a=allocatePointer_complex(NZ_GHOST);
  a2_half=allocatePointer(NZ_GHOST);
  dr_a2_half=allocatePointer(NZ_GHOST);
  
  if (laserEnvelope==NULL)
    laserEnvelope=zeroLaserField;

  for (i=0; i<Nz; i++)
    {
      for (j=0; j<Nr; j++)
	{
	  a_old[i][j]=laserEnvelope(zmin+i*dz, j*dr, smin-ds);
	  a[i][j]    =laserEnvelope(zmin+i*dz, j*dr, smin   );  
	  a2_half[i][j]=0.5*a[i][j]*conj(a[i][j]);
	}
      
      for (j=0; j<Nr-1; j++)
	dr_a2_half[i][j]=(a2_half[i][j+1]-a2_half[i][abs(j-1)])/(2*dr);
    }
}

void allocateParticles()
{
  int j, k, n;
  double r_min, r0;

#ifdef FLUID
  printf("... WARNING: function [allocateParticles()] ignored!!\n\n");
  return;
#endif

  Nptot=Nppc*(Nr-1);
  active  =(int*)malloc(Nptot*sizeof(int));
  active_p=(int*)malloc(Nptot*sizeof(int));
  rp  =(double*)malloc(Nptot*sizeof(double));
  urp =(double*)malloc(Nptot*sizeof(double));
  psip=(double*)malloc(Nptot*sizeof(double));
  rp_p  =(double*)malloc(Nptot*sizeof(double));
  urp_p =(double*)malloc(Nptot*sizeof(double));
  psip_p=(double*)malloc(Nptot*sizeof(double));
  q0  =(double*)malloc(Nptot*sizeof(double));
  rp0 =(double*)malloc(Nptot*sizeof(double));
  tmp_d_urp  =(double*)malloc(Nptot*sizeof(double));

  n=0;
  for (j=0; j<Nr-1; j++)
    {
      r_min=j*dr;
      r0=r_min+0.5*dr;
      
      for (k=0; k<Nppc; k++)
	{
	  rp0[n]=r_min+dr*(k+0.5)/(double)Nppc;
	  n++;
	}
    }
  
  for (k=0; k<AB_order; k++)
    {
      d_rp[k]  =(double*)malloc(Nptot*sizeof(double));
      d_urp[k] =(double*)malloc(Nptot*sizeof(double));
      d_psip[k]=(double*)malloc(Nptot*sizeof(double));
      
      memset((void*)d_rp[k]  , 0, Nptot*sizeof(double));
      memset((void*)d_urp[k] , 0, Nptot*sizeof(double));
      memset((void*)d_psip[k], 0, Nptot*sizeof(double));
    }
}

void solveTridiag(double *a, double *b, double *c, double *r, double *y, int N)
{
  int i;
  double buffer;

  for (i=1; i<N; i++)
    {
      buffer=a[i]/b[i-1];
      
      b[i]-=buffer*c[i-1];
      r[i]-=buffer*r[i-1];
    }

  y[N-1]=r[N-1]/b[N-1];
  for (i=N-2; i>=0; i--)
    y[i]=(r[i]-c[i]*y[i+1])/b[i];
}

void solveScalarField(double *field, double *source)
{
  /*
    for a given zeta location, this routine solves for fied, 
    such that 

       \nabla^2 field = source
  */
  
  int j;
  static int FIRST_TIME=1;
  static double *a, *b, *c, *r;

  if (FIRST_TIME)
    {
      a=(double*)malloc(Nr*sizeof(double));
      b=(double*)malloc(Nr*sizeof(double));
      c=(double*)malloc(Nr*sizeof(double));
      r=(double*)malloc(Nr*sizeof(double));
      
      FIRST_TIME=0;
    }

  a[0]=0; // not used
  b[0]=-1;
  c[0]=1;
  r[0]=dr*dr*source[0]/4.;
  for (j=1; j<Nr; j++)
    {
      a[j]=1.-1./(2.*j);
      b[j]=-2.;
      c[j]=1.+1./(2.*j);
      r[j]=dr*dr*source[j];
    }

  solveTridiag(a, b, c, r, field, Nr);
}

void solveBesselEquation(double *field, double *source)
{
  int j, k;
  static int FIRST_TIME=1;
  static double *a, *b, *c, *r;

  /*
    for a given zeta location, this routine solves for fied, 
    such that 

       (\nabla^2 - EPSI) field = source
  */
  
  if (FIRST_TIME)
    {
      a=(double*)malloc(Nr*sizeof(double));
      b=(double*)malloc(Nr*sizeof(double));
      c=(double*)malloc(Nr*sizeof(double));
      r=(double*)malloc(Nr*sizeof(double));

      FIRST_TIME=0;
    }

  a[0]=0; // not used
  b[0]=-1-dr*dr*EPSI/4.;
  c[0]=1;
  r[0]=dr*dr*source[0]/4.;
  for (j=1; j<Nr; j++)
    {
      a[j]=1.-1./(2.*j);
      b[j]=-2-dr*dr*EPSI;
      c[j]=1.+1./(2.*j);
      r[j]=dr*dr*source[j];
    }
  
  solveTridiag(a, b, c, r, field, Nr);
}

#ifdef PIC
#ifdef _LINEAR_SHAPE_FUNCTION_
void compute_densityAndCurrent(double *d, double *rho_, double *Jz_, double *Jr_, int i, double *rp_, double *urp_, double *psip_, int *active_)
{
  /*
    Calculation for density and current at a given zeta location.
    Version /w linear shape functions.
  */

  int n, j;
  double w, gammap, uzp, buffer, r0, dvol;
  double sgn_r, int_a2_half, buf;

  memset((void*)d   , 0, Nr*sizeof(double));
  memset((void*)rho_, 0, Nr*sizeof(double));
  memset((void*)Jr_ , 0, Nr*sizeof(double));
  memset((void*)Jz_ , 0, Nr*sizeof(double));
  
  for (n=0; n<Nptot; n++)
    {
      if (!active_[n]) continue;
      
      sgn_r=SGN(rp_[n]);
      w=sgn_r*rp_[n]/dr;
      j=(int)w;
      if (j>=Nr) continue;
      w-=(double)j;

      int_a2_half=(1-w)*a2_half[i][j]+w*a2_half[i][j+1];
      uzp=(1+int_a2_half+urp_[n]*urp_[n]-pow(1+psip_[n],2))/(2*(1+psip_[n]));
      gammap=uzp+1+psip_[n];

      buffer=gammap*q0[n]/(1+psip_[n]);
      d[j  ]+=(1-w)*buffer;
      d[j+1]+=w*buffer;
      rho_[j  ]+=(1-w)*buffer/gammap;
      rho_[j+1]+=w*buffer/gammap;
      /*
	buf=(1+a2_half[i][j]  +urp_[n]*urp_[n]+pow(1+psip_[n],2))/(2*(1+psip_[n]));
	rho_[j  ]+=(1-w)*buffer/buf;
	buf=(1+a2_half[i][j+1]+urp_[n]*urp_[n]+pow(1+psip_[n],2))/(2*(1+psip_[n]));
	rho_[j+1]+=w*buffer/buf;
      */

      buffer=uzp*q0[n]/(1+psip_[n]);
      Jz_[j  ]-=(1-w)*buffer;
      Jz_[j+1]-=w*buffer;

      buffer=sgn_r*urp_[n]*q0[n]/(1+psip_[n]);
      Jr_[j  ]-=(1-w)*buffer;
      Jr_[j+1]-=w*buffer;
    }
    
  for (j=0; j<Nr; j++)
    {
      r0=j*dr;
      
      if (j==0)
	dvol=(1./6.)*2*M_PI*dr*dr*dz;
      else
	dvol=2*M_PI*r0*dr*dz;
      
      d[j]   /=dvol;
      rho_[j]/=dvol;
      Jz_[j] /=dvol;
      Jr_[j] /=dvol;
    }

  /*
  d[0]=(4./3.)*d[1]-(1./3.)*d[2];
  Jz_[0]=(4./3.)*Jz_[1]-(1./3.)*Jz_[2];
  */

  d[Nr-1]=rho_[Nr-1]=backgroundDensityProfile(s, rmax);
  Jz_[Nr-1]=0;
  Jr_[0]=Jr_[Nr-1]=0;

  for (j=0; j<Nr; j++)
    Jz_[j]+=rhob[i][j]+Jz_ionization[i][j];
  
}
#endif

#ifdef _QUADRATIC_SHAPE_FUNCTION_
void compute_densityAndCurrent(double *d, double *rho_, double *Jz_, double *Jr_, int i, double *rp_, double *urp_, double *psip_, int *active_)
{
  /*
    Calculation for density and current at a given zeta location.
    Version /w quadratic shape functions.
  */

  int n, j, jl, j1, jj;
  double w, wl, w2, w_[3], gammap, uzp, buffer1, buffer2, buffer3, r0, dvol, dvol2;
  double sgn_r, int_a2_half;
  
  memset((void*)d   , 0, Nr*sizeof(double));
  memset((void*)rho_, 0, Nr*sizeof(double));
  memset((void*)Jr_ , 0, Nr*sizeof(double));
  memset((void*)Jz_ , 0, Nr*sizeof(double));
  
  for (n=0; n<Nptot; n++)
    {
      if (!active_[n]) continue;
      
      sgn_r=SGN(rp_[n]);
      w=sgn_r*rp_[n]/dr;
      jl=(int)w;
      j=(int)(w+0.5);
      if (j<0 || j>=Nr) continue;
      wl=w-(double)jl;
      w-=(double)j;
      w2=w*w;
      w_[1]=0.75-w2;
      w_[2]=0.5*(0.25+w+w2);
      w_[0]=1-w_[1]-w_[2];

      int_a2_half=(1-wl)*a2_half[i][jl]+wl*a2_half[i][jl+1];
      uzp=(1+int_a2_half+urp_[n]*urp_[n]-pow(1+psip_[n],2))/(2*(1+psip_[n]));
      gammap=uzp+1+psip_[n];

      buffer1=gammap       *q0[n]/(1+psip_[n]);
      buffer2=uzp          *q0[n]/(1+psip_[n]);
      buffer3=sgn_r*urp_[n]*q0[n]/(1+psip_[n]);

      for (j1=0; j1<3; j1++)
	{
	  jj=MIN(j+j1-1, Nr-1);

	  d[abs(jj)]   +=w_[j1]*buffer1;
	  rho_[abs(jj)]+=w_[j1]*buffer1/gammap;
	  Jz_[abs(jj)]-=w_[j1]*buffer2;
	  Jr_[abs(jj)]-=w_[j1]*buffer3*SGN(jj);
	}
    }
    
  for (j=0; j<Nr; j++)
    {
      r0=j*dr;
      
      if (j==0)
	{
	  dvol=VOL0_QUADRATIC*2*M_PI*dr*dr*dz;
	  dvol2=VOL0_GEO*2*M_PI*dr*dr*dz;
	}
      else if (j==1)
	{
	  dvol=VOL1_QUADRATIC*2*M_PI*r0*dr*dz;
	  dvol2=2*M_PI*r0*dr*dz;
	}
      else
	dvol=dvol2=2*M_PI*r0*dr*dz;
      
      d[j]   /=dvol;
      rho_[j]/=dvol;
      Jz_[j] /=dvol;
      Jr_[j] /=dvol2;
    }

  /*
  d[0]=(4./3.)*d[1]-(1./3.)*d[2];
  Jz_[0]=(4./3.)*Jz_[1]-(1./3.)*Jz_[2];
  */

  d[Nr-1]=backgroundDensityProfile(s, rmax  );
  d[Nr-2]=backgroundDensityProfile(s, rmax-dr);
  Jz_[Nr-1]=Jz_[Nr-2]=0;
  Jr_[0]=Jr_[Nr-1]=Jr_[Nr-2]=0;

  for (j=0; j<Nr; j++)
    Jz_[j]+=rhob[i][j]+Jz_ionization[i][j];

}
#endif

void compute_dtParticle(int i, double *rp_, double *urp_, double *psip_, int *active_)
{
  int j, j1, jl, jj, n;
  double w, wl, w2, w_[3], uzp, gammap;
  double int_Ez, int_Er_Bphi, int_Bphi, buffer;
  double sgn_r, int_a2_half, int_dr_a2_half;

  for (n=0; n<Nptot; n++)
    {
      if (!active_[n]) continue;

      sgn_r=SGN(rp_[n]);

#ifdef  _LINEAR_SHAPE_FUNCTION_
      w=sgn_r*rp_[n]/dr;
      j=(int)w;
      if (j>=Nr) continue;
      w-=(double)j;
     
      int_Ez     =(1-w)*Ez[i][j]     +w*Ez[i][j+1];
      int_Er_Bphi=(1-w)*Er_Bphi[i][j]+w*Er_Bphi[i][j+1];
      int_Bphi   =(1-w)*Bphi[i][j]   +w*Bphi[i][j+1];      
      int_a2_half   =(1-w)*a2_half[i][j]+w*a2_half[i][j+1];
      int_dr_a2_half=(1-w)*dr_a2_half[i][j]+w*dr_a2_half[i][j+1];
      
#endif

#ifdef _QUADRATIC_SHAPE_FUNCTION_
      w=sgn_r*rp_[n]/dr;
      jl=(int)w;
      j=(int)(w+0.5);
      if (j<0 || j>=Nr) continue;
      wl=w-(double)jl;
      w-=(double)j;
      w2=w*w;
      w_[1]=0.75-w2;
      w_[2]=0.5*(0.25+w+w2);
      w_[0]=1-w_[1]-w_[2];

      int_Ez=int_Er_Bphi=int_Bphi=0;
      for (j1=0; j1<3; j1++)
	{
	  jj=MIN(j+j1-1, Nr-1);
	  
	  int_Ez     +=Ez[i][abs(jj)]*w_[j1];
	  int_Er_Bphi+=Er_Bphi[i][abs(jj)]*w_[j1]*SGN(jj);
	  int_Bphi   +=Bphi[i][abs(jj)]*w_[j1]*SGN(jj);
	}

      int_a2_half   =(1-wl)*a2_half[i][jl]+wl*a2_half[i][jl+1];
      int_dr_a2_half=(1-wl)*dr_a2_half[i][jl]+wl*dr_a2_half[i][jl+1];

#endif      
      int_Er_Bphi*=sgn_r;
      int_Bphi*=sgn_r;
      int_dr_a2_half*=sgn_r;

      uzp=(1+int_a2_half+urp_[n]*urp_[n]-pow(1+psip_[n],2))/(2*(1+psip_[n]));
      gammap=uzp+1+psip_[n];
      
      buffer=urp_[n]/(1+psip_[n]);
      d_rp[0][n]=-buffer;
      d_urp[0][n]=(0.5*int_dr_a2_half+gammap*int_Er_Bphi)/(1+psip_[n])+int_Bphi;
      d_psip[0][n]=buffer*(int_Er_Bphi)-int_Ez;
    }
}

void prepare_compute_dtParticle(int i, double *rp_, double *urp_, double *psip_, int *active_)
{
  int j, j1, jl, jj, n;
  double w, wl, w2, w_[3], uzp, gammap;
  double int_Ez, int_Er_Bphi, int_Bphi, buffer;
  double sgn_r, int_a2_half, int_dr_a2_half;

  for (n=0; n<Nptot; n++)
    {
      if (!active_[n]) continue;

      sgn_r=SGN(rp_[n]);

#ifdef  _LINEAR_SHAPE_FUNCTION_
      w=sgn_r*rp_[n]/dr;
      j=(int)w;
      if (j>=Nr) continue;
      w-=(double)j;
     
      int_Ez     =(1-w)*Ez[i][j]     +w*Ez[i][j+1];
      int_Er_Bphi=(1-w)*Er_Bphi[i][j]+w*Er_Bphi[i][j+1];
      int_Bphi   =(1-w)*Bphi[i][j]   +w*Bphi[i][j+1];      
      int_a2_half   =(1-w)*a2_half[i][j]+w*a2_half[i][j+1];
      int_dr_a2_half=(1-w)*dr_a2_half[i][j]+w*dr_a2_half[i][j+1];
      
#endif

#ifdef _QUADRATIC_SHAPE_FUNCTION_
      w=sgn_r*rp_[n]/dr;
      jl=(int)w;
      j=(int)(w+0.5);
      if (j<0 || j>=Nr) continue;
      wl=w-(double)jl;
      w-=(double)j;
      w2=w*w;
      w_[1]=0.75-w2;
      w_[2]=0.5*(0.25+w+w2);
      w_[0]=1-w_[1]-w_[2];

      int_Ez=int_Er_Bphi=int_Bphi=0;
      for (j1=0; j1<3; j1++)
	{
	  jj=MIN(j+j1-1, Nr-1);
	  
	  int_Ez     +=Ez[i][abs(jj)]*w_[j1];
	  int_Er_Bphi+=Er_Bphi[i][abs(jj)]*w_[j1]*SGN(jj);
	  int_Bphi   +=Bphi[i][abs(jj)]*w_[j1]*SGN(jj);
	}

      int_a2_half   =(1-wl)*a2_half[i][jl]+wl*a2_half[i][jl+1];
      int_dr_a2_half=(1-wl)*dr_a2_half[i][jl]+wl*dr_a2_half[i][jl+1];

#endif      
      int_Er_Bphi*=sgn_r;
      int_Bphi*=sgn_r;
      int_dr_a2_half*=sgn_r;

      uzp=(1+int_a2_half+urp_[n]*urp_[n]-pow(1+psip_[n],2))/(2*(1+psip_[n]));
      gammap=uzp+1+psip_[n];
      
      buffer=urp_[n]/(1+psip_[n]);
      d_rp[0][n]  =-buffer;
      tmp_d_urp[n]=(0.5*int_dr_a2_half+gammap*int_Er_Bphi)/(1+psip_[n]);
      d_psip[0][n]=buffer*(int_Er_Bphi)-int_Ez;
    }
}

void addBphi_compute_dtParticle(int i, double *rp_, double *urp_, double *psip_, int *active_)
{
  int j, j1, jj, n;
  double w, w2, w_[3], sgn_r, int_Bphi, buffer;
 
  for (n=0; n<Nptot; n++)
    {
      if (!active_[n]) continue;

      sgn_r=SGN(rp_[n]);

#ifdef  _LINEAR_SHAPE_FUNCTION_
      w=sgn_r*rp_[n]/dr;
      j=(int)w;
      if (j>=Nr) continue;
      w-=(double)j;
     
      int_Bphi=(1-w)*Bphi[i][j]+w*Bphi[i][j+1];      
#endif

#ifdef _QUADRATIC_SHAPE_FUNCTION_
      w=sgn_r*rp_[n]/dr;
      j=(int)(w+0.5);
      if (j<0 || j>=Nr) continue;
      w-=(double)j;
      w2=w*w;
      w_[1]=0.75-w2;
      w_[2]=0.5*(0.25+w+w2);
      w_[0]=1-w_[1]-w_[2];

      int_Bphi=0;
      for (j1=0; j1<3; j1++)
	{
	  jj=MIN(j+j1-1, Nr-1);
	  int_Bphi+=Bphi[i][abs(jj)]*w_[j1]*SGN(jj);
	}
#endif      
      int_Bphi*=sgn_r;
      d_urp[0][n]=tmp_d_urp[n]+int_Bphi;
    }
}
#endif

#define DF_DR(F, R, S) (-(F)/(R)+(S))
void integroB(double *f, double *source)
{
  int j;
  double fbuf, df[4];

  /*
    df/dr + f/r = source --> df/dr = -f/r + source
    
    f(r=0)=0
    df/dr(r=0)=source(0)/2
    
  */

  f[0]=0;
  for (j=0; j<Nr-1; j++)
    {
      if (j==0)
	df[0]=source[0]/2;
      else
	df[0]=DF_DR(f[j], j*dr, source[j]);
      
      fbuf=f[j]+df[0]*dr/2.;
      df[1]=DF_DR(fbuf, j*dr+0.5*dr, 0.5*(source[j]+source[j+1]));
      
      fbuf=f[j]+df[1]*dr/2.;
      df[2]=DF_DR(fbuf, j*dr+0.5*dr, 0.5*(source[j]+source[j+1]));

      fbuf=f[j]+df[2]*dr;
      df[3]=DF_DR(fbuf, j*dr+dr, source[j+1]);
     
      f[j+1]=f[j]+(df[0]+2*df[1]+2*df[2]+df[3])*dr/6.;
    }
  
}

// calcolo B con tridiagonale
void solveBTridiag(double *field, double *source)
{
  /*
    df/dr + f/r = source 
    
    f(r=0)=0
    df/dr(r=0)=source(0)/2
    
  */

  int j;
  static int FIRST_TIME=1;
  static double *a, *b, *c, *r;

  if (FIRST_TIME)
    {
      a=(double*)malloc(Nr*sizeof(double));
      b=(double*)malloc(Nr*sizeof(double));
      c=(double*)malloc(Nr*sizeof(double));
      r=(double*)malloc(Nr*sizeof(double));

      FIRST_TIME=0;
    }

  a[0]=0; // not used                                                                                                                                                              
  b[0]=1.0;
  c[0]=0.5;
  r[0]=dr*source[1];
  for (j=1; j<Nr-1; j++)
    {
      a[j]=-0.5;
      b[j]=1./(j+1);
      c[j]= 0.5;
      r[j]=dr*source[j+1];
    }

  field[0]=0;
  solveTridiag(a, b, c, r, field+1, Nr-1);
}

// calcolo B con recurrence
void solveBRecurrence(double *field, double *source)
{
  /*
    df/dr + f/r = source 
    
    f(r=0)=0
    df/dr(r=0)=source(0)/2
    
  */

  int j;

  field[0]=0;
  field[1]=0.5*dr*source[0];
  for (j=1; j<Nr-1; j++)
    field[j+1]=j*field[j]/(j+1.)+0.5*(source[j]+source[j+1])*(j+0.5)*dr/(j+1.);
}

void divergence(double *div_f, double *f)
{
  int j;

  div_f[0]=2*f[1]/dr;
  for(j=1; j<Nr-1; j++)
    div_f[j]=(f[j+1]-f[j-1])/(2*dr)+(0.25*f[j+1]+0.5*f[j]+0.25*f[j-1])/(j*dr);  // with transverse filtering
  div_f[Nr-1]=0;
}

void divergenceNoFilter(double *div_f, double *f)
{
  int j;

  div_f[0]=2*f[1]/dr;
  for(j=1; j<Nr-1; j++)
    div_f[j]=(f[j+1]-f[j-1])/(2*dr)+f[j]/(j*dr);  // with transverse filtering
  div_f[Nr-1]=0;
}

#ifdef PIC
void swipeParticle(int iter_max, double error_max)
{
  int i, j, k, n;
  double d_rp_tot, d_urp_tot, d_psip_tot, error, psi_min;
  static int FIRST_TIME=1;
  static double *source, *source2;

  if (FIRST_TIME)
    {
      source=(double*)malloc(Nr*sizeof(double));
      source2=(double*)malloc(Nr*sizeof(double));
      
      FIRST_TIME=0;
    }

  psi_min=psi_min_0;      
  
 start_iteration:
  
  for (n=0; n<Nptot; n++)
    {
      active[n]=1;
      rp[n]=rp0[n];
      urp[n]=0;
      psip[n]=0;
      q0[n]=(2*M_PI*dr*dz*backgroundDensityProfile(s, rp[n])/Nppc)*rp[n];
    }
  
  for (k=0; k<AB_order; k++)
    {      
      memset((void*)d_rp[k]  , 0, Nptot*sizeof(double));
      memset((void*)d_urp[k] , 0, Nptot*sizeof(double));
      memset((void*)d_psip[k], 0, Nptot*sizeof(double));
    }

  for (i=0; i<Nz+NZ_GHOST; i++)
    {
      memset((void*)Ez[i]     , 0, Nr*sizeof(double));
      memset((void*)Er_Bphi[i], 0, Nr*sizeof(double));
      memset((void*)Bphi[i]   , 0, Nr*sizeof(double));
      memset((void*)B2[i]     , 0, Nr*sizeof(double));      
    }
  
  compute_densityAndCurrent(density[Nz-1], rho[Nz-1], Jz[Nz-1], Jr[Nz-1], Nz-1, rp, urp, psip, active);
  if (filterTransverseCurrent==YES)
    filterNxB_C_odd(Jr[Nz-1], Nr, Npass_filterTransverseCurrent);

  for (i=Nz-1; i>=1; i--)
    { 
#ifdef _VERBOSE_
      printf("grid point: %i/%i, psi_min=%.4f \n", i, Nz-1, psi_min);
#endif      
      // Er-Bphi      
      for (j=0; j<Nr; j++)
	Er_Bphi[i][j]=4*Er_Bphi[i+1][j]/3-Er_Bphi[i+2][j]/3-2*dz*Jr[i][j]/3;
          
      memcpy((void*)Bphi[i], Bphi[i+1], Nr*sizeof(double));
      //memcpy((void*)B2[i]  , Bphi[i+1], Nr*sizeof(double)); // XXX serve? NO!

      prepare_compute_dtParticle(i, rp, urp, psip, active);
      
      for (n=0; n<Nptot; n++)
	{
	  if (!active[n]) continue;
	  
	  d_rp_tot=d_psip_tot=0;
	  for (j=0; j<AB_order; j++)
	    {
	      d_rp_tot+=AB_coefficients[j]*d_rp[j][n];
	      d_psip_tot+=AB_coefficients[j]*d_psip[j][n];
	    }
	  
	  rp_p[n]  =  rp[n]-d_rp_tot*dz;
	  psip_p[n]=psip[n]-d_psip_tot*dz;
	}
      
      for (k=0; k<iter_max; k++)
	{
	  for (j=0; j<Nr; j++)
	    Bphi[i][j]=alpha_B*Bphi[i][j]+(1-alpha_B)*B2[i][j];
	  
	  addBphi_compute_dtParticle(i, rp, urp, psip, active);

	  memcpy((void*)active_p, active, Nptot*sizeof(int));
	  
	  for (n=0; n<Nptot; n++)
	    {
	      if (!active[n]) continue;
	      
	      d_urp_tot=0;
	      for (j=0; j<AB_order; j++)
		d_urp_tot+=AB_coefficients[j]*d_urp[j][n];
	      
	      urp_p[n]=urp[n]-d_urp_tot*dz;
	  
	      if (!isfinite(urp_p[n]))
		{
		  psi_min+=d_psi_min;
		  printf("\n\n!!!!!!!! ERROR, NaN detected: redefining psi_min=%.3f\n\n\n", psi_min);
		  goto start_iteration;
		}

	      if (fabs(rp_p[n])>=rmax || psip_p[n]<psi_min)
		active_p[n]=0;
	    }
	  
	  compute_densityAndCurrent(density[i-1], rho[i-1], Jz[i-1], Jr[i-1], i-1, rp_p, urp_p, psip_p, active_p); //XXXX i-->i-1 in posizione 4
	  if (filterTransverseCurrent==YES)
	    filterNxB_C_odd(Jr[i-1], Nr, Npass_filterTransverseCurrent);

	  for (j=0; j<Nr; j++)
	    //source[j]=EPSI*Bphi[i][j]+(-3*Jr[i-1][j]+4*Jr[i][j]-Jr[i+1][j])/(2*dz);
	    source[j]=EPSI*(2*Bphi[i][j]-Bphi[i+1][j])+(-3*Jr[i-1][j]+4*Jr[i][j]-Jr[i+1][j])/(2*dz);
	  divergence(source2, source);    
	  for (j=0; j<Nr; j++)
	    source2[j]-=EPSI*Jz[i-1][j];
	  solveBesselEquation(Ez[i-1], source2);
	  for (j=0; j<Nr; j++)
	    Ez[i-1][j]=(4*Ez[i][j]-Ez[i+1][j]-2*dz*Ez[i-1][j])/3.;
	  
	  for (j=0; j<Nr; j++)
	    source[j]=(0.25*Jz[i-1][j]+0.50*Jz[i][j]+0.25*Jz[i+1][j])-(Ez[i+1][j]-Ez[i-1][j])/(2*dz);
	  //integroB(B2[i], source);      
	  //solveBTridiag(B2[i], source);
	  solveBRecurrence(B2[i], source);

	  error=-1;
	  for (j=0; j<Nr; j++)
	    error=MAX(error, fabs(Bphi[i][j]-B2[i][j]));
	  
	  if (error<error_max) break;
	}

#ifdef _VERBOSE_
      printf(" ... iterations: %i/%i (%e)\n", k, iter_max, error);
#endif      

#ifndef _DISABLE_WARNING_DIAGNOSTIC_
      if (k>=iter_max)
	printf("WARNING: for i=%i number of max iterations reached (iter_max=%i, error=%e)\n", i, iter_max, error);
#endif

      for (n=0; n<Nptot; n++)
	{
	  rp[n]=rp_p[n];
	  urp[n]=urp_p[n];
	  psip[n]=psip_p[n];
	  active[n]=active_p[n];
	  
	  if (!active[n]) continue;
	  
	  if (rp[n]<0)
	    {
	      rp[n]=-rp[n];
	      urp[n]=-urp[n];
		  
	      for (j=0; j<AB_order; j++)
		{
		  d_rp[j][n]*=-1;
		  d_urp[j][n]*=-1;
		}
	    }
	  else if (fabs(rp[n])>rmax)
	    active[n]=0;
	}

      for (j=4; j>0; j--)
	{
	  memcpy((void*)d_rp[j], d_rp[j-1], Nptot*sizeof(double));
	  memcpy((void*)d_urp[j], d_urp[j-1], Nptot*sizeof(double));
	  memcpy((void*)d_psip[j], d_psip[j-1], Nptot*sizeof(double));
	}
    }
  
  // --- values of the fields at z=zmin

  compute_densityAndCurrent(density[0], rho[0], Jz[0], Jr[0], 0, rp, urp, psip, active);
  if (filterTransverseCurrent==YES)
    filterNxB_C_odd(Jr[0], Nr, Npass_filterTransverseCurrent);

  for (j=0; j<Nr; j++)
    source[j]=EPSI*(2*Bphi[1][j]-Bphi[2][j])+(-3*Jr[0][j]+4*Jr[1][j]-Jr[2][j])/(2*dz);
  divergence(source2, source);    
  for (j=0; j<Nr; j++)
    source2[j]-=EPSI*Jz[0][j];
  solveBesselEquation(Ez[0], source2);
  for (j=0; j<Nr; j++)
    Ez[0][j]=(4*Ez[1][j]-Ez[2][j]-2*dz*Ez[0][j])/3.;

  for (j=0; j<Nr; j++)
    Er_Bphi[0][j]=4*Er_Bphi[1][j]/3-Er_Bphi[2][j]/3-2*dz*Jr[0][j]/3;

  for (j=0; j<Nr; j++)
    source[j]=Jz[0][j]-(-3*Ez[0][j]+4*Ez[1][j]-Ez[2][j])/(2*dz);
  //integroB(Bphi[0], source);      
  //solveBTridiag(Bphi[0], source);
  solveBRecurrence(Bphi[0], source);

  if (filterDensityForOutput==YES)
    for (i=0; i<Nz; i++)
      {
	binomialFilter(density[i], Nr, 0.5);
	binomialFilter(density[i], Nr, 1.5);
      }
}
#endif

#ifdef FLUID
#ifndef LINEAR
void swipeParticle(int iter_max, double error_max)
{
  int i, j, k;
  double buffer, error;
  static int FIRST_TIME=1;
  static double *source, *source2, *source_density, *source_density2, *source_ur;
  
  if (FIRST_TIME)
    {
      source=(double*)malloc(Nr*sizeof(double));
      source2=(double*)malloc(Nr*sizeof(double));
      source_density=(double*)malloc(Nr*sizeof(double));
      source_density2=(double*)malloc(Nr*sizeof(double));
      source_ur=(double*)malloc(Nr*sizeof(double));

      FIRST_TIME=0;
    }

  for (i=0; i<Nz+NZ_GHOST; i++)
    {
      memset((void*)Ez[i]     , 0, Nr*sizeof(double));
      memset((void*)Bphi[i]   , 0, Nr*sizeof(double));
      memset((void*)B2[i]     , 0, Nr*sizeof(double));      
      memset((void*)ur[i]     , 0, Nr*sizeof(double));      
      memset((void*)uz[i]     , 0, Nr*sizeof(double));      
      memset((void*)psi[i]    , 0, Nr*sizeof(double));      
      memset((void*)Jr[i]     , 0, Nr*sizeof(double));      
      memset((void*)Jz[i]     , 0, Nr*sizeof(double));      
      memset((void*)density[i], 0, Nr*sizeof(double));      
      memset((void*)gamma_[i] , 0, Nr*sizeof(double));      
    }

  i=Nz-1;
  for (j=0; j<Nr; j++)
    {
      density[i][j]=backgroundDensityProfile(s, j*dr);
      gamma_[i][j]=1+psi[i][j]+uz[i][j];
    }

  memcpy((void*)density[Nz  ], density[Nz-1], Nr*sizeof(double));
  memcpy((void*)density[Nz+1], density[Nz-1], Nr*sizeof(double));
  memcpy((void*)gamma_[Nz  ], gamma_[Nz-1], Nr*sizeof(double));
  memcpy((void*)gamma_[Nz+1], gamma_[Nz-1], Nr*sizeof(double));
  
  for (i=Nz-1; i>0; i--)
    {
#ifdef _VERBOSE_
      printf("grid point: %i/%i (%f)\n", i, Nz-1, density[i][0]);
#endif      
      // Er-Bphi      
      for (j=0; j<Nr; j++)
	Er_Bphi[i][j]=4*Er_Bphi[i+1][j]/3.-Er_Bphi[i+2][j]/3.-2*dz*Jr[i][j]/3.;

      for (j=0; j<Nr; j++)
	{
	  source[j]=-Jr[i][j];
	  source_density2[j]=density[i+1][j]*(1-uz[i+1][j]/gamma_[i+1][j]);
	}
      divergenceNoFilter(source_density, source); // source_density is the divergence of source (=beta_r*density)
      
      memcpy((void*)Bphi[i], Bphi[i+1], Nr*sizeof(double));
      memcpy((void*)B2[i]  , Bphi[i+1], Nr*sizeof(double));

      for (j=0; j<Nr; j++)
	psi[i-1][j]=psi[i+1][j]+2*dz*Ez[i][j];
      
      for (j=0; j<Nr; j++)
	source_ur[j]=(ur[i][j]*(ur[i][MIN(j+1, Nr-1)]-ur[i][abs(j-1)])/(2*dr*(j==Nr-1?0.5:1))+0.5*dr_a2_half[i][j]+gamma_[i][j]*Er_Bphi[i][j])/(1+psi[i][j]);

      for (k=0; k<iter_max; k++)
	{
	  for (j=0; j<Nr; j++)
	    {
	      Bphi[i][j]=alpha_B*Bphi[i][j]+(1-alpha_B)*B2[i][j];
	      
	      //buffer=(ur[i][j]*(ur[i][MIN(j+1, Nr-1)]-ur[i][abs(j-1)])/(2*dr*(j==Nr-1?0.5:1))+0.5*dr_a2_half[i][j]+gamma_[i][j]*Er_Bphi[i][j])/(1+psi[i][j])+Bphi[i][j];
	      buffer=source_ur[j]+Bphi[i][j];
	      ur[i-1][j]=ur[i+1][j]-2*dz*buffer;
	      
	      uz[i-1][j]=(1+a2_half[i-1][j]+ur[i-1][j]*ur[i-1][j]-pow(1+psi[i-1][j],2))/(2*(1+psi[i-1][j]));
	      gamma_[i-1][j]=1+psi[i-1][j]+uz[i-1][j];
	      
	      //density[i-1][j]=(density[i+1][j]*(1-uz[i+1][j]/gamma_[i+1][j])-2*dz*source_density[j])/(1-uz[i-1][j]/gamma_[i-1][j]);
	      density[i-1][j]=(source_density2[j]-2*dz*source_density[j])/(1-uz[i-1][j]/gamma_[i-1][j]);

	      buffer=-density[i-1][j]/gamma_[i-1][j];
	      Jr[i-1][j]=buffer*ur[i-1][j];
	      Jz[i-1][j]=buffer*uz[i-1][j]+rhob[i-1][j]+Jz_ionization[i-1][j];
	    }
	  
	  if (filterTransverseCurrent==YES)
	    filterNxB_C_odd(Jr[i-1], Nr, Npass_filterTransverseCurrent);
	  	  
	  for (j=0; j<Nr; j++)
	    //source[j]=EPSI*Bphi[i][j]+(-3*Jr[i-1][j]+4*Jr[i][j]-Jr[i+1][j])/(2*dz);
	    source[j]=EPSI*(2*Bphi[i][j]-Bphi[i+1][j])+(-3*Jr[i-1][j]+4*Jr[i][j]-Jr[i+1][j])/(2*dz);
	  divergence(source2, source);    
	  for (j=0; j<Nr; j++)
	    source2[j]-=EPSI*Jz[i-1][j];
	  solveBesselEquation(Ez[i-1], source2);
	  for (j=0; j<Nr; j++)
	    Ez[i-1][j]=(4*Ez[i][j]-Ez[i+1][j]-2*dz*Ez[i-1][j])/3.;
	  /*
	    instab se la box e' grande
	    divergence(source, Jr[i-1]);
	    solveScalarField(Ez[i-1], source);
	  */
	  for (j=0; j<Nr; j++)
	    source[j]=(0.25*Jz[i-1][j]+0.50*Jz[i][j]+0.25*Jz[i+1][j])-(Ez[i+1][j]-Ez[i-1][j])/(2*dz);
	  solveBRecurrence(B2[i], source);
	  
	  error=-1;
	  for (j=0; j<Nr; j++)
	    error=MAX(error, fabs(Bphi[i][j]-B2[i][j]));
	  
	  if (error<error_max) break;
	}

#ifdef _VERBOSE_
      printf(" ... iterations: %i/%i (%e)\n", k, iter_max, error);
#endif      

#ifndef _DISABLE_WARNING_DIAGNOSTIC_
      if (k>=iter_max)
	printf("WARNING: for i=%i number of max iterations reached (iter_max=%i, error=%e)\n", i, iter_max, error);
#endif
    }
  
  for (j=0; j<Nr; j++)
    Er_Bphi[0][j]=4*Er_Bphi[1][j]/3-Er_Bphi[2][j]/3-2*dz*Jr[0][j]/3;
  
  for (j=0; j<Nr; j++)
    source[j]=Jz[0][j]-(-3*Ez[0][j]+4*Ez[1][j]-Ez[2][j])/(2*dz);
  //integroB(Bphi[0], source);      
  //solveBTridiag(Bphi[0], source);
  solveBRecurrence(Bphi[0], source);

  for (i=0; i<Nz; i++)
    for (j=0; j<Nr; j++)
      rho[i][j]=density[i][j]/gamma_[i][j];

  if (filterDensityForOutput==YES)
    for (i=0; i<Nz; i++)
      {
	binomialFilter(density[i], Nr, 0.5);
	binomialFilter(density[i], Nr, 1.5);
      } 
}
#endif
//------------------------------ SOLVER FLUID LINEAR
#ifdef LINEAR
#ifndef ACCURATE
void solve_delta()
{
  int i, j;
  double buffer;

  // nabla(a^2) --> nabla(|a|^2/2)
  for (i=1; i<Nz-1; i++)
    {
      B2[i][0]=4*(a2_half[i][1]-a2_half[i][0])/(dr*dr)+(a2_half[i-1][0]-2*a2_half[i][0]+a2_half[i+1][0])/(dz*dz);
      for (j=1; j<Nr-1; j++)
	B2[i][j]=(a2_half[i][j+1]-2*a2_half[i][j]+a2_half[i][j-1])/(dr*dr)+(a2_half[i][j+1]-a2_half[i][j-1])/(2*j*dr*dr)+(a2_half[i+1][j]-2*a2_half[i][j]+a2_half[i-1][j])/(dz*dz);
    }

  for (j=0; j<Nr; j++)
    density[Nz][j]=density[Nz+1][j]=backgroundDensityProfile(s, j*dr);

  for (i=Nz; i>0; i--)
    for (j=0; j<Nr; j++)
      {
	density[i-1][j]=2*density[i][j]-density[i+1][j]+backgroundDensityProfile(s, j*dr)*(rhob[i][j]+Jz_ionization[i][j]+B2[i][j]/2.+backgroundDensityProfile(s, j*dr)-density[i][j])*dz*dz;
	//buffer=(backgroundDensityProfile(s, (j+1)*dr)-backgroundDensityProfile(s, fabs((j-1)*dr)))/(2*dr);
	//density[i-1][j]=2*density[i][j]-density[i+1][j]+dz*dz*(backgroundDensityProfile(s, j*dr)*(backgroundDensityProfile(s, j*dr)+rhob[i][j]+B2[i][j]/2.-density[i][j])+buffer*dr_a2_half[i][j]/2.);
      }
}

void solve_psi_noAnal()
{
  int i, j, k;
  static int FIRST_TIME=1;
  static double *a, *b, *c, *r;

  /*
    first solve for density perturbation (density)
  */
  
  if (FIRST_TIME)
    {
      a=(double*)malloc(Nr*sizeof(double));
      b=(double*)malloc(Nr*sizeof(double));
      c=(double*)malloc(Nr*sizeof(double));
      r=(double*)malloc(Nr*sizeof(double));

      FIRST_TIME=0;
    }
  
  for (i=0; i<Nz; i++)
    { 
      for (j=0; j<Nr; j++)
	if (j==0)
	  {
	    a[j]=0;
	    b[j]=-1-backgroundDensityProfile(s, j*dr)*dr*dr/4.;
	    c[j]=1;
	    r[j]=(density[i][j]-backgroundDensityProfile(s, j*dr)*(1+a2_half[i][j]/2.))*dr*dr/4.;
	  }
	else
	  {
	    a[j]=1.-1./(2.*j);
	    b[j]=-2-backgroundDensityProfile(s, j*dr)*dr*dr;
	    c[j]=1.+1./(2.*j);
	    r[j]=(density[i][j]-backgroundDensityProfile(s, j*dr)*(1+a2_half[i][j]/2.))*dr*dr;
	  }
      
      solveTridiag(a, b, c, r, psi[i], Nr);
    }
}

void swipeParticle(int iter_max, double error_max)
{
  /*
    - fluid/linear calculation: if the transverse density profile is rapidly varying the method is not consistent
    - iter_max and error_max are not used in this case but are kept for homogeneity
  */
  int i, j;
  static int FIRST_TIME=1;
  static double *source;

  if (FIRST_TIME)
    {
      printf("this is:: swipeParticle FLUID LINEAR [computeAllTheFields=%i]\n", computeAllTheFields);

      source=(double*)malloc(Nr*sizeof(double));
      FIRST_TIME=0;
    }

  solve_delta(); 
  solve_psi_noAnal();   

  // compute only Ez, Er_Bphi and proper density for driver(s) evolutions
  for (j=0; j<Nr-1; j++)
    {
      Ez[0][j]=-(-3*psi[0][j]+4*psi[1][j]-psi[2][j])/(2*dz);
      Er_Bphi[0][j]=-(psi[0][j+1]-psi[0][abs(j-1)])/(2*dr);
    }
  for (i=1; i<Nz; i++)
    for (j=0; j<Nr; j++)
      {
	Ez[i][j]=-(psi[i+1][j]-psi[i-1][j])/(2*dz);
	Er_Bphi[i][j]=-(psi[i][j+1]-psi[i][abs(j-1)])/(2*dr);
      }

  for (i=0; i<Nz; i++)
    for (j=0; j<Nr; j++)
      rho[i][j]=density[i][j]*(1-a2_half[i][j]/2.);
  
  if (computeAllTheFields==NO) return;
  
  for (i=0; i<Nz; i++)
    for (j=0; j<Nr; j++)
      Jz[i][j]=backgroundDensityProfile(s, j*dr)*(psi[i][j]-a2_half[i][j]/2)+rhob[i][j]+Jz_ionization[i][j];
  
  for (i=0; i<Nz; i++)
    {
      if (i==0)					
	for (j=0; j<Nr; j++)
	  source[j]=Jz[0][j]-(-3*Ez[0][j]+4*Ez[1][j]-Ez[2][j])/(2*dz);
      else
	for (j=0; j<Nr; j++)
	  source[j]=Jz[i][j]-(Ez[i+1][j]-Ez[i-1][j])/(2*dz);

      solveBRecurrence(Bphi[i], source);
    }

  for (j=0; j<Nr; j++)
    Jr[0][j]=(-3*Er_Bphi[0][j]+4*Er_Bphi[1][j]-Er_Bphi[2][j])/(2*dz);
  for (i=1; i<Nz; i++)
    for (j=0; j<Nr; j++)
      Jr[i][j]=(Er_Bphi[i+1][j]-Er_Bphi[i-1][j])/(2*dz);
}
#endif

#ifdef ACCURATE
  
void solve_psi_noAnal(int i)
{
  int j, k;
  static int FIRST_TIME=1;
  static double *a, *b, *c, *r;

  /*
    first solve for density perturbation (density)
  */
  
  if (FIRST_TIME)
    {
      a=(double*)malloc(Nr*sizeof(double));
      b=(double*)malloc(Nr*sizeof(double));
      c=(double*)malloc(Nr*sizeof(double));
      r=(double*)malloc(Nr*sizeof(double));

      FIRST_TIME=0;
    }
  
  for (j=0; j<Nr; j++)
    if (j==0)
      {
	a[j]=0;
	b[j]=-1-backgroundDensityProfile(s, j*dr)*dr*dr/4.;
	c[j]=1;
	r[j]=(density[i][j]-backgroundDensityProfile(s, j*dr)*(1+a2_half[i][j]/2.))*dr*dr/4.;
      }
    else
      {
	a[j]=1.-1./(2.*j);
	b[j]=-2-backgroundDensityProfile(s, j*dr)*dr*dr;
	c[j]=1.+1./(2.*j);
	r[j]=(density[i][j]-backgroundDensityProfile(s, j*dr)*(1+a2_half[i][j]/2.))*dr*dr;
      }
  
  solveTridiag(a, b, c, r, psi[i], Nr);
}

void swipeParticle(int iter_max, double error_max)
{
  /*
    - fluid/linear calculation: if the transverse density profile is rapidly varying the method is not consistent
    - iter_max and error_max are not used in this case but are kept for homogeneity
  */
  int i, j;
  double buffer;
  static int FIRST_TIME=1;
  static double *source;

  if (FIRST_TIME)
    {
      printf("this is:: swipeParticle FLUID LINEAR ACCURATE [computeAllTheFields=%i]\n", computeAllTheFields);
      source=(double*)malloc(Nr*sizeof(double));
      FIRST_TIME=0;
    }

  // nabla(a^2) = nabla(|a|^2/2) --> stored in B2
  for (i=1; i<Nz-1; i++)
    {
      B2[i][0]=4*(a2_half[i][1]-a2_half[i][0])/(dr*dr)+(a2_half[i-1][0]-2*a2_half[i][0]+a2_half[i+1][0])/(dz*dz);
      for (j=1; j<Nr-1; j++)
	B2[i][j]=(a2_half[i][j+1]-2*a2_half[i][j]+a2_half[i][j-1])/(dr*dr)+(a2_half[i][j+1]-a2_half[i][j-1])/(2*j*dr*dr)+(a2_half[i+1][j]-2*a2_half[i][j]+a2_half[i-1][j])/(dz*dz);
    }
  
  for (j=0; j<Nr; j++)
    density[Nz-1][j]=density[Nz][j]=density[Nz+1][j]=backgroundDensityProfile(s, j*dr);
  
  for (i=Nz-1; i>0; i--)
    {
      solve_psi_noAnal(i);

      for (j=0; j<Nr; j++)
	{
	  Ez[i][j]=-(-3*psi[i][j]+4*psi[i+1][j]-psi[i+2][j])/(2*dz);
	  Er_Bphi[i][j]=-(psi[i][MIN(j+1, Nr-1)]-psi[i][abs(j-1)])/(2*dr);
	}

      // compute dEz/dxi --> stored in uz
      for (j=0; j<Nr; j++)
	uz[i][j]=(-3*Ez[i][j]+4*Ez[i+1][j]-Ez[i+2][j])/(2*dz);

      // compute Er --> sored in ur
      for (j=0; j<Nr; j++)
	source[j]=-density[i][j]+backgroundDensityProfile(s, j*dr)+rhob[i][j]+Jz_ionization[i][j]-uz[i][j];
      solveBRecurrence(ur[i], source);
      
      // density
      for (j=0; j<Nr; j++)
	{
	  buffer=(backgroundDensityProfile(s, (j+1)*dr)-backgroundDensityProfile(s, fabs((j-1)*dr)))/(2*dr);
	  density[i-1][j]=2*density[i][j]-density[i+1][j]+dz*dz*(backgroundDensityProfile(s, j*dr)*(backgroundDensityProfile(s, j*dr)+rhob[i][j]+Jz_ionization[i][j]+B2[i][j]/2.-density[i][j])+buffer*(ur[i][j]+dr_a2_half[i][j]/2.));
	}
    }

  solve_psi_noAnal(0);
  
  for (j=0; j<Nr; j++)
    {
      Ez[0][j]=-(-3*psi[0][j]+4*psi[1][j]-psi[2][j])/(2*dz);
      Er_Bphi[0][j]=-(psi[0][MIN(j+1, Nr-1)]-psi[0][abs(j-1)])/(2*dr);
    }
  
  for (i=0; i<Nz; i++)
    for (j=0; j<Nr; j++)
      rho[i][j]=density[i][j]*(1-a2_half[i][j]/2.);

  if (computeAllTheFields==NO) return;
  
  for (i=0; i<Nz; i++)
    for (j=0; j<Nr; j++)
      Jz[i][j]=backgroundDensityProfile(s, j*dr)*(psi[i][j]-a2_half[i][j]/2)+rhob[i][j]+Jz_ionization[i][j];
  
  for (i=0; i<Nz; i++)
    {
      if (i==0)					
	for (j=0; j<Nr; j++)
	  source[j]=Jz[0][j]-(-3*Ez[0][j]+4*Ez[1][j]-Ez[2][j])/(2*dz);
      else
	for (j=0; j<Nr; j++)
	  source[j]=Jz[i][j]-(Ez[i+1][j]-Ez[i-1][j])/(2*dz);

      solveBRecurrence(Bphi[i], source);
    }

  for (j=0; j<Nr; j++)
    Jr[0][j]=(-3*Er_Bphi[0][j]+4*Er_Bphi[1][j]-Er_Bphi[2][j])/(2*dz);
  for (i=1; i<Nz; i++)
    for (j=0; j<Nr; j++)
      Jr[i][j]=(Er_Bphi[i+1][j]-Er_Bphi[i-1][j])/(2*dz);
}
#endif

#endif
#endif

void indexOfRefractionNoWake(int nonlinearIndexOfRefraction)
{
  /*
    If nonlinearIndexOfRefraction = YES we consider the contribution from |a|
  */

  int i, j;
  
  for (i=0; i<Nz; i++)
    for (j=0; j<Nr; j++)
      {
	if (nonlinearIndexOfRefraction==NO)
	  rho[i][j]=backgroundDensityProfile(s, j*dr);
	else
	  rho[i][j]=backgroundDensityProfile(s, j*dr)*(1-a2_half[i][j]/2.);
      }
}

void solveTridiag_complex(complex *a, complex *b, complex *c, complex *r, complex *y, int N)
{
  int i;
  complex buffer;

  for (i=1; i<N; i++)
    {
      buffer=a[i]/b[i-1];
      
      b[i]-=buffer*c[i-1];
      r[i]-=buffer*r[i-1];
    }

  y[N-1]=r[N-1]/b[N-1];
  for (i=N-2; i>=0; i--)
    y[i]=(r[i]-c[i]*y[i+1])/b[i];
  
  
}

void evolveLaserPulse()
{
  // COMPLETISSIMO: trucco per catturare la fase + derivata secomda temporale
  // CN2 per equazione completa per envelope (con derivata seconda). La parte con la corrente e' impliticizzata. 
  // usata media (0.5,0,0.5)

  /*
    Mettere l'evoluzione DOPO i calcoli laser plasma!
    Questa funzione, dati in ingresso i campi laser a_old=a(s-ds), a=a(s) e rho=rho(s),
    restituisce il campo al tempo s+ds e rinomina i campi nel seguente modo: a_old=a(s), a=a(s+ds)
  */

  int i, j;
  double d_theta1, d_theta2, d_theta;
  static double *phase;
  static complex *A, *B, *C, *r, *b_1, *b_2;
  static int FIRST_TIME=1;
  complex a_m, a_0, a_p;
  // double S0, S2, S4, A0, A2;

  if (self_consistent_laser==YES && FIRST_TIME==NO && RESET_evolveLaserPulse==YES)
    {
      FIRST_TIME=1;
      printf("evolveLaserPulse [phaseTracking].. resetting the laser solver\n");
      return;
    }
  
  if (self_consistent_laser==NO)
    {
      // assuming the propagation distance "s" is properly set
      // a_old and tmp_a are not set
     
      printf("this is:: evolveLaserPulse [prescribed_laser] s=%f\n", s);
 
      for (i=0; i<Nz; i++)
	{
	  for (j=0; j<Nr; j++)
	    {
	      a[i][j]=laserEnvelope(zmin+i*dz, j*dr, s);  
	      a2_half[i][j]=0.5*a[i][j]*conj(a[i][j]);
	    }
	  
	  for (j=0; j<Nr-1; j++)
	    dr_a2_half[i][j]=(a2_half[i][j+1]-a2_half[i][abs(j-1)])/(2*dr);
	  dr_a2_half[i][Nr-1]=0;
	}

      return;
    }

  if (filterProperDensity==YES)
    for (i=0; i<Nz; i++)
      filter4xB_C(rho[i], Nr, Npass_filterProperDensity);

  if (FIRST_TIME)
    {
      printf("this is:: evolveLaserPulse [phaseTracking]\n");

      A=(complex*)malloc(Nr*sizeof(complex));
      B=(complex*)malloc(Nr*sizeof(complex));
      C=(complex*)malloc(Nr*sizeof(complex));
      r=(complex*)malloc(Nr*sizeof(complex));
 
      b_1=(complex*)malloc(Nr*sizeof(complex));
      b_2=(complex*)malloc(Nr*sizeof(complex));
      
      phase=(double*)malloc((Nz+NZ_GHOST)*sizeof(double));
           
      for (j=0; j<Nr; j++)
	if (j==0)
	  {
	    A[0]  =0.;
	    b_1[0]=-1/ds+I*k0_kp-_NONPARAXIAL_TERMS_*1.5/dz-2*ds/(dr*dr);
	    b_2[0]= 1/ds+I*k0_kp-_NONPARAXIAL_TERMS_*1.5/dz+2*ds/(dr*dr);
	    C[0]  =2*ds/(dr*dr);
	  }
	else 
	  {
	    A[j]  =0.5*ds*(1.-0.5/(double)j)/(dr*dr);
	    b_1[j]=-1/ds+I*k0_kp-_NONPARAXIAL_TERMS_*1.5/dz-ds/(dr*dr);
	    b_2[j]= 1/ds+I*k0_kp-_NONPARAXIAL_TERMS_*1.5/dz+ds/(dr*dr);
	    C[j]  =0.5*ds*(1.+0.5/(double)j)/(dr*dr);
	  }

      C[Nr-1]=0.;
      
      FIRST_TIME=0;
    }

  for (i=0; i<Nz+NZ_GHOST; i++)
    phase[i]=atan2(cimag(a[i][0]), creal(a[i][0]));

  for (i=Nz-1; i>=0; i--)
    {
      // if (skipPoint(i)) continue;
      
      d_theta1=phase[i+1]-phase[i  ];
      d_theta2=phase[i+2]-phase[i+1];

      // d_theta>0 e vicino al bordo superiore
      if (d_theta1<-1.5*M_PI)
	d_theta1+=2*M_PI;

      if (d_theta2<-1.5*M_PI)
	d_theta2+=2*M_PI;

      // d_theta<0 e vicino al bordo inferiore
      if (d_theta1>1.5*M_PI)
	d_theta1-=2*M_PI;

      if (d_theta2>1.5*M_PI)
	d_theta2-=2*M_PI;

      d_theta=1.5*d_theta1/dz-0.5*d_theta2/dz;

      for (j=0; j<Nr; j++)
	{
	   B[j]=b_1[j]-0.5*ds*rho[i][j]+_NONPARAXIAL_TERMS_*I*d_theta;
	  
	  r[j]=-2*a[i][j]/ds
	    -A[j]*a_old[i][abs(j-1)]+(b_2[j]+0.5*ds*rho[i][j]+_NONPARAXIAL_TERMS_*I*d_theta)*a_old[i][j]-C[j]*a_old[i][MIN(j+1, Nr-1)] // A[0]=0, C[Nr-1]=0 
	    -_NONPARAXIAL_TERMS_*4*(tmp_a[i+1][j]-a_old[i+1][j])*cexp(I*(phase[i]-phase[i+1]))/(2*dz)
	    +_NONPARAXIAL_TERMS_*  (tmp_a[i+2][j]-a_old[i+2][j])*cexp(I*(phase[i]-phase[i+2]))/(2*dz);
	}

      solveTridiag_complex(A, B, C, r, tmp_a[i], Nr);

      tmp_a[i][Nr-1]=0;
    }
  
  for (i=0; i<Nz; i++)
    {
      //memcpy((void*)a_old[i], a[i], Nr*sizeof(complex));
      //memcpy((void*)a[i], tmp_a[i], Nr*sizeof(complex));
      for (j=0; j<Nr; j++)
	{
	  a_m=a_old[i][j];
	  a_0=a[i][j];
	  a_p=tmp_a[i][j];
	  
	  a_old[i][j]=a_0;
	  //a[i][j]=a_p; // no interpolazioni a t=k+1 (tridiagonale puro)
	  a[i][j]=alpha_A*a_p+(1-alpha_A)*(2*a_0-a_m); // estrapolo k+1 da k e k-1 e poi combino linearmente con quello ottenuto dal tridiagonale 
	}

      for (j=0; j<Nr; j++)
	a2_half[i][j]=0.5*creal(a[i][j]*conj(a[i][j]));

      if (filterLaserIntensity==YES)
	filter4xB_C(a2_half[i], Nr, Npass_filterLaserIntensity);
      
      for (j=0; j<Nr-1; j++)
	dr_a2_half[i][j]=(a2_half[i][j+1]-a2_half[i][abs(j-1)])/(2*dr);
      dr_a2_half[i][Nr-1]=0;
    }
}

double uniformDensityProfile(double z, double r)
{
  return 1.;
}

complex gaussian_envelope(double z, double r, double s)
{
  double zz=s-kpzf;
  double w0=kpW;
  double w;
  double Zr=0.5*k0_kp*w0*w0;

  w=w0*sqrt(1+zz*zz/(Zr*Zr));
  return a0*(w0/w)*exp(-r*r/(w*w))*cexp(-I*atan(zz/Zr))*cexp(I*(zz/Zr)*r*r/(w*w))*exp(-pow((z-kpz0-delta_g*s)/kpL,2));
}

complex gaussian_envelope_simplified(double z, double r, double s)
{
  return a0*exp(-r*r/(kpW*kpW))*exp(-pow((z-kpz0-delta_g*s)/kpL,2));
}

complex supergaussian_envelope_simplified(double z, double r, double s)
{
  /*
    n_supergaussian=1: Gaussian
    n_supergaussian=2,3,4,....: super-Gaussian    
  */

  return a0*exp(-pow(r/kpW, 2*n_supergaussian))*exp(-pow((z-kpz0-delta_g*s)/kpL,2));
}

void dumpFieldASCII(double **field, char *fileName)
{
  int i, j;
  FILE *f;
  char buffer[200];
  
  sprintf(buffer, "%s/%s_%.1f", DIRECTORY_OUTPUT, fileName, s);
  f=fopen(buffer, "w");
  fprintf(f, "%i %i 1\n", Nz, Nr);
  fprintf(f, "%e 0 %f %f\n", zmin, zmax, rmax);
  
  for (i=0; i<Nz; i++)
    for (j=0; j<Nr; j++)
      fprintf(f, "%e %e %e\n", zmin+i*dz, j*dr, field[i][j]);
  fclose(f);
}

void dumpBeamASCII(char *fileName, int indx0, int Np)
{
  int i, indx0_, Np_;
  FILE *f;
  char buffer[200];
  
  if (beamDriver==NO)
    return;

  if (indx0==ALL_PARTICLES && Np==ALL_PARTICLES)
    {
      indx0_=0;
      Np_=Npart;
    }
  else
    {
      indx0_=indx0;
      Np_=Np;
    }

  sprintf(buffer, "%s/%s_%.1f", DIRECTORY_OUTPUT, fileName, s);
  f=fopen(buffer, "w");
  
  for (i=0; i<Np_; i++)
    fprintf(f, "%e %e %e %e %e %e %e %e %i %i\n", xpb[indx0_+i], ypb[indx0_+i], zpb[indx0_+i], uxb[indx0_+i], uyb[indx0_+i], uzb[indx0_+i], q[indx0_+i], me_mb[indx0_+i], particle_beam_active[indx0_+i], self_consistent[indx0_+i]);
  
  fclose(f);
}

void binaryDump_beam(char *fileName, int indx0, int Np)
{
  int i, ibuf[2], indx0_, Np_;
  double buf[8];
  FILE *f;
  char buffer[200];
  
  if (beamDriver==NO)
    return;

  if (indx0==ALL_PARTICLES && Np==ALL_PARTICLES)
    {
      indx0_=0;
      Np_=Npart;
    }
  else
    {
      indx0_=indx0;
      Np_=Np;
    }
  
  sprintf(buffer, "%s/%s_%.1f", DIRECTORY_OUTPUT, fileName, s);
  f=fopen(buffer, "w");
  
  for (i=0; i<Np_; i++)
    {
      buf[0]=xpb[indx0_+i];
      buf[1]=ypb[indx0_+i];
      buf[2]=zpb[indx0_+i];
      buf[3]=uxb[indx0_+i];
      buf[4]=uyb[indx0_+i];
      buf[5]=uzb[indx0_+i];
      buf[6]=q[indx0_+i];
      buf[7]=me_mb[indx0_+i];
      fwrite(buf, sizeof(double), 8, f);
      
      ibuf[0]=particle_beam_active[indx0_+i];
      ibuf[1]=self_consistent[indx0_+i];
      fwrite(ibuf, sizeof(int), 2, f);
    }
  
  fclose(f);
}

int dumpNow(double time, double dtime, double *time_dump, int N_time_dump)
{
  int i;

  for (i=0; i<N_time_dump; i++)
    if (fabs(time_dump[i]-time)<dtime/2)
      break;

  if (i<N_time_dump)
    return 1;
  else
    return 0;
}

void dumpLaserEnvelope(char *fileName)
{
  int i, j;
  FILE *f;
  char buffer[200];
  
  sprintf(buffer, "%s/%s_%.1f", DIRECTORY_OUTPUT, fileName, s);
  f=fopen(buffer, "w");
  fprintf(f, "%i %i 1\n", Nz, Nr);
  fprintf(f, "%e 0 %f %f\n", zmin, zmax, rmax);
  
  for (i=0; i<Nz; i++)
    for (j=0; j<Nr; j++)
      fprintf(f, "%e %e %e\n", zmin+i*dz, j*dr, sqrt(2*a2_half[i][j]));
  fclose(f);
}


void dumpLineoutASCII(double **field, double r, char *fileName)
{
  int i, j;
  FILE *f;
  char buffer[200];
  
  sprintf(buffer, "%s/%s_%.1f", DIRECTORY_OUTPUT, fileName, s);
  f=fopen(buffer, "w");
  
  for (i=0; i<Nz; i++)
    fprintf(f, "%e %e\n", zmin+i*dz, field[i][(int)(r/dr+0.5)]);
  
  fclose(f);
}

double a0_max()
{
  int i, j;
  double a2_half_max=-1;
  
  for (i=0; i<Nz; i++)
    for (j=0; j<Nr; j++)
      a2_half_max=MAX(a2_half[i][j], a2_half_max);
  
  return sqrt(2*a2_half_max);
}

double intensityCentroid()
{
  int i, j;
  double w, z, zc1, norm;
  
  zc1=norm=0;
  for (i=0; i<Nz; i++)
    for (j=0; j<Nr; j++)
      {
	z=zmin+i*dz;

	w=a2_half[i][j]*j;	
	norm+=w;
	zc1+=z*w;
      }
  zc1/=norm;
  
  return zc1;
}

double laserEnergy(int normalizeEnergy)
{
  int i, j;
  double energy;
  complex buf1, buf2;
  static int FIRST_TIME=YES;
  static double energy0;
  
  energy=0;
  for (i=1; i<Nz-1; i++)
    for (j=0; j<Nr-1; j++)
      {
	buf1=k0_kp*a[i][j]-I*(a[i+1][j]-a[i-1][j])/(2*dz);
	buf2=(a[i][j+1]-a[i][abs(j-1)])/(2*dr);
	
	energy+=j*creal(buf1*conj(buf1)+0.5*buf2*conj(buf2));
      }

  if (FIRST_TIME==YES)
    {
      if (normalizeEnergy==YES)
	energy0=energy;
      else
	energy0=1;
      
      FIRST_TIME=NO;
    }

  return energy/energy0;
}

void binaryDump_field(char *fileName)
{
  /*
    Same structure as in INF&RNO for (density, Ez, Er, Bphi, Re[a], Im[a]) + (Re[a_old], Im[a_old], Jz, Jr) = 6 + 4 = 10 data fields
  */

  int i, j, ibuf[3];
  double dbuf[4], *buffer_data=NULL;
  FILE *f;
  char buffer[200];
  
  sprintf(buffer, "%s/%s_%.1f", DIRECTORY_OUTPUT, fileName, s);
  f=fopen(buffer, "w");
  
  ibuf[0]=Nz;
  ibuf[1]=Nr;
  ibuf[2]=10;
  fwrite(ibuf, sizeof(int), 3, f);

  dbuf[0]=zmin;
  dbuf[1]=0;
  dbuf[2]=zmax;
  dbuf[3]=rmax;
  fwrite(dbuf, sizeof(double), 4, f);

  if (buffer_data==NULL)
    buffer_data=(double*)malloc(Nr*sizeof(double));
  
  for (i=0; i<Nz; i++)
    {
      fwrite(density[i], sizeof(double), Nr, f);
      fwrite(Ez[i], sizeof(double), Nr, f);
      
      for (j=0; j<Nr; j++)
	buffer_data[j]=Er_Bphi[i][j]+Bphi[i][j];
      fwrite(buffer_data, sizeof(double), Nr, f);
      
      fwrite(Bphi[i], sizeof(double), Nr, f);

      for (j=0; j<Nr; j++)
	buffer_data[j]=creal(a[i][j]);
      fwrite(buffer_data, sizeof(double), Nr, f);
      
      for (j=0; j<Nr; j++)
	buffer_data[j]=cimag(a[i][j]);
      fwrite(buffer_data, sizeof(double), Nr, f);

      for (j=0; j<Nr; j++)
	buffer_data[j]=creal(a_old[i][j]);
      fwrite(buffer_data, sizeof(double), Nr, f);
      
      for (j=0; j<Nr; j++)
	buffer_data[j]=cimag(a_old[i][j]);
      fwrite(buffer_data, sizeof(double), Nr, f);
      
      fwrite(Jz[i], sizeof(double), Nr, f);
      fwrite(Jr[i], sizeof(double), Nr, f);
    }
  
  fclose(f);
  free(buffer_data);
}

void beamParticleAllocate()
{
  self_consistent=(int*)realloc(self_consistent, Npart*sizeof(int));
  particle_beam_active=(int*)realloc(particle_beam_active, Npart*sizeof(int));
    
  xpb=(double*)realloc(xpb, Npart*sizeof(double));
  ypb=(double*)realloc(ypb, Npart*sizeof(double));
  zpb=(double*)realloc(zpb, Npart*sizeof(double));
  uxb=(double*)realloc(uxb, Npart*sizeof(double));
  uyb=(double*)realloc(uyb, Npart*sizeof(double));
  uzb=(double*)realloc(uzb, Npart*sizeof(double));

  q=(double*)realloc(q, Npart*sizeof(double));
  me_mb=(double*)realloc(me_mb, Npart*sizeof(double));
}

void depositBeamCharge_linearShapeFunction()
{
  /*
    - linear shape function
    - no check of the boundaries!
  */
  int i, j, k;
  double r, buf_z, buf_r, w00, w01, w10, w11, buffer1, buffer2;
  
  if (beamDriver==NO)
    return;

  for(i=0; i<Nz; i++)
    memset(rhob[i], 0, Nr*sizeof(double));
    
  for (k=0; k<Npart; k++)
    {
      if (particle_beam_active[k]==NO || self_consistent[k]==NO) continue;
      
      buf_z=(zpb[k]-zmin)/dz;
      i=(int)buf_z;

      r=hypot(xpb[k], ypb[k]);      
      buf_r=r/dr;
      j=(int)buf_r;
          
      buf_z-=(double)i;
      buf_r-=(double)j;
      
      w00=(1-buf_z)*(1-buf_r);
      w01=(1-buf_z)*buf_r;
      w10=buf_z*(1-buf_r);
      w11=buf_z*buf_r;

      buffer1=q[k]/(j==0?(1./6.):j); 
      buffer2=q[k]/(j+1.);	  
      
      rhob[i  ][j  ]+=w00*buffer1;
      rhob[i  ][j+1]+=w01*buffer2;
      rhob[i+1][j  ]+=w10*buffer1;
      rhob[i+1][j+1]+=w11*buffer2;
    }
  
  for (i=0; i<Nz; i++)
    for (j=0; j<Nr; j++)
      rhob[i][j]/=(2*M_PI*dz*dr*dr);
}

void depositBeamCharge_quadraticShapeFunction()
{
  /*
    - quadratic shape function
    - no check of the boundaries!
  */
  int i, j, i1, j1, ii, jj, k;
  double r, buf, buf2, w_z[3], w_r[3], r0, dvol;
  
  if (beamDriver==NO)
    return;

  for(i=0; i<Nz; i++)
    memset(rhob[i], 0, Nr*sizeof(double));
    
  for (k=0; k<Npart; k++)
    {
      if (particle_beam_active[k]==NO || self_consistent[k]==NO) continue;
      
      buf=(zpb[k]-zmin)/dz;
      i=(int)(buf+0.5);
      buf-=(double)i;
      buf2=buf*buf;
      w_z[1]=0.75-buf2;
      w_z[2]=0.5*(0.25+buf+buf2);
      w_z[0]=1.-w_z[1]-w_z[2];
	
      r=hypot(xpb[k], ypb[k]);      
      buf=r/dr;
      j=(int)(buf+0.5);
      buf-=(double)j;
      buf2=buf*buf;
      w_r[1]=0.75-buf2;
      w_r[2]=0.5*(0.25+buf+buf2);
      w_r[0]=1.-w_r[1]-w_r[2];
          
      for (i1=0; i1<3; i1++)
	for (j1=0; j1<3; j1++)
	  {
	    ii=i+i1-1;
	    jj=j+j1-1;
	    
	    rhob[ii][abs(jj)]+=q[k]*w_z[i1]*w_r[j1];
	  }
    }

  for (i=0; i<Nz; i++)
    for (j=0; j<Nr; j++)
      {
	r0=j*dr;
	
	if (j==0)
	  dvol=VOL0_QUADRATIC*2*M_PI*dr*dr*dz;
	else if (j==1)
	  dvol=VOL1_QUADRATIC*2*M_PI*r0*dr*dz;
	else
	  dvol=2*M_PI*r0*dr*dz;
	
	rhob[i][j]/=dvol;
      }
}

void pushParticles_linearShapeFunction()
{
  int i, j, k;
  double r, buf_z, buf_r, w00, w01, w10, w11;
  double xp_, yp_, zp_, Fx, Fy, Fz, phi, gamma_inv, qq;
  double int_Ez, int_Er_Bphi;

  if (beamDriver==NO)
    return;

  for (k=0; k<Npart; k++)
    {
      if (particle_beam_active[k]==NO) continue;

      gamma_inv=1./sqrt(1+uxb[k]*uxb[k]+uyb[k]*uyb[k]+uzb[k]*uzb[k]);
      xp_=xpb[k]+ uxb[k]*gamma_inv   *(ds/2.);
      yp_=ypb[k]+ uyb[k]*gamma_inv   *(ds/2.);
      zp_=zpb[k]+(uzb[k]*gamma_inv-1)*(ds/2.);
      
      if (zp_<=zmin+2*dz || hypot(xp_, yp_)>=rmax-2*dr)
	{
	  particle_beam_active[k]=NO;
	  continue;
	}

      buf_z=(zp_-zmin)/dz;
      i=(int)buf_z;

      r=hypot(xp_, yp_);      
      buf_r=r/dr;
      j=(int)buf_r;
          
      buf_z-=(double)i;
      buf_r-=(double)j;

      w00=(1-buf_z)*(1-buf_r);
      w01=(1-buf_z)*buf_r;
      w10=buf_z*(1-buf_r);
      w11=buf_z*buf_r;
      
      int_Ez=
	w00*Ez[i  ][j  ]+
	w01*Ez[i  ][j+1]+
	w10*Ez[i+1][j  ]+
	w11*Ez[i+1][j+1];
      
      int_Er_Bphi=
	w00*Er_Bphi[i  ][j  ]+
	w01*Er_Bphi[i  ][j+1]+
	w10*Er_Bphi[i+1][j  ]+
	w11*Er_Bphi[i+1][j+1];

      phi=atan2(yp_, xp_);

      qq=SGN(q[k]);
      Fx=qq*me_mb[k]*int_Er_Bphi*cos(phi);
      Fy=qq*me_mb[k]*int_Er_Bphi*sin(phi);
      Fz=qq*me_mb[k]*int_Ez;

      uxb[k]+=Fx*ds;
      uyb[k]+=Fy*ds;
      uzb[k]+=Fz*ds;

      gamma_inv=1./sqrt(1+uxb[k]*uxb[k]+uyb[k]*uyb[k]+uzb[k]*uzb[k]);
      xpb[k]=xp_+ uxb[k]*gamma_inv   *(ds/2.);
      ypb[k]=yp_+ uyb[k]*gamma_inv   *(ds/2.);
      zpb[k]=zp_+(uzb[k]*gamma_inv-1)*(ds/2.);
    }

  for (k=0; k<Npart; k++)
    {
      if (particle_beam_active[k]==NO) continue;
      
      if (hypot(xpb[k], ypb[k])>=rmax-2*dr || zpb[k]<=zmin+2*dz)
	{
	  particle_beam_active[k]=NO;
	  continue;
	}
    }
}

void pushParticles_quadraticShapeFunction()
{
  int i, j, i1, j1, ii, jj, k;
  double r, buf, buf2, w_z[3], w_r[3];
  double xp_, yp_, zp_, Fx, Fy, Fz, phi, gamma_inv, qq;
  double int_Ez, int_Er_Bphi;

  if (beamDriver==NO)
    return;

  for (k=0; k<Npart; k++)
    {
      if (particle_beam_active[k]==NO) continue;

      gamma_inv=1./sqrt(1+uxb[k]*uxb[k]+uyb[k]*uyb[k]+uzb[k]*uzb[k]);
      xp_=xpb[k]+ uxb[k]*gamma_inv   *(ds/2.);
      yp_=ypb[k]+ uyb[k]*gamma_inv   *(ds/2.);
      zp_=zpb[k]+(uzb[k]*gamma_inv-1)*(ds/2.);
      
      if (zp_<=zmin+3*dz || hypot(xp_, yp_)>=rmax-3*dr)
	{
	  particle_beam_active[k]=NO;
	  continue;
	}
      
      buf=(zp_-zmin)/dz;
      i=(int)(buf+0.5);
      buf-=(double)i;
      buf2=buf*buf;
      w_z[1]=0.75-buf2;
      w_z[2]=0.5*(0.25+buf+buf2);
      w_z[0]=1.-w_z[1]-w_z[2];

      r=hypot(xp_, yp_);      
      buf=r/dr;
      j=(int)(buf+0.5);
      buf-=(double)j;
      buf2=buf*buf;
      w_r[1]=0.75-buf2;
      w_r[2]=0.5*(0.25+buf+buf2);
      w_r[0]=1.-w_r[1]-w_r[2];
      
      int_Ez=int_Er_Bphi=0;
      for (i1=0; i1<3; i1++)
	for (j1=0; j1<3; j1++)
	  {
	    ii=i+i1-1;
	    jj=j+j1-1;

	    int_Ez+=Ez[ii][abs(jj)]*w_z[i1]*w_r[j1];
	    int_Er_Bphi+=Er_Bphi[ii][abs(jj)]*w_z[i1]*w_r[j1]*SGN(jj);
	  }

      phi=atan2(yp_, xp_);

      qq=SGN(q[k]);
      Fx=qq*me_mb[k]*int_Er_Bphi*cos(phi);
      Fy=qq*me_mb[k]*int_Er_Bphi*sin(phi);
      Fz=qq*me_mb[k]*int_Ez;

      uxb[k]+=Fx*ds;
      uyb[k]+=Fy*ds;
      uzb[k]+=Fz*ds;

      gamma_inv=1./sqrt(1+uxb[k]*uxb[k]+uyb[k]*uyb[k]+uzb[k]*uzb[k]);
      xpb[k]=xp_+ uxb[k]*gamma_inv   *(ds/2.);
      ypb[k]=yp_+ uyb[k]*gamma_inv   *(ds/2.);
      zpb[k]=zp_+(uzb[k]*gamma_inv-1)*(ds/2.);
    }

  for (k=0; k<Npart; k++)
    {
      if (particle_beam_active[k]==NO) continue;
      
      if (hypot(xpb[k], ypb[k])>=rmax-3*dr || zpb[k]<=zmin+3*dz)
	{
	  particle_beam_active[k]=NO;
	  continue;
	}
    }
}

int fileType(char *filename)
{
  int j, value;
  char cmd[1000];
  FILE *f;
  
  sprintf(cmd, "file %s > /tmp/fileDataType", filename);
  system(cmd);

  f=fopen("/tmp/fileDataType", "r");
  fgets(cmd, 1000, f);
  fclose(f);

  for (j=0; j<strlen(cmd)-5; j++)
    if (!strncmp(cmd+j, "ASCII", 5))
      return ASCII;
  
  return BINARY;
}

void read_doubleBuffer(double *buffer, int nitems, FILE *f, int type)
{
  int k;

  if (type==BINARY)
    fread(buffer, sizeof(double), nitems, f);
  else if (type==ASCII)
    for (k=0; k<nitems; k++)
      fscanf(f, "%lf", &buffer[k]);
}

void read_intBuffer(int *buffer, int nitems, FILE *f, int type)
{
  int k;

  if (type==BINARY)
    fread(buffer, sizeof(int), nitems, f);
  else if (type==ASCII)
    for (k=0; k<nitems; k++)
      fscanf(f, "%i", &buffer[k]);  
}

void loadLaser(char *filename, double s_start)
{
  /*
    filename with the complete path!
    
    - reads the format written by binaryDump_field (or by the NEW version of binary_dump_allFields in INFERN): a is in position 4-5, a_old is in position 6-7 (counting from 0)
    - Nz, Nr, zmin, zmax, rmax have to be consistent
    - this function overwrites the laser field acording to the values read from the file, and properly set a_old, a, a2_half and dr_a2_half 
    - since the variables "s" and "smin" are reset, the function should be placed right before the temporal cycle to avoid inconsistencies
  */

  int i, j, Nd, ft, ibuf[3];
  double *dbuf; 
  FILE *f;

  s=smin=s_start;

  ft=fileType(filename);  
  f=fopen(filename, "r");
  read_intBuffer(ibuf, 3, f, ft);
  Nd=ibuf[2];
  dbuf=(double*)malloc(MAX(Nd*Nr, 4)*sizeof(double));
  read_doubleBuffer(dbuf, 4, f, ft); 

  for (i=0; i<Nz; i++)
    {
      read_doubleBuffer(dbuf, Nd*Nr, f, ft); 
      for (j=0; j<Nr; j++)
	{
	  a[i][j]    =dbuf[Nr*4+j]+I*dbuf[Nr*5+j];
	  a_old[i][j]=dbuf[Nr*6+j]+I*dbuf[Nr*7+j];
	}
    }
  fclose(f);
  free(dbuf);

  for (i=0; i<Nz; i++)
    for (j=0; j<Nr; j++)
      a2_half[i][j]=0.5*a[i][j]*conj(a[i][j]);

  for (i=0; i<Nz; i++)
    {
      for (j=0; j<Nr-1; j++)
	dr_a2_half[i][j]=(a2_half[i][j+1]-a2_half[i][abs(j-1)])/(2*dr);
      dr_a2_half[i][Nr-1]=0;
    }
}

void loadLaserFromINFERNO(char *filename, double s_start)
{
  /*  
    filename with the complete path!
    
    - reads the BINARY file written by INFERNO: a is in position 4-5 (counting from 0), a_old is not stored (a_old=a is assumed)
    - Nz, Nr, zmin, zmax, rmax have to be consistent
    - this function overwrites the laser field according to the values read from the file, and properly set a_old, a, a2_half and dr_a2_half 
    - since the variables "s" and "smin" are reset, the function should be placed right before the temporal cycle to avoid inconsistencies
  */

  int i, j, Nd, ft, ibuf[3];
  double *dbuf; 
  FILE *f;

  s=smin=s_start;

  ft=fileType(filename);
  f=fopen(filename, "r");
  read_intBuffer(ibuf, 3, f, ft);
  Nd=ibuf[2];
  dbuf=(double*)malloc(MAX(Nd*Nr, 4)*sizeof(double));
  read_doubleBuffer(dbuf, 4, f, ft); 
  
  for (i=0; i<Nz; i++)
    {
      read_doubleBuffer(dbuf, Nd*Nr, f, ft); 
      for (j=0; j<Nr; j++)
	a_old[i][j]=a[i][j]=dbuf[Nr*4+j]+I*dbuf[Nr*5+j];
    }
  fclose(f);
  free(dbuf);
    
  for (i=0; i<Nz; i++)
    for (j=0; j<Nr; j++)
      a2_half[i][j]=0.5*a[i][j]*conj(a[i][j]);

  for (i=0; i<Nz; i++)
    {
      for (j=0; j<Nr-1; j++)
	dr_a2_half[i][j]=(a2_half[i][j+1]-a2_half[i][abs(j-1)])/(2*dr);
      dr_a2_half[i][Nr-1]=0;
    }
}

void loadExtInjBeam(char *filename, double z0, int indx0, int Np)
{
  /*
    filename with the complete path!
    z0: new centroid position (if z0==ORIGINAL_LOCATION the original location is preserved)
    
    - The particles read from file are ADDED to the list of beam particles (if they already exist)
    - If indx0 and Np are specified the number of particles read from file must be Np, otherwise error may occur
    - in case z0 is specified the calculation of the beam center si done taking into account ONLY the active particles (i.e., particle_beam_active[]==YES)
  */

  int i, indx0_, ft, ibuf[2];
  double buf[8], z0_0, norm;
  FILE *f;

  beamDriver=YES;
  
  ft=fileType(filename);

  f=fopen(filename, "r");

  z0_0=norm=0;
  if (indx0==ALL_PARTICLES && Np==ALL_PARTICLES)
    {
      indx0_=Npart;
      while(1)
	{
	  read_doubleBuffer(buf, 8, f, ft); 
	  read_intBuffer(ibuf, 8, f, ft);
	  if (feof(f)) break;
	  Npart++;
	  beamParticleAllocate();
	  i=Npart-1;
	  xpb[i]=buf[0];
	  ypb[i]=buf[1];
	  zpb[i]=buf[2];
	  uxb[i]=buf[3];
	  uyb[i]=buf[4];
	  uzb[i]=buf[5];
	  q[i]=buf[7];
	  me_mb[i]=buf[8];
	  
	  particle_beam_active[i]=ibuf[0];
	  self_consistent[i]=ibuf[1];
	  
	  if (particle_beam_active[i]==YES)
	    {
	      z0_0+=q[i]*zpb[i];
	      norm+=q[i];
	    }
	}
    }
  else
    {
      Npart+=Np;
      beamParticleAllocate();
      
      for (i=0; i<Np; i++)
	fscanf(f, "%lf %lf %lf %lf %lf %lf %lf %lf %i %i\n", &xpb[indx0+i], &ypb[indx0+i], &zpb[indx0+i], &uxb[indx0+i], &uyb[indx0+i], &uzb[indx0+i], &q[indx0+i], &me_mb[indx0+i], &particle_beam_active[indx0+i], &self_consistent[indx0+i]);

      for (i=0; i<Np; i++)
	if (particle_beam_active[indx0+i]==YES)
	  {
	    z0_0+=q[indx0+i]*zpb[indx0+i];
	    norm+=q[indx0+i];
	  }
    }  
  
  fclose(f);
  
  if (z0!=ORIGINAL_LOCATION)
    {
      z0_0/=norm;

      if (indx0==ALL_PARTICLES && Np==ALL_PARTICLES)
	for (i=indx0_; i<Npart; i++)
	  zpb[i]=(zpb[i]-z0_0)+z0;
      else
	for (i=0; i<Np; i++)
	  zpb[indx0+i]=(zpb[indx0+i]-z0_0)+z0;
    }
}

void loadExtInjBeamFromINFERNO(char *filename, double z0, int beamSelfConsistent, double n0)
{
  /*
    - filename = file (with the complete path) from explicit INF&RNO simulation (z, r, uz, ur, q, z0, r0, label) x # macroparticles
    - z0 = new beam centroid position (if z0==ORIGINAL_LOCATION the original location is preserved)
    - beamSelfConsistent = YES, NO
    - n0 = reference plasma density [in e/cm^3]

    N.B.1. the normalization of the quantities must (e.g., z, r, etc..) be consistent between explicit and QS simulation
    N.B.2. the mass of the particles is assumed to be _MASS_ELECTRON_, the charge is assumed to be _CHARGE_ELECTRON_

    This routine transforms the r-z data into x-y-z data (arbitrary azimuthal angle)
    This routine append the beam at the first available position (i.e., at the end)
    Loading the beam from INF&RNO we have: N_physPartBeam = 2*M_PI*dr^2*dz*Qtot_fromSimINF&RNO*(n_0/k_p^3) where dr, dz are normalized to k_p
    The charge in the simulation is normalized according to: Q_sym = (k_p^3/n_0)xN_physPartBeam = 2*M_PI*dr^2*dz*Qtot_fromSimINF&RNO
  */

  int i, ft, indx0;
  double z0_0, norm, buf[8], theta;
  FILE *f;

  beamDriver=YES;
  
  ft=fileType(filename);
  f=fopen(filename, "r");

  z0_0=norm=0;
  indx0=Npart;
  while(1)
    {
      read_doubleBuffer(buf, 8, f, ft); 
      if (feof(f)) break;
      Npart++;
      beamParticleAllocate();
      i=Npart-1;

#ifdef _USE_RAND_
      theta=2.*M_PI*rand()/(double)RAND_MAX;
#else
      theta=2.*M_PI*random()/(double)RAND_MAX;
#endif
      xpb[i]=buf[1]*cos(theta);
      ypb[i]=buf[1]*sin(theta);
      zpb[i]=buf[0];
      uxb[i]=buf[3]*cos(theta);
      uyb[i]=buf[3]*sin(theta);
      uzb[i]=buf[2];
      q[i]=_CHARGE_ELECTRON_*2*M_PI*dr*dr*dz*buf[4];
      me_mb[i]=1.;

      particle_beam_active[i]=YES;
      self_consistent[i]=beamSelfConsistent;

      z0_0+=q[i]*zpb[i];
      norm+=q[i];
    }

  fclose(f);
  
  if (z0!=ORIGINAL_LOCATION)
    {
      z0_0/=norm;
      for (i=indx0; i<Npart; i++)
	zpb[i]=(zpb[i]-z0_0)+z0;
    }
}

void loadBeamFromINFERNO_setCharge(char *filename, double z0, double N_physPartBeam, int beamParticleCharge, int beamSelfConsistent, double n0)
{
  /*
    - filename = file (with the complete path) from explicit INF&RNO simulation (z, r, uz, ur, q, z0, r0, label) x # macroparticles
    - z0 = new beam centroid position (if z0==ORIGINAL_LOCATION the original location is preserved)
    - N_physPartBeam number of physical particles in the bunch
    - beamParticleCharge = _CHARGE_ELECTRON_, _CHARGE_POSITRON_
    - beamSelfConsistent = YES, NO
    - n0 = reference plasma density [in e/cm^3]

    N.B.1. the normalization of the quantities must (e.g., z, r, etc..) be consistent between explicit and QS simulation
    N.B.2. the mass of the particles is assumed to be _MASS_ELECTRON_, the charge has to be specified in beamParticleCharge (background electrons in INF&RNO have a positive charge)

    This routine transforms the r-z data into x-y-z data (arbitrary azimuthal angle)
    This routine append the beam at the first available position (i.e., at the end)
    The charge in the simulation is normalized according to: Q_sym=(k_p^3/n_0)xN_physPartBeam
  */

  int i, ft, indx0;
  double z0_0, norm, buf[8], theta;
  double Qtot=beamParticleCharge*8*M_PI*N_physPartBeam*sqrt(M_PI*n0*R0_cgs*R0_cgs*R0_cgs); // normalized charge of the beam
  FILE *f;

  beamDriver=YES;
  
  ft=fileType(filename);
  f=fopen(filename, "r");

  z0_0=norm=0;
  indx0=Npart;
  while(1)
    {
      read_doubleBuffer(buf, 8, f, ft); 
      if (feof(f)) break;
      Npart++;
      beamParticleAllocate();
      i=Npart-1;

#ifdef _USE_RAND_
      theta=2.*M_PI*rand()/(double)RAND_MAX;
#else
      theta=2.*M_PI*random()/(double)RAND_MAX;
#endif
      xpb[i]=buf[1]*cos(theta);
      ypb[i]=buf[1]*sin(theta);
      zpb[i]=buf[0];
      uxb[i]=buf[3]*cos(theta);
      uyb[i]=buf[3]*sin(theta);
      uzb[i]=buf[2];
      q[i]=buf[4];
      me_mb[i]=1.;

      particle_beam_active[i]=YES;
      self_consistent[i]=beamSelfConsistent;

      z0_0+=q[i]*zpb[i];
      norm+=q[i];
    }

  fclose(f);

  for (i=indx0; i<Npart; i++)
    q[i]=Qtot*q[i]/norm;
  
  if (z0!=ORIGINAL_LOCATION)
    {
      z0_0/=norm;
      for (i=indx0; i<Npart; i++)
	zpb[i]=(zpb[i]-z0_0)+z0;
    }
}

void restoreAold()
{
  // calcolo a_old usando l'approx parassiale
  int i, j;

  for (i=0; i<Nz; i++)
    {

      a_old[i][0]=a[i][0]+0.5*I*ds*(rho[i][0]*a[i][0]-4*(a[i][1]-a[i][0])/(dr*dr))/k0_kp;
      
      for (j=1; j<Nr-1; j++)
	a_old[i][j]=a[i][j]+0.5*I*ds*(rho[i][j]*a[i][j]-(a[i][j+1]-2*a[i][j]+a[i][j-1])/(dr*dr)-(a[i][j+1]-a[i][j-1])/(2*j*dr*dr))/k0_kp;
    }
}

