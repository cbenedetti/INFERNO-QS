/*
*** Copyright Notice ***

INtegrated Fluid & paRticles simulatioN cOde (INF&RNO) Copyright (c) 2026,
The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy). All rights reserved.

If you have questions about your rights to use or distribute this software,
please contact Berkeley Lab's Intellectual Property Office at
IPO@lbl.gov.

NOTICE.  This Software was developed under funding from the U.S. Department
of Energy and the U.S. Government consequently retains certain rights.  As
such, the U.S. Government has been granted for itself and others acting on
its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the
Software to reproduce, distribute copies to the public, prepare derivative 
works, and perform publicly and display publicly, and to permit others to do so.


****************************

*** License Agreement ***

INtegrated Fluid & paRticles simulatioN cOde (INF&RNO) Copyright (c) 2026,
The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy). All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

(1) Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.

(2) Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

(3) Neither the name of the University of California, Lawrence Berkeley
National Laboratory, U.S. Dept. of Energy nor the names of its contributors
may be used to endorse or promote products derived from this software
without specific prior written permission.


THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

You are under no obligation whatsoever to provide any bug fixes, patches,
or upgrades to the features, functionality or performance of the source
code ("Enhancements") to anyone; however, if you choose to make your
Enhancements available either publicly, or directly to Lawrence Berkeley
National Laboratory, without imposing a separate written license agreement
for such Enhancements, then you hereby grant the following license: a
non-exclusive, royalty-free perpetual license to install, use, modify,
prepare derivative works, incorporate into other computer software,
distribute, and sublicense such enhancements or derivative works thereof,
in binary and source code form.
 */

#include"source_QS_laserIndependentGrid.c"
#include"source_beamInitialDistribution.c"

// user-defined density profile
double matched_radius, z_plasma_start, L_ramp_entrance, L_plasma, L_ramp_exit, z_plasma_end;
double plasma(double z, double r)
{
  double buf, n;
  
  n=1.+4*r*r/pow(matched_radius,4);
  
  if (z<z_plasma_start)
    buf=0;
  else if (z>=z_plasma_start && z<z_plasma_start+L_ramp_entrance)
    buf=(z-z_plasma_start)/L_ramp_entrance;
  else if (z>=z_plasma_start+L_ramp_entrance && z<z_plasma_start+L_ramp_entrance+L_plasma)
    buf=1.;
  else if (z>=z_plasma_start+L_ramp_entrance+L_plasma && z<z_plasma_start+L_ramp_entrance+L_plasma+L_ramp_exit)
    buf=1-(z-(z_plasma_start+L_ramp_entrance+L_plasma))/L_ramp_exit;
  else
    buf=0;

  return n*buf;
}

int main(int narg, char **args)
{
  int iter, k, N_dumpList, Npb, selfConsb;
  double err, dumpList[10000], Zr, ss, s_dumpList;
  double kpz0b, kpsigmazb, kpsigmaxb, kpepsib, nb_over_n0b, gammab, dgammab;
  FILE *f, *g;
  char buffer2[1000];

  // settings of the quasi-static solver
  iter=300;
  err=1e-4;
  alpha_A=1;
  alpha_B=0.98;
  psi_min_0=-0.99;
  d_psi_min=0.01;
  
  // definition of the plasma
  z_plasma_start=6.;
  L_ramp_entrance=100.;
  L_plasma=1000.;
  matched_radius=3.;
  L_ramp_exit=100.;
  z_plasma_end=1300.;
  backgroundDensityProfile=plasma;

  // definition of the laser parameters
  self_consistent_laser=YES;
  laserEnvelope=gaussian_envelope;
  a0=2.;
  kpL=1.5;
  kpW=3.;
  k0_kp=100.;
  kpz0=0.;
  kpzf=0.;
  Zr=0.5*k0_kp*kpW*kpW;

  // definition of the computational box/resolution/numerical parameters
  zmax=6.;
  zmin=-8.;
  dz=1./30.;
  dzl=1./200.;
  compute_Nzl();
  compute_Nz();
  rmax=16.;
  dr=1./20.;
  drl=1./10.;
  compute_Nr();
  compute_Nrl();
  s=smin=0;
  smax=z_plasma_end;
  ds=Zr/300.;
  Nppc=2;
  
  DIRECTORY_OUTPUT="./";
    
  numericalParameters();
  allocateFields();
  allocateParticles();
        
  // definition of the externally injected bunch
  Npb=20000;
  kpz0b=-5.;
  kpsigmazb=0.3;
  kpsigmaxb=0.3;
  kpepsib=0.5;
  nb_over_n0b=1.;
  gammab=2000.;
  dgammab=0.;
  selfConsb=YES;
  generateBeam_GaussianXGaussian(kpz0b, kpsigmazb, kpsigmaxb, kpepsib, gammab, dgammab, nb_over_n0b, _MASS_ELECTRON_, _CHARGE_ELECTRON_, selfConsb, 0, Npb);

  // definition of the output dumps 
  N_dumpList=131;
  s_dumpList=10.;
  for (k=0; k<N_dumpList; k++)
    dumpList[k]=s_dumpList*k;

  // output files
  sprintf(buffer2, "%s/densityProfile", DIRECTORY_OUTPUT);
  f=fopen(buffer2, "w");
  for(ss=smin; ss<=smax; ss+=.1)
    fprintf(f, "%e %e\n", ss, backgroundDensityProfile(ss, 0.));
  fclose(f);
  
  sprintf(buffer2, "%s/a0_max", DIRECTORY_OUTPUT);
  f=fopen(buffer2, "w");

  sprintf(buffer2, "%s/laserEnergy", DIRECTORY_OUTPUT);
  g=fopen(buffer2, "w");

  // temporal loop
  for (s=smin; s<=smax+ds; s+=ds)
    {
      printf("Propagation distance, kp*s=%f\n", s);
      
      depositBeamCharge_linearShapeFunction();
      swipeParticle(iter, err);
      	
      fprintf(f, "%e %e\n", s, a0_max()); fflush(f);
      fprintf(g, "%e %e\n", s, laserEnergy(YES)); fflush(g);
      
      if (dumpNow(s, ds, dumpList, N_dumpList))
	{
	  binaryDump_field("fieldData");
	  dumpLineoutASCII(Ez, 0., "cutEz");
	  dumpBeamASCII("bunch", 0, Npb);
	}

      evolveLaserPulse();
      pushParticles_linearShapeFunction();

    }

  fclose(f);
  fclose(g);
}
