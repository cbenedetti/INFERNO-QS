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

/*
  - all "metric" quantities are normalized with kp (chosen from the simulation)
  - normalized emittance is normalized with kp
  - dgamma_over_gamma0 is NOT a percentage (e.g., dgamma_over_gamma0=0.1 -> 10%)
  - beamParticleMass = _MASS_ELECTRON_, _MASS_PROTON_ 
  - beamParticleCharge = _CHARGE_ELECTRON_, _CHARGE_POSITRON_
  - beamSelfConsistent = YES, NO
  - indx0 = index of the first particle
  - Np = how many particles for the beam created
  
  Do not update Npart (# of total beam particles)! It's done automatically.
*/

#define _MASS_ELECTRON_ 0.511
#define _MASS_PROTON_ 938.
#define _CHARGE_ELECTRON_ (-1)
#define _CHARGE_POSITRON_ (+1)

void random_gaussian(double *eta1, double *eta2)
{
  double rnd1, rnd2;
  
  rnd1=random()/(double)RAND_MAX;
  rnd2=random()/(double)RAND_MAX;

  *eta1=sqrt(-2*log(rnd1))*cos(2*M_PI*rnd2);
  *eta2=sqrt(-2*log(rnd1))*sin(2*M_PI*rnd2);
}

void generateBeam_GaussianXGaussian(double z0, double sigma_z, double sigma_x, double epsilon_x_norm, double gamma0, double dgamma_over_gamma0, double nb_over_n0, double beamParticleMass, int beamParticleCharge, int beamSelfConsistent, int indx0, int Np)
{
  /*
    Longitudinal Gaussian x Transverse Gaussian
  */
  
  int i, ii;
  double eta1, eta2;
  double sigma_ux=epsilon_x_norm/sigma_x, qq=beamParticleCharge*nb_over_n0*pow(2*M_PI, 1.5)*sigma_z*sigma_x*sigma_x/Np;
  
  beamDriver=YES;
  Npart+=Np;
  beamParticleAllocate();
  
  for (i=0; i<Np; i++)
    {
      ii=i+indx0;

      self_consistent[ii]=beamSelfConsistent;
      particle_beam_active[ii]=YES;
      
      random_gaussian(&eta1, &eta2); 
      xpb[ii]=sigma_x *eta1;
      uxb[ii]=sigma_ux*eta2;
      
      random_gaussian(&eta1, &eta2); 
      ypb[ii]=sigma_x *eta1;
      uyb[ii]=sigma_ux*eta2;
      
      random_gaussian(&eta1, &eta2); 
      zpb[ii]=z0+sigma_z*eta1;
      uzb[ii]=gamma0*(1+dgamma_over_gamma0*eta2);
      
      q[ii]=qq;
      
      me_mb[ii]=_MASS_ELECTRON_/beamParticleMass;
    } 
  
  for (i=0; i<Np; i++)
    {
      ii=i+indx0;

      if (particle_beam_active[ii]==NO) continue;
      
      if (hypot(xpb[ii], ypb[ii])>=rmax-dr || zpb[ii]<=zmin+dz || zpb[ii]>=zmax-dz)
	{
	  particle_beam_active[ii]=NO;
	  continue;
	}
    }
}

void generateBeam_FlattopSymXGaussian(double z0, double Lz, double Rz, double sigma_x, double epsilon_x_norm, double gamma0, double dgamma_over_gamma0, double nb_over_n0, double beamParticleMass, int beamParticleCharge, int beamSelfConsistent, int indx0, int Np)
{
  /*
    Longitudinal Flattop + ramp(symmetric) x Transverse Gaussian
  */

  int i, ii;
  double eta1, eta2, eta3;
  double f1=(Rz/2)/(Lz+Rz), f2=f1+Lz/(Lz+Rz);
  double sigma_ux=epsilon_x_norm/sigma_x, qq=beamParticleCharge*nb_over_n0*(2*M_PI*sigma_x*sigma_x)*(Lz+Rz)/Np;

  beamDriver=YES;
  Npart+=Np; 
  beamParticleAllocate();
  
  for (i=0; i<Np; i++)
    {
      ii=i+indx0;

      self_consistent[ii]=beamSelfConsistent;
      particle_beam_active[ii]=YES;
     
      random_gaussian(&eta1, &eta2); 
      xpb[ii]=sigma_x *eta1;
      uxb[ii]=sigma_ux*eta2;
      
      random_gaussian(&eta1, &eta2); 
      ypb[ii]=sigma_x *eta1;
      uyb[ii]=sigma_ux*eta2;
      
      eta3=random()/(double)RAND_MAX;
      if (eta3<f1)
	{
	  while(1)
	    {
	      eta1=random()/(double)RAND_MAX;
	      eta2=random()/(double)RAND_MAX;
	      if (eta2<pow(sin(eta1*M_PI/2),2))
		break;
	    }
	  
	  zpb[ii]=z0-Lz/2-Rz+Rz*eta1;
	}   
      else if (eta3<f2)
	{
	  eta1=random()/(double)RAND_MAX;
	  zpb[ii]=z0-Lz/2+Lz*eta1;
	}
      else
	{
	  while(1)
	    {
	      eta1=random()/(double)RAND_MAX;
	      eta2=random()/(double)RAND_MAX;
	      if (eta2<pow(cos(eta1*M_PI/2),2))
		break;
	    }
	  zpb[ii]=z0+Lz/2+Rz*eta1;
	}   
	

      random_gaussian(&eta1, &eta2); 
      uzb[ii]=gamma0*(1+dgamma_over_gamma0*eta2);
      
      q[ii]=qq;
      
      me_mb[ii]=_MASS_ELECTRON_/beamParticleMass;
    } 

  for (i=0; i<Np; i++)
    {
      ii=i+indx0;

      if (particle_beam_active[ii]==NO) continue;
      
      if (hypot(xpb[ii], ypb[ii])>=rmax-dr || zpb[ii]<=zmin+dz || zpb[ii]>=zmax-dz)
	{
	  particle_beam_active[ii]=NO;
	  continue;
	}
    }
}

void generateBeam_FlattopAsymXGaussian(double z0, double Lz, double R_front, double R_back, double sigma_x, double epsilon_x_norm, double gamma0, double dgamma_over_gamma0, double nb_over_n0, double beamParticleMass, int beamParticleCharge, int beamSelfConsistent, int indx0, int Np)
{
  /*
    Longitudinal Flattop + ramp(asymmetric) x Transverse Gaussian
  */

  int i, ii;
  double eta1, eta2, eta3;
  double f1=(R_back/2)/(Lz+0.5*(R_front+R_back)), f2=f1+Lz/(Lz+0.5*(R_front+R_back));
  double sigma_ux=epsilon_x_norm/sigma_x, qq=beamParticleCharge*nb_over_n0*(2*M_PI*sigma_x*sigma_x)*(Lz+0.5*(R_front+R_back))/Np;

  beamDriver=YES;
  Npart+=Np;
  beamParticleAllocate();
  
  for (i=0; i<Np; i++)
    {
      ii=i+indx0;

      self_consistent[ii]=beamSelfConsistent;
      particle_beam_active[ii]=YES;
     
      random_gaussian(&eta1, &eta2); 
      xpb[ii]=sigma_x *eta1;
      uxb[ii]=sigma_ux*eta2;
      
      random_gaussian(&eta1, &eta2); 
      ypb[ii]=sigma_x *eta1;
      uyb[ii]=sigma_ux*eta2;
      
      eta3=random()/(double)RAND_MAX;
      if (eta3<f1)
	{
	  while(1)
	    {
	      eta1=random()/(double)RAND_MAX;
	      eta2=random()/(double)RAND_MAX;
	      if (eta2<pow(sin(eta1*M_PI/2),2))
		break;
	    }
	  
	  zpb[ii]=z0-Lz/2-R_back*(1-eta1);
	}   
      else if (eta3<f2)
	{
	  eta1=random()/(double)RAND_MAX;
	  zpb[ii]=z0-Lz/2+Lz*eta1;
	}
      else
	{
	  while(1)
	    {
	      eta1=random()/(double)RAND_MAX;
	      eta2=random()/(double)RAND_MAX;
	      if (eta2<pow(cos(eta1*M_PI/2),2))
		break;
	    }
	  zpb[ii]=z0+Lz/2+R_front*eta1;
	}   
	

      random_gaussian(&eta1, &eta2); 
      uzb[ii]=gamma0*(1+dgamma_over_gamma0*eta2);
      
      q[ii]=qq;
      
      me_mb[ii]=_MASS_ELECTRON_/beamParticleMass;
    } 

  for (i=0; i<Np; i++)
    {
      ii=i+indx0;

      if (particle_beam_active[ii]==NO) continue;
      
      if (hypot(xpb[ii], ypb[ii])>=rmax-dr || zpb[ii]<=zmin+dz || zpb[ii]>=zmax-dz)
	{
	  particle_beam_active[ii]=NO;
	  continue;
	}
    }
}

void generateBeam_TriangularXGaussian(double z0, double L_z, double sigma_x, double epsilon_x_norm, double gamma0, double dgamma_over_gamma0, double nb_over_n0, double beamParticleMass, int beamParticleCharge, int beamSelfConsistent, int indx0, int Np)
{
  /*
    Longitudinal Triangular (peak in the front) x Transverse Gaussian.
    
    z0 = position of the head [i.e., peak of the density] of the beam. 
  */
  
  int i, ii;
  double eta1, eta2;
  double sigma_ux=epsilon_x_norm/sigma_x, qq=beamParticleCharge*nb_over_n0*(2*M_PI*sigma_x*sigma_x)*(L_z/2.)/Np;
  
  beamDriver=YES;
  Npart+=Np;
  beamParticleAllocate();
  
  for (i=0; i<Np; i++)
    {
      ii=i+indx0;

      self_consistent[ii]=beamSelfConsistent;
      particle_beam_active[ii]=YES;
      
      random_gaussian(&eta1, &eta2); 
      xpb[ii]=sigma_x *eta1;
      uxb[ii]=sigma_ux*eta2;
      
      random_gaussian(&eta1, &eta2); 
      ypb[ii]=sigma_x *eta1;
      uyb[ii]=sigma_ux*eta2;
      
      eta1=random()/(double)RAND_MAX;
      zpb[ii]=z0+L_z*(sqrt(eta1)-1);
      random_gaussian(&eta1, &eta2); // eta1 scartato
      uzb[ii]=gamma0*(1+dgamma_over_gamma0*eta2);
      
      q[ii]=qq;
      
      me_mb[ii]=_MASS_ELECTRON_/beamParticleMass;
    } 
  
  for (i=0; i<Np; i++)
    {
      ii=i+indx0;

      if (particle_beam_active[ii]==NO) continue;
      
      if (hypot(xpb[ii], ypb[ii])>=rmax-dr || zpb[ii]<=zmin+dz || zpb[ii]>=zmax-dz)
	{
	  particle_beam_active[ii]=NO;
	  continue;
	}
    }
}

void generateBeam_TrapezoidalXGaussian(double sigma_x, double epsilon_x_norm, double gamma0, double dgamma_over_gamma0, double z_head, double nb_over_n0_head, double z_tail, double nb_over_n0_tail, double beamParticleMass, int beamParticleCharge, int beamSelfConsistent, int indx0, int Np)
{
  /*
    Trapezoidal profile (peak in the front!!! otherwise it does not work) x Transverse Gaussian.    
  */
  
  int i, ii;
  double eta1, eta2, buf, hh;
  double sigma_ux=epsilon_x_norm/sigma_x, qq=beamParticleCharge*(nb_over_n0_head+nb_over_n0_tail)*(2*M_PI*sigma_x*sigma_x)*(z_head-z_tail)/2./Np;

  
  beamDriver=YES;
  Npart+=Np;
  beamParticleAllocate();

  hh=(z_head-z_tail)/(Np-1);
  
  for (i=0; i<Np; i++)
    {
      ii=i+indx0;

      self_consistent[ii]=beamSelfConsistent;
      particle_beam_active[ii]=YES;
      
      random_gaussian(&eta1, &eta2); 
      xpb[ii]=sigma_x *eta1;
      uxb[ii]=sigma_ux*eta2;
      
      random_gaussian(&eta1, &eta2); 
      ypb[ii]=sigma_x *eta1;
      uyb[ii]=sigma_ux*eta2;

      buf=i/(double)(Np-1);
      
      zpb[ii]=z_tail*(1-buf)+z_head*buf;
      random_gaussian(&eta1, &eta2); // eta1 scartato
      uzb[ii]=gamma0*(1+dgamma_over_gamma0*eta2);
      
      q[ii]=beamParticleCharge*(2*M_PI*sigma_x*sigma_x)*(nb_over_n0_tail*(1-buf)+nb_over_n0_head*buf)*hh;
      
      me_mb[ii]=_MASS_ELECTRON_/beamParticleMass;
    } 
  
  for (i=0; i<Np; i++)
    {
      ii=i+indx0;

      if (particle_beam_active[ii]==NO) continue;
      
      if (hypot(xpb[ii], ypb[ii])>=rmax-dr || zpb[ii]<=zmin+dz || zpb[ii]>=zmax-dz)
	{
	  particle_beam_active[ii]=NO;
	  continue;
	}
    }
}

void generateBeam_TriangularXFlattop(double z0, double L_z, double R, double epsilon_x_norm, double gamma0, double dgamma_over_gamma0, double nb_over_n0, double beamParticleMass, int beamParticleCharge, int beamSelfConsistent, int indx0, int Np)
{
  /*
    Longitudinal Triangular (peak in the front) x Transverse uniform (with radius R, in this case sigma_x=R/2).
    
    z0 = position of the head [i.e., peak of the density] of the beam. 
  */
  
  int i, ii;
  double eta1, eta2;
  double sigma_ux=epsilon_x_norm/(0.5*R), qq=beamParticleCharge*nb_over_n0*(M_PI/2)*L_z*R*R/Np;
  
  beamDriver=YES;
  Npart+=Np;
  beamParticleAllocate();
  
  for (i=0; i<Np; i++)
    {
      ii=i+indx0;

      self_consistent[ii]=beamSelfConsistent;
      particle_beam_active[ii]=YES;

      eta1=random()/(double)RAND_MAX;
      eta2=random()/(double)RAND_MAX;
      xpb[ii]=R*sqrt(eta1)*cos(2*M_PI*eta2);
      ypb[ii]=R*sqrt(eta1)*sin(2*M_PI*eta2);
      
      random_gaussian(&eta1, &eta2); 
      uxb[ii]=sigma_ux*eta1;
      uyb[ii]=sigma_ux*eta2;  
      
      eta1=random()/(double)RAND_MAX;
      zpb[ii]=z0+L_z*(sqrt(eta1)-1);
      random_gaussian(&eta1, &eta2); // eta1 scartato
      uzb[ii]=gamma0*(1+dgamma_over_gamma0*eta2);
      
      q[ii]=qq;
      
      me_mb[ii]=_MASS_ELECTRON_/beamParticleMass;
    } 
  
  for (i=0; i<Np; i++)
    {
      ii=i+indx0;

      if (particle_beam_active[ii]==NO) continue;
      
      if (hypot(xpb[ii], ypb[ii])>=rmax-dr || zpb[ii]<=zmin+dz || zpb[ii]>=zmax-dz)
	{
	  particle_beam_active[ii]=NO;
	  continue;
	}
    }
}

void generateBeam_TruncatedtriangularXFlattop(double z0, double L0, double L_z, double R, double epsilon_x_norm, double gamma0, double dgamma_over_gamma0, double nb_over_n0, double beamParticleMass, int beamParticleCharge, int beamSelfConsistent, int indx0, int Np)
{
  /*
    Longitudinal Triangular (peak in the front, L0 would be the total ideal length, L_z is the length [measured from beam head] over which the charge is confined) x Transverse uniform (with radius R, in this case sigma_x=R/2).
    
    z0 = position of the head [i.e., peak of the density] of the beam. 
  */
  
  int i, ii;
  double eta1, eta2;
  double sigma_ux=epsilon_x_norm/(0.5*R), qq=beamParticleCharge*nb_over_n0*(M_PI*R*R)*((2-L_z/L0)*(L_z/2.))/Np;

  if (L_z>L0)
    {
      printf("error [generateBeam_TruncatedtriangularXFlattop]: L_z > L0!! \n\n");
      exit(0);
    }

  beamDriver=YES;
  Npart+=Np;
  beamParticleAllocate();
  
  for (i=0; i<Np; i++)
    {
      ii=i+indx0;

      self_consistent[ii]=beamSelfConsistent;
      particle_beam_active[ii]=YES;

      eta1=random()/(double)RAND_MAX;
      eta2=random()/(double)RAND_MAX;
      xpb[ii]=R*sqrt(eta1)*cos(2*M_PI*eta2);
      ypb[ii]=R*sqrt(eta1)*sin(2*M_PI*eta2);
      
      random_gaussian(&eta1, &eta2); 
      uxb[ii]=sigma_ux*eta1;
      uyb[ii]=sigma_ux*eta2;  
      
      eta1=random()/(double)RAND_MAX;
      zpb[ii]=z0-L0+sqrt(L_z*(2*L0-L_z)*eta1+pow(L0-L_z,2));
      random_gaussian(&eta1, &eta2); // eta1 scartato
      uzb[ii]=gamma0*(1+dgamma_over_gamma0*eta2);
      
      q[ii]=qq;
      
      me_mb[ii]=_MASS_ELECTRON_/beamParticleMass;
    } 
  
  for (i=0; i<Np; i++)
    {
      ii=i+indx0;

      if (particle_beam_active[ii]==NO) continue;
      
      if (hypot(xpb[ii], ypb[ii])>=rmax-dr || zpb[ii]<=zmin+dz || zpb[ii]>=zmax-dz)
	{
	  particle_beam_active[ii]=NO;
	  continue;
	}
    }
}

void generateTestParticleLineBeam(double z_tail, double z_head, double gamma0, int indx0, int Np)
{
  /*
    Line of particles (passive) from z_tail ot z_head and initial energy gamma0 
  */
  
  int i, ii;
  
  beamDriver=YES;
  Npart+=Np;
  beamParticleAllocate();
  
  for (i=0; i<Np; i++)
    {
      ii=i+indx0;

      self_consistent[ii]=NO;
      particle_beam_active[ii]=YES;
      
      xpb[ii]=0;
      uxb[ii]=0;
      
      ypb[ii]=0;
      uyb[ii]=0;
      
      zpb[ii]=z_tail+(z_head-z_tail)*(i/(double)(Np-1));
      uzb[ii]=gamma0;
      
      q[ii]=-1;
      
      me_mb[ii]=1;
    } 
  
  for (i=0; i<Np; i++)
    {
      ii=i+indx0;

      if (particle_beam_active[ii]==NO) continue;
      
      if (hypot(xpb[ii], ypb[ii])>=rmax-dr || zpb[ii]<=zmin+dz || zpb[ii]>=zmax-dz)
	{
	  particle_beam_active[ii]=NO;
	  continue;
	}
    }
}
