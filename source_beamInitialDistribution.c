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
