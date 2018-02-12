/*	
   - Works reliably only for a mixture of gasses (ionizatiojn is a small contribution).	
   - To be included after the other include files
   - Sample time loop:
   
   computeCurrentFromIonizedParticles();

   // wakefield calculation (swipeParticle);

   // output ...
   
   // calculation of E_laserField_ionization [+ other laser field components] (using an user-provided function, e.g., computeLaserField)

   solveIonizationRateEquation(); 

   for (jj=0; jj<NsubCycle_ionization; j++)
     {
       s_subCycle_ionization=s+jj*ds_subCycle_ionization; // this has to be computed explicitely 
       createIonizedElectrons();
       pushIonizedElectrons();
     }

    // laser driver push

    ** GRID STUCTURE;
    - ionization grid: Nzl_ionization x Nrl_ionization --> ionized particles dynamics is computed here (laser field is defined on this grid)
    - ions grid: Nzl_ionization x Nrl_ions_ionization --> the ionization fraction arrays (f6, f7) are defined on this grid, it describes the creation of new particles
    - wake grid
*/

// sub-cycling
int NsubCycle_ionization;  // how many sub-steps in a quasi-static step
int Nzl_ionization, Nrl_ionization, Nrl_ions_ionization;
double s_subCycle_ionization, ds_subCycle_ionization; 
double zmin_ionization, zmax_ionization, rmax_ionization, dzl_ionization, drl_ionization, drl_ions_ionization;

// ionization data
double ratio_densityIonsToBeIonizedComparedToPreionizedBackground; // e.g., for mixture H2(1-eta)+N2(eta), with N pre-ionized to +5: n5+/n0 = epsilon/(1+4*epsilon)
double kp_cgsUnits_ionization; 
double **f6_map_ionization, **f7_map_ionization, **f7_2_map_ionization, **E_laserField_ionization, **E_laserFieldLongitudinalComponent_ionization; // Nzl_ionization x Nrl_ionization
complex **a_buffer;

void (*ionizationPhysicsDetails_ionization) (double *, double *, double, double); // gives W56, W67 at a particular location

// electrons created via ionization
int Nelectron_ionization, *label_ionization;
double *x_ionization, *y_ionization, *z_ionization, *ux_ionization, *uy_ionization, *uz_ionization, *w_ionization;
float *x0_ionization, *y0_ionization, *z0_ionization, *zeta0_ionization;

void finalizeIonization()
{
  /*
    - to be placed after all the other initializations
  */
  int i;

  NsubCycle_ionization=(int)(ds/ds_subCycle_ionization+0.5);
  ds_subCycle_ionization=ds/(NsubCycle_ionization-1);

  if (zmin_ionization<zmin || zmax_ionization>zmax || rmax_ionization>rmax)
    {
      printf("ERROR: ionization grid MUST be contained in the standard grid!\n");
      exit(0);
    }
  
  Nzl_ionization=(int)((zmax_ionization-zmin_ionization)/dzl_ionization+0.5);
  dzl_ionization=(zmax_ionization-zmin_ionization)/(Nzl_ionization-1);

  Nrl_ionization=(int)(rmax_ionization/drl_ionization+0.5);
  drl_ionization=rmax_ionization/(Nrl_ionization-1);

  Nrl_ions_ionization=(int)(rmax_ionization/drl_ions_ionization+0.5);
  drl_ions_ionization=rmax_ionization/(Nrl_ions_ionization-1);
  
  f6_map_ionization=(double**)malloc(Nzl_ionization*sizeof(double*));
  f7_map_ionization=(double**)malloc(Nzl_ionization*sizeof(double*));
  f7_2_map_ionization=(double**)malloc(Nzl_ionization*sizeof(double*));

  E_laserField_ionization=(double**)malloc(Nzl_ionization*sizeof(double*));
  E_laserFieldLongitudinalComponent_ionization=(double**)malloc(Nzl_ionization*sizeof(double*));

  a_buffer=(complex**)malloc(Nzl_ionization*sizeof(complex*));
  
  for (i=0; i<Nzl_ionization; i++)
    {
      f6_map_ionization[i]=(double*)malloc(Nrl_ions_ionization*sizeof(double));
      f7_map_ionization[i]=(double*)malloc(Nrl_ions_ionization*sizeof(double));
      f7_2_map_ionization[i]=(double*)malloc(Nrl_ions_ionization*sizeof(double));

      E_laserField_ionization[i]=(double*)malloc(Nrl_ionization*sizeof(double));
      E_laserFieldLongitudinalComponent_ionization[i]=(double*)malloc(Nrl_ionization*sizeof(double));

      a_buffer[i]=(complex*)malloc(Nrl_ionization*sizeof(complex));
    }
  
  if (ionizationPhysicsDetails_ionization==NULL)
    printf("ERROR: ionizationPhysicsDetails_ionization was not secified!!\n\n");

  Nelectron_ionization=0;
  x_ionization=NULL;
  y_ionization=NULL;
  z_ionization=NULL;

  ux_ionization=NULL;
  uy_ionization=NULL;
  uz_ionization=NULL;
  
  w_ionization=NULL;

  label_ionization=NULL;

  x0_ionization=NULL;
  y0_ionization=NULL;
  z0_ionization=NULL;
  zeta0_ionization=NULL;
  
}

void interpolateOntoIonizationGrid(complex **a_laser, double zmin_laser, double zmax_laser, double rmax_laser, double dz_laser, double dr_laser, int Nz_laser, int Nr_laser)
{
  /*
    This function interpolates a_laser onto the _ionization grid.
    N.B.The _ionization grid MUST be contained inside the _laser grid!! 
  */
  
  int i, j, ii, jj;
  double zz, rr, buf_zz, buf_rr, w00, w01, w10, w11;

  for (i=0; i<Nzl_ionization; i++)
    for (j=0; j<Nrl_ionization; j++)
      {
	zz=zmin_ionization+i*dzl_ionization;
	rr=j*drl_ionization;

	buf_zz=(zz-zmin_laser)/dz_laser;
	buf_rr=rr/dr_laser;

	ii=(int)buf_zz;
	jj=(int)buf_rr;

	buf_zz-=(double)ii;
	buf_rr-=(double)jj;

	w00=(1-buf_zz)*(1-buf_rr);
	w01=(1-buf_zz)*buf_rr;
	w10=buf_zz*(1-buf_rr);
	w11=buf_zz*buf_rr;
	
	a_buffer[i][j]=w00*a_laser[ii][jj]+w01*a_laser[ii][MIN(jj+1, Nr_laser-1)]+w10*a_laser[MIN(ii+1, Nz_laser-1)][jj]+w11*a_laser[MIN(ii+1, Nz_laser-1)][MIN(jj+1, Nr_laser-1)];
      }
}

void computeLaserField(complex **a_laser, double k0_kp_laser, double zmin_laser, double zmax_laser, double rmax_laser, double dz_laser, double dr_laser, int Nz_laser, int Nr_laser, int FLAG_LONG_FIELD)
{
  /*
    This function is an example (one complex laser field only).
    
    Calculation of the laser electric field from envelope, a, to be used in the formula for Wi->j and for the equation of motion.

    This function takes the a_laser field defined on the "laser grid", interpolates it onto the "ionization grid"
    
    The longitudinal fields can be computed as follows:

    Ez=E_laserFieldLongitudinalComponent_ionization*cos(theta)
    Bz=E_laserFieldLongitudinalComponent_ionization*sin(theta)
    
    where theta=atan2(y, x)
  */
  int i, j;
  double zeta;
  complex buffer;

  // intepolate the field "a_laser" onto the _ionization grid
  interpolateOntoIonizationGrid(a_laser, zmin_laser, zmax_laser, rmax_laser, dz_laser, dr_laser, Nz_laser, Nr_laser);
  
  for (i=1; i<Nzl_ionization-1; i++)
    for (j=0; j<Nrl_ionization-1; j++)
      {
	zeta=zmin_ionization+i*dzl_ionization;
	buffer=0.5*((a_buffer[i+1][j]-a_buffer[i-1][j])/(2*dzl_ionization)+I*k0_kp_laser*a_buffer[i][j])*cexp(I*k0_kp_laser*zeta);
	E_laserField_ionization[i][j]=creal(buffer+conj(buffer));

	if (FLAG_LONG_FIELD==YES)
	  {
	    if (j==0)
	      E_laserFieldLongitudinalComponent_ionization[i][j]=0;
	    else
	      {
		buffer=-0.5*((a_buffer[i][j+1]-a_buffer[i][j-1])/(2*drl_ionization))*cexp(I*k0_kp_laser*zeta);
		E_laserFieldLongitudinalComponent_ionization[i][j]=creal(buffer+conj(buffer));
	      }
	  }
	else
	  E_laserFieldLongitudinalComponent_ionization[i][j]=0;
      }

  // setting the values in j=Nr_ionization-1
  for (i=1; i<Nzl_ionization-1; i++)
    {
      E_laserField_ionization[i][Nrl_ionization-1]=2*E_laserField_ionization[i][Nrl_ionization-2]-E_laserField_ionization[i][Nrl_ionization-3];
      E_laserFieldLongitudinalComponent_ionization[i][Nrl_ionization-1]=2*E_laserFieldLongitudinalComponent_ionization[i][Nrl_ionization-2]-E_laserFieldLongitudinalComponent_ionization[i][Nrl_ionization-3];
    }

  // setting the values in i=0, =Nzl_ionization-1
  for (j=0; j<Nrl_ionization; j++)
    {
      E_laserField_ionization[0][j]=2*E_laserField_ionization[1][j]-E_laserField_ionization[2][j];
      E_laserFieldLongitudinalComponent_ionization[0][j]=2*E_laserFieldLongitudinalComponent_ionization[1][j]-E_laserFieldLongitudinalComponent_ionization[2][j];
      E_laserField_ionization[Nzl_ionization-1][j]=2*E_laserField_ionization[Nzl_ionization-2][j]-E_laserField_ionization[Nzl_ionization-3][j];
      E_laserFieldLongitudinalComponent_ionization[Nzl_ionization-1][j]=2*E_laserFieldLongitudinalComponent_ionization[Nzl_ionization-2][j]-E_laserFieldLongitudinalComponent_ionization[Nzl_ionization-3][j];
    }
}

void compute_tmp_ionization(double *f6_tmp, double *f7_tmp, double *d_f6, double *d_f7, double h)
{
  int j;
  
  for (j=0; j<Nrl_ions_ionization; j++)
    {
      f6_tmp[j]=h*d_f6[j];
      f7_tmp[j]=h*d_f7[j];
    }
}

void compute_df_ionization(double *d_f6, double *d_f7, double *f6_tmp, double *f7_tmp, double zeta)
{
  int j;
  double W56, W67;
  
  for (j=0; j<Nrl_ions_ionization; j++)
    {
      ionizationPhysicsDetails_ionization(&W56, &W67, zeta, j*drl_ions_ionization);

      d_f6[j]=-W56*(1-f6_tmp[j]-f7_tmp[j])+W67*f6_tmp[j];
      d_f7[j]=-W67*f6_tmp[j];
    }
}

void solveIonizationRateEquation()
{
  /*
    N.B.1. The probability is not accumulated form zeta=+inf to zeta=-inf, the probability of the transition is computed for only one step:
    
    1. f5(zeta)=1, f6(zeta)=0, f7(zeta)=0
    2. perform one RK4 step from zeta to zeta-dzeta
    3. generate particle at zeta-dzeta according to f5(zeta-dzeta), f6(zeta-dzeta), f7(zeta-dzeta)
 
    N.B.2. Before calling this function the user MUST compute E_laserField_ionization (this depends on the details of the physics)

    calcoli con RK4
    
  */
  
  int i, j, jj;
  double zeta, r, ww, W56, W67;
  static int FIRST_TIME=1;
  static double *d_f6[4], *d_f7[4], *f6_tmp, *f7_tmp;

  if (FIRST_TIME==1)
    {
      for (i=0; i<4; i++)
	{
	  d_f6[i]=malloc(Nrl_ions_ionization*sizeof(double));
	  d_f7[i]=malloc(Nrl_ions_ionization*sizeof(double));
	}

      f6_tmp=malloc(Nrl_ions_ionization*sizeof(double));
      f7_tmp=malloc(Nrl_ions_ionization*sizeof(double));
      
      FIRST_TIME=0;
    }
  
  memset((void*)f6_map_ionization[Nzl_ionization-1], 0, Nrl_ions_ionization*sizeof(double));  
  memset((void*)f7_map_ionization[Nzl_ionization-1], 0, Nrl_ions_ionization*sizeof(double));  

  for (i=Nzl_ionization-1; i>0; i--)
    {
      zeta=zmin_ionization+i*dzl_ionization;
      
      compute_df_ionization(d_f6[0], d_f7[0], f6_map_ionization[Nzl_ionization-1], f7_map_ionization[Nzl_ionization-1], zeta); // values equal to zero
      compute_tmp_ionization(f6_tmp, f7_tmp, d_f6[0], d_f7[0], -0.5*dzl_ionization);

      compute_df_ionization(d_f6[1], d_f7[1], f6_tmp, f7_tmp, zeta-0.5*dzl_ionization);
      compute_tmp_ionization(f6_tmp, f7_tmp, d_f6[1], d_f7[1], -0.5*dzl_ionization);
      
      compute_df_ionization(d_f6[2], d_f7[2], f6_tmp, f7_tmp, zeta-0.5*dzl_ionization);
      compute_tmp_ionization(f6_tmp, f7_tmp, d_f6[2], d_f7[2], -dzl_ionization);

      compute_df_ionization(d_f6[3], d_f7[3], f6_tmp, f7_tmp, zeta-dzl_ionization);
      
      for (j=0; j<Nrl_ions_ionization; j++)
	{
	  f6_map_ionization[i-1][j]=-(d_f6[0][j]+2*d_f6[1][j]+2*d_f6[2][j]+d_f6[3][j])*dzl_ionization/6.;
	  f7_map_ionization[i-1][j]=-(d_f7[0][j]+2*d_f7[1][j]+2*d_f7[2][j]+d_f7[3][j])*dzl_ionization/6.;
	}
    }

  for (i=0; i<Nzl_ionization; i++)
    for (j=0; j<Nrl_ions_ionization; j++)
      {
	ionizationPhysicsDetails_ionization(&W56, &W67, zmin_ionization+i*dzl_ionization, j*drl_ions_ionization);
	f7_2_map_ionization[i][j]=1-exp(-W67*dzl_ionization);
      }
}

void solveIonizationRateEquation_2()
{
  /*
    N.B.1. The probability is not accumulated form zeta=+inf to zeta=-inf, the probability of the transition is computed for only one step:
    
    1. f5(zeta)=1, f6(zeta)=0, f7(zeta)=0
    2. perform one RK2 step from zeta to zeta-dzeta
    3. generate particle at zeta-dzeta according to f5(zeta-dzeta), f6(zeta-dzeta), f7(zeta-dzeta)
 
    N.B.2. Before calling this function the user MUST compute E_laserField_ionization (this depends on the details of the physics)
    
    calcoli con RK2
    
  */
  
  int i, j, jj;
  double zeta, r, ww, W56, W67;
  static int FIRST_TIME=1;
  static double *d_f6, *d_f7, *f6_tmp, *f7_tmp;

  if (FIRST_TIME==1)
    {
      d_f6=malloc(Nrl_ions_ionization*sizeof(double));
      d_f7=malloc(Nrl_ions_ionization*sizeof(double));

      f6_tmp=malloc(Nrl_ions_ionization*sizeof(double));
      f7_tmp=malloc(Nrl_ions_ionization*sizeof(double));
      
      FIRST_TIME=0;
    }
  
  memset((void*)f6_map_ionization[Nzl_ionization-1], 0, Nrl_ions_ionization*sizeof(double));  
  memset((void*)f7_map_ionization[Nzl_ionization-1], 0, Nrl_ions_ionization*sizeof(double));  

  for (i=Nzl_ionization-1; i>0; i--)
    {
      zeta=zmin_ionization+i*dzl_ionization;
      
      compute_df_ionization(d_f6, d_f7, f6_map_ionization[Nzl_ionization-1], f7_map_ionization[Nzl_ionization-1], zeta); // values equal to zero
      compute_tmp_ionization(f6_tmp, f7_tmp, d_f6, d_f7, -0.5*dzl_ionization);

      compute_df_ionization(d_f6, d_f7, f6_tmp, f7_tmp, zeta-0.5*dzl_ionization);
      
      for (j=0; j<Nrl_ions_ionization; j++)
	{
	  f6_map_ionization[i-1][j]=-d_f6[j]*dzl_ionization;
	  f7_map_ionization[i-1][j]=-d_f7[j]*dzl_ionization;
	}
    }

  for (i=0; i<Nzl_ionization; i++)
    for (j=0; j<Nrl_ions_ionization; j++)
      {
	ionizationPhysicsDetails_ionization(&W56, &W67, zmin_ionization+i*dzl_ionization, j*drl_ions_ionization);
	f7_2_map_ionization[i][j]=1-exp(-W67*dzl_ionization);
      }
}

void solveIonizationRateEquation_3()
{
  /*
    N.B.1. The probability is not accumulated form zeta=+inf to zeta=-inf, the probability of the transition is computed for only one step:
    
    1. f5(zeta)=1, f6(zeta)=0, f7(zeta)=0
    2. perform one RK2 step from zeta to zeta-dzeta
    3. generate particle at zeta-dzeta according to f5(zeta-dzeta), f6(zeta-dzeta), f7(zeta-dzeta)
 
    N.B.2. Before calling this function the user MUST compute E_laserField_ionization (this depends on the details of the physics)

    calcoli smart che risolvono il problema del rate negativo
    
  */
  
  int i, j, jj;
  double zeta, r, ww, W56, W67;
  double u, d_f7, f7_tmp;
    
  for (i=0; i<Nzl_ionization; i++)
    for (j=0; j<Nrl_ions_ionization; j++)
      {
	ionizationPhysicsDetails_ionization(&W56, &W67, zmin_ionization+i*dzl_ionization, j*drl_ions_ionization);
	
	f6_map_ionization[i][j]  =1-exp(-W56*dzl_ionization); // --> u=f6+f7
	f7_2_map_ionization[i][j]=1-exp(-W67*dzl_ionization);
      }
  
  memset((void*)f7_map_ionization[Nzl_ionization-1], 0, Nrl_ions_ionization*sizeof(double));  

  for (i=Nzl_ionization-1; i>0; i--)
    for (j=0; j<Nrl_ions_ionization; j++)
      {
	ionizationPhysicsDetails_ionization(&W56, &W67, zmin_ionization+i*dzl_ionization, j*drl_ions_ionization);
	u=f6_map_ionization[i][j];
	d_f7=-W67*u;
	f7_tmp=-0.5*dzl_ionization*d_f7;

	u=0.5*(f6_map_ionization[i-1][j]+f6_map_ionization[i][j]);
	if (u<=f7_tmp)
	  d_f7=0;
	else
	  {
	    ionizationPhysicsDetails_ionization(&W56, &W67, zmin_ionization+i*dzl_ionization-0.5*dzl_ionization, j*drl_ions_ionization);
	    d_f7=-W67*(u-f7_tmp);
	  }
	
	f7_map_ionization[i-1][j]=-d_f7*dzl_ionization;   
      }
  
  for (i=0; i<Nzl_ionization; i++)
    for (j=0; j<Nrl_ions_ionization; j++)
      f6_map_ionization[i][j]-=f7_map_ionization[i][j];
  
  for (i=0; i<Nzl_ionization; i++)
    for (j=0; j<Nrl_ions_ionization; j++)
      f6_map_ionization[i][j]=MAX(f6_map_ionization[i][j], 0);
}

void addIonizedElectron(double zeta_0, double r_0, double charge_0)
{
  int i;
  double phi, ww;

  ww=charge_0*(2*M_PI*fabs(r_0)*drl_ions_ionization*ds_subCycle_ionization)*ratio_densityIonsToBeIonizedComparedToPreionizedBackground*backgroundDensityProfile(s, fabs(r_0));
  if (ww==0.) return;
  
  Nelectron_ionization++;

  x_ionization=(double*)realloc(x_ionization, Nelectron_ionization*sizeof(double));
  y_ionization=(double*)realloc(y_ionization, Nelectron_ionization*sizeof(double));
  z_ionization=(double*)realloc(z_ionization, Nelectron_ionization*sizeof(double));

  ux_ionization=(double*)realloc(ux_ionization, Nelectron_ionization*sizeof(double));
  uy_ionization=(double*)realloc(uy_ionization, Nelectron_ionization*sizeof(double));
  uz_ionization=(double*)realloc(uz_ionization, Nelectron_ionization*sizeof(double));

  w_ionization=(double*)realloc(w_ionization, Nelectron_ionization*sizeof(double));

  label_ionization=(int*)realloc(label_ionization, Nelectron_ionization*sizeof(int));

  x0_ionization=(float*)realloc(x0_ionization, Nelectron_ionization*sizeof(float));
  y0_ionization=(float*)realloc(y0_ionization, Nelectron_ionization*sizeof(float));
  z0_ionization=(float*)realloc(z0_ionization, Nelectron_ionization*sizeof(float));
  zeta0_ionization=(float*)realloc(zeta0_ionization, Nelectron_ionization*sizeof(float));
  
  i=Nelectron_ionization-1;

  phi=2*M_PI*random()/(double)RAND_MAX;
  x_ionization[i]=r_0*cos(phi);
  y_ionization[i]=r_0*sin(phi);
  z_ionization[i]=zeta_0;
  ux_ionization[i]=uy_ionization[i]=uz_ionization[i]=0;

  w_ionization[i]=ww;

  label_ionization[i]=i;

  x0_ionization[i]=x_ionization[i];
  y0_ionization[i]=y_ionization[i];
  z0_ionization[i]=z_ionization[i]+s_subCycle_ionization;
  zeta0_ionization[i]=z_ionization[i];
}
  
void createIonizedElectrons()
{
  int k, i, ii, n, flag_56;
  double rnd, x_01, x_12, rnd_dspl1, rnd_dspl2;
  
  for (k=0; k<Nrl_ions_ionization; k++)
    {
      flag_56=0;
      for (i=Nzl_ionization-1; i>=0; i--)
	{
	  x_01=1-f6_map_ionization[i][k]-f7_map_ionization[i][k];
	  x_12=1-f7_map_ionization[i][k];

	  rnd=random()/(double)RAND_MAX;
	  
	  if (rnd<=x_01)
	    continue;
	  else if (rnd<=x_12)
	    {
	      rnd_dspl1=(double)random()/(double)RAND_MAX-0.5;
	      rnd_dspl2=(double)random()/(double)RAND_MAX-0.5;
	      
	      addIonizedElectron(zmin_ionization+dzl_ionization*(i+rnd_dspl1), drl_ions_ionization*(k+rnd_dspl2), 1.); // transition 5->6, 1 electron
	      flag_56=1;    
	      break;
	    }
	  else if (rnd>x_12)
	    {
	      rnd_dspl1=(double)random()/(double)RAND_MAX-0.5;
	      rnd_dspl2=(double)random()/(double)RAND_MAX-0.5;
	      addIonizedElectron(zmin_ionization+dzl_ionization*(i+rnd_dspl1), drl_ions_ionization*(k+rnd_dspl2), 1.); // transition 5->7, 1st electron 

	      rnd_dspl1=(double)random()/(double)RAND_MAX-0.5;
	      rnd_dspl2=(double)random()/(double)RAND_MAX-0.5;
	      addIonizedElectron(zmin_ionization+dzl_ionization*(i+rnd_dspl1), drl_ions_ionization*(k+rnd_dspl2), 1.); // transition 5->7, 2nd electron
	      break;
	    }
	}

      // transition 6->7 (after 5->6) and generation of 1 electorn
      if (flag_56)
	for (ii=i-1; ii>=0; ii--)
	  {
	    rnd=random()/(double)RAND_MAX;
	    
	    if (rnd<f7_2_map_ionization[ii][k])
	      {
		rnd_dspl1=(double)random()/(double)RAND_MAX-0.5;
		rnd_dspl2=(double)random()/(double)RAND_MAX-0.5;
		addIonizedElectron(zmin_ionization+dzl_ionization*(ii+rnd_dspl1), drl_ions_ionization*(k+rnd_dspl2), 1.); // transition 5->6, 1 electron
		break;
	      }
	  }
    }
}

int interpolateIonizedParticle(double *int_Ex, double *int_Ey, double *int_Ez, double *int_Bx, double *int_By, double *int_Bz, double tmp_x, double tmp_y, double tmp_z)
{
  int i, j;
  double buf_zz, buf_rr, w00, w01, w10, w11;
  double r, theta, Ez_w, Er_w, Bphi_w, Ftrans_l, Flong_l;

  r=hypot(tmp_x, tmp_y);
  theta=atan2(tmp_y, tmp_x);
  
  *int_Ex=*int_Ey=*int_Ez=*int_Bx=*int_By=*int_Bz=0;

  // wake component
  buf_zz=(tmp_z-zmin)/dz;
  i=(int)buf_zz;
  buf_zz-=(double)i;

  buf_rr=r/dr;
  j=(int)buf_rr;
  buf_rr-=(double)j;

  if (i<0 || j>=Nr-1) // particle outside simulation box --> error message
    return 1;
  
  w00=(1-buf_zz)*(1-buf_rr);
  w01=(1-buf_zz)*buf_rr;
  w10=buf_zz*(1-buf_rr);
  w11=buf_zz*buf_rr;

  Bphi_w=w00*Bphi[i][j]+w01*Bphi[i][j+1]+w10*Bphi[i+1][j]+w11*Bphi[i+1][j+1];
  Ez_w  =w00*Ez[i][j]  +w01*Ez[i][j+1]  +w10*Ez[i+1][j]  +w11*Ez[i+1][j+1];
  Er_w  =(w00*Er_Bphi[i][j]  +w01*Er_Bphi[i][j+1]  +w10*Er_Bphi[i+1][j]  +w11*Er_Bphi[i+1][j+1])+Bphi_w;
  
  /* 
     r^=[ cos(theta), sin(theta), 0]
     n^=[-sin(theta), cos(theta), 0]
     z^=[           0,         0, 1]
  */

  *int_Ex=cos(theta)*Er_w;
  *int_Ey=sin(theta)*Er_w;
  *int_Ez=Ez_w;
  *int_Bx=-sin(theta)*Bphi_w;
  *int_By= cos(theta)*Bphi_w;
  *int_Bz=0;
  
  // laser component
  buf_zz=(tmp_z-zmin_ionization)/dzl_ionization;
  i=(int)buf_zz;
  buf_zz-=(double)i;

  buf_rr=r/drl_ionization;
  j=(int)buf_rr;
  buf_rr-=(double)j;

  if (i<0 || i>=Nzl_ionization-1 || j>=Nrl_ionization-1) // particle outside of the domain where the laser field is defined
    return 0;

  w00=(1-buf_zz)*(1-buf_rr);
  w01=(1-buf_zz)*buf_rr;
  w10=buf_zz*(1-buf_rr);
  w11=buf_zz*buf_rr;
  
  Ftrans_l=w00*E_laserField_ionization[i][j]+w01*E_laserField_ionization[i][j+1]+w10*E_laserField_ionization[i+1][j]+w11*E_laserField_ionization[i+1][j+1];
  Flong_l=w00*E_laserFieldLongitudinalComponent_ionization[i][j]+w01*E_laserFieldLongitudinalComponent_ionization[i][j+1]+w10*E_laserFieldLongitudinalComponent_ionization[i+1][j]+w11*E_laserFieldLongitudinalComponent_ionization[i+1][j+1];

  *int_Ex+=Ftrans_l;
  *int_By+=Ftrans_l;
  *int_Ez+=Flong_l*cos(theta);
  *int_Bz+=Flong_l*sin(theta);
  
  return 0;
}

void pushIonizedElectrons()
{
  int n, lost, flag;
  double gamma, int_Ex, int_Ey, int_Ez, int_Bx, int_By, int_Bz;
  double dx_ionization[4], dy_ionization[4], dz_ionization[4], dux_ionization[4], duy_ionization[4], duz_ionization[4],
    beta_x, beta_y, beta_z, tmp_x, tmp_y, tmp_z, tmp_ux, tmp_uy, tmp_uz;
  
  lost=0;
  for (n=0; n<Nelectron_ionization; n++)
    {
      // step 1
      flag=interpolateIonizedParticle(&int_Ex, &int_Ey, &int_Ez, &int_Bx, &int_By, &int_Bz, x_ionization[n], y_ionization[n], z_ionization[n]);
      gamma=sqrt(1+ux_ionization[n]*ux_ionization[n]+uy_ionization[n]*uy_ionization[n]+uz_ionization[n]*uz_ionization[n]);
      dx_ionization[0]=beta_x=ux_ionization[n]/gamma;
      dy_ionization[0]=beta_y=uy_ionization[n]/gamma;
      beta_z=uz_ionization[n]/gamma;
      dz_ionization[0]=beta_z-1;
      dux_ionization[0]=-(int_Ex+beta_y*int_Bz-beta_z*int_By);
      duy_ionization[0]=-(int_Ey+beta_z*int_Bx-beta_x*int_Bz);
      duz_ionization[0]=-(int_Ez+beta_x*int_By-beta_y*int_Bx);
      
      // step 2
      tmp_x=x_ionization[n]+dx_ionization[0]*ds_subCycle_ionization/2.;
      tmp_y=y_ionization[n]+dy_ionization[0]*ds_subCycle_ionization/2.;
      tmp_z=z_ionization[n]+dz_ionization[0]*ds_subCycle_ionization/2.;
      tmp_ux=ux_ionization[n]+dux_ionization[0]*ds_subCycle_ionization/2.;
      tmp_uy=uy_ionization[n]+duy_ionization[0]*ds_subCycle_ionization/2.;
      tmp_uz=uz_ionization[n]+duz_ionization[0]*ds_subCycle_ionization/2.;
      flag=interpolateIonizedParticle(&int_Ex, &int_Ey, &int_Ez, &int_Bx, &int_By, &int_Bz, tmp_x, tmp_y, tmp_z);
      if (flag)
	{
	  lost++;
	  continue;
	}
      gamma=sqrt(1+tmp_ux*tmp_ux+tmp_uy*tmp_uy+tmp_uz*tmp_uz);
      dx_ionization[1]=beta_x=tmp_ux/gamma;
      dy_ionization[1]=beta_y=tmp_uy/gamma;
      beta_z=tmp_uz/gamma;
      dz_ionization[1]=beta_z-1;
      dux_ionization[1]=-(int_Ex+beta_y*int_Bz-beta_z*int_By);
      duy_ionization[1]=-(int_Ey+beta_z*int_Bx-beta_x*int_Bz);
      duz_ionization[1]=-(int_Ez+beta_x*int_By-beta_y*int_Bx);
       
      // step 3
      tmp_x=x_ionization[n]+dx_ionization[1]*ds_subCycle_ionization/2.;
      tmp_y=y_ionization[n]+dy_ionization[1]*ds_subCycle_ionization/2.;
      tmp_z=z_ionization[n]+dz_ionization[1]*ds_subCycle_ionization/2.;
      tmp_ux=ux_ionization[n]+dux_ionization[1]*ds_subCycle_ionization/2.;
      tmp_uy=uy_ionization[n]+duy_ionization[1]*ds_subCycle_ionization/2.;
      tmp_uz=uz_ionization[n]+duz_ionization[1]*ds_subCycle_ionization/2.;
      flag=interpolateIonizedParticle(&int_Ex, &int_Ey, &int_Ez, &int_Bx, &int_By, &int_Bz, tmp_x, tmp_y, tmp_z);
      if (flag)
	{
	  lost++;
	  continue;
	}
      gamma=sqrt(1+tmp_ux*tmp_ux+tmp_uy*tmp_uy+tmp_uz*tmp_uz);
      dx_ionization[2]=beta_x=tmp_ux/gamma;
      dy_ionization[2]=beta_y=tmp_uy/gamma;
      beta_z=tmp_uz/gamma;
      dz_ionization[2]=beta_z-1;
      dux_ionization[2]=-(int_Ex+beta_y*int_Bz-beta_z*int_By);
      duy_ionization[2]=-(int_Ey+beta_z*int_Bx-beta_x*int_Bz);
      duz_ionization[2]=-(int_Ez+beta_x*int_By-beta_y*int_Bx);

      // step 4
      tmp_x=x_ionization[n]+dx_ionization[2]*ds_subCycle_ionization;
      tmp_y=y_ionization[n]+dy_ionization[2]*ds_subCycle_ionization;
      tmp_z=z_ionization[n]+dz_ionization[2]*ds_subCycle_ionization;
      tmp_ux=ux_ionization[n]+dux_ionization[2]*ds_subCycle_ionization;
      tmp_uy=uy_ionization[n]+duy_ionization[2]*ds_subCycle_ionization;
      tmp_uz=uz_ionization[n]+duz_ionization[2]*ds_subCycle_ionization;
      flag=interpolateIonizedParticle(&int_Ex, &int_Ey, &int_Ez, &int_Bx, &int_By, &int_Bz, tmp_x, tmp_y, tmp_z);
      if (flag)
	{
	  lost++;
	  continue;
	}
      gamma=sqrt(1+tmp_ux*tmp_ux+tmp_uy*tmp_uy+tmp_uz*tmp_uz);
      dx_ionization[3]=beta_x=tmp_ux/gamma;
      dy_ionization[3]=beta_y=tmp_uy/gamma;
      beta_z=tmp_uz/gamma;
      dz_ionization[3]=beta_z-1;
      dux_ionization[3]=-(int_Ex+beta_y*int_Bz-beta_z*int_By);
      duy_ionization[3]=-(int_Ey+beta_z*int_Bx-beta_x*int_Bz);
      duz_ionization[3]=-(int_Ez+beta_x*int_By-beta_y*int_Bx);

      x_ionization[n]+=(dx_ionization[0]+2*dx_ionization[1]+2*dx_ionization[2]+dx_ionization[3])*ds_subCycle_ionization/6.;
      y_ionization[n]+=(dy_ionization[0]+2*dy_ionization[1]+2*dy_ionization[2]+dy_ionization[3])*ds_subCycle_ionization/6.;
      z_ionization[n]+=(dz_ionization[0]+2*dz_ionization[1]+2*dz_ionization[2]+dz_ionization[3])*ds_subCycle_ionization/6.;
      ux_ionization[n]+=(dux_ionization[0]+2*dux_ionization[1]+2*dux_ionization[2]+dux_ionization[3])*ds_subCycle_ionization/6.;
      uy_ionization[n]+=(duy_ionization[0]+2*duy_ionization[1]+2*duy_ionization[2]+duy_ionization[3])*ds_subCycle_ionization/6.;
      uz_ionization[n]+=(duz_ionization[0]+2*duz_ionization[1]+2*duz_ionization[2]+duz_ionization[3])*ds_subCycle_ionization/6.;

      if (z_ionization[n]<zmin || hypot(x_ionization[n], y_ionization[n])>rmax-dr) // XXX occhio alla uscita
	{
	  lost++;
	  continue;
	}

      x_ionization[n-lost]=x_ionization[n];
      y_ionization[n-lost]=y_ionization[n];
      z_ionization[n-lost]=z_ionization[n];

      ux_ionization[n-lost]=ux_ionization[n];
      uy_ionization[n-lost]=uy_ionization[n];
      uz_ionization[n-lost]=uz_ionization[n];

      w_ionization[n-lost]=w_ionization[n];
      
      label_ionization[n-lost]=label_ionization[n];
      
      x0_ionization[n-lost]=x0_ionization[n];
      y0_ionization[n-lost]=y0_ionization[n];
      z0_ionization[n-lost]=z0_ionization[n];
      zeta0_ionization[n-lost]=zeta0_ionization[n];

    }
  
  Nelectron_ionization-=lost;
}

void pushIonizedElectrons_nomotion()
{
  int n, lost;
  
  for (n=0; n<Nelectron_ionization; n++)
    z_ionization[n]-=ds_subCycle_ionization;
  
  lost=0;
  for (n=0; n<Nelectron_ionization; n++)
    {
      if (z_ionization[n]<zmin || hypot(x_ionization[n], y_ionization[n])>rmax-dr)
	{
	  lost++;
	  continue;
	}

      x_ionization[n-lost]=x_ionization[n];
      y_ionization[n-lost]=y_ionization[n];
      z_ionization[n-lost]=z_ionization[n];

      ux_ionization[n-lost]=ux_ionization[n];
      uy_ionization[n-lost]=uy_ionization[n];
      uz_ionization[n-lost]=uz_ionization[n];

      w_ionization[n-lost]=w_ionization[n];

      label_ionization[n-lost]=label_ionization[n];

      x0_ionization[n-lost]=x0_ionization[n];
      y0_ionization[n-lost]=y0_ionization[n];
      z0_ionization[n-lost]=z0_ionization[n];
      zeta0_ionization[n-lost]=zeta0_ionization[n];
    }
  
  Nelectron_ionization-=lost;
}

double **densityIonized=NULL;
void densityIonizedElectrons()
{
  /*
    Define the density of ionized particles in the same grid as the wakefield
  */
  int i, j, n;
  double w_z, w_r, w00, w01, w10, w11, buffer1, buffer2;
  
  if (densityIonized==NULL)
    densityIonized=allocatePointer(NZ_GHOST);

  for(i=0; i<Nz; i++)
    memset(densityIonized[i], 0, Nr*sizeof(double));

   for (n=0; n<Nelectron_ionization; n++)
     {
       w_z=(z_ionization[n]-zmin)/dz;
       i=(int)w_z;
       w_z-=(double)i;

       w_r=hypot(x_ionization[n], y_ionization[n])/dr;
       j=(int)w_r;
       w_r-=(double)j;

       w00=(1-w_z)*(1-w_r);
       w01=(1-w_z)*w_r;
       w10=w_z*(1-w_r);
       w11=w_z*w_r;

       buffer1=w_ionization[n]/(j==0?(1./6.):j); 
       buffer2=w_ionization[n]/(j+1.);	  

       densityIonized[i  ][j  ]+=w00*buffer1;
       densityIonized[i  ][j+1]+=w01*buffer2;
       densityIonized[i+1][j  ]+=w10*buffer1;
       densityIonized[i+1][j+1]+=w11*buffer2;
     }
   
   for (i=0; i<Nz; i++)
     for (j=0; j<Nr; j++)
      densityIonized[i][j]/=(2*M_PI*dz*dr*dr);
}

void computeCurrentFromIonizedParticles()
{
  /*
    - linear shape function
    - no check of the boundaries!
    - Jz_ionization grid is the wakefield grid
  */
  int i, j, n;
  double w_z, w_r, w00, w01, w10, w11, buffer1, buffer2;
  double beta_z, gamma;
  
  for(i=0; i<Nz; i++)
    memset(Jz_ionization[i], 0, Nr*sizeof(double));
  
  for (n=0; n<Nelectron_ionization; n++)
    {
      w_z=(z_ionization[n]-zmin)/dz;
      i=(int)w_z;
      w_z-=(double)i;
      
      w_r=hypot(x_ionization[n], y_ionization[n])/dr;
      j=(int)w_r;
      w_r-=(double)j;
      
      w00=(1-w_z)*(1-w_r);
      w01=(1-w_z)*w_r;
      w10=w_z*(1-w_r);
      w11=w_z*w_r;
      
      gamma=sqrt(1+ux_ionization[n]*ux_ionization[n]+uy_ionization[n]*uy_ionization[n]+uz_ionization[n]*uz_ionization[n]);       
      beta_z=uz_ionization[n]/gamma;

      buffer1=w_ionization[n]*beta_z/(j==0?(1./6.):j); 
      buffer2=w_ionization[n]*beta_z/(j+1.);	  
      
      Jz_ionization[i  ][j  ]+=w00*buffer1;
      Jz_ionization[i  ][j+1]+=w01*buffer2;
      Jz_ionization[i+1][j  ]+=w10*buffer1;
      Jz_ionization[i+1][j+1]+=w11*buffer2;
    }
  
  for (i=0; i<Nz; i++)
    for (j=0; j<Nr; j++)
      Jz_ionization[i][j]*=(-1./(2.*M_PI*dz*dr*dr)); // the factor "-1" takes into account the electron charge 
}

// ==================================================================================================================================================

#define ALPHA_FINE_STRUCTURE_CONSTANT (1./137.)
#define U_H 13.598 // eV

// Nitrogen 
#define U_N0to1_ionization 14.534 // eV
#define U_N5to6_ionization 552.057 // eV
#define U_N6to7_ionization 667.029 // eV

void ionizationPhysicsDetails_Nitrogen(double *W56, double *W67, double zz, double rr)
{
  /* 
     This function sets the values of W56 and/or W67 starting from the amplitude of the ionizing electric field specified in E_laserField_ionization.
     N.B.1. E_laserField_ionization is normalized with E0, so it might be required to denormalize it using E0_amplitudeNormalizationField_ionization (if needed)
  */
  int ii, jj;
  double buf_zz, buf_rr, w00, w01, w10, w11;

  double omega_a_norm, E_atomic_norm, C;
  double n_star, n0_star, l, l_star, Eionization;
  double buf1, buf2;

  if (zz<zmin_ionization || zz>zmax_ionization-dzl_ionization || rr>rmax_ionization-drl)
    {
      *W56=0;
      *W67=0;
      return;
    }

  buf_zz=(zz-zmin_ionization)/dzl_ionization;
  buf_rr=rr/drl_ionization;
  
  ii=(int)buf_zz;
  jj=(int)buf_rr;
  
  buf_zz-=(double)ii;
  buf_rr-=(double)jj;
  
  w00=(1-buf_zz)*(1-buf_rr);
  w01=(1-buf_zz)*buf_rr;
  w10=buf_zz*(1-buf_rr);
  w11=buf_zz*buf_rr;

  Eionization=w00*E_laserField_ionization[ii][jj]+w01*E_laserField_ionization[ii][jj+1]+w10*E_laserField_ionization[ii+1][jj]+w11*E_laserField_ionization[ii+1][jj+1];
  
  n0_star=sqrt(U_H/U_N0to1_ionization);
  l_star=n0_star-1;
  l=0; // electrons #5 and #6 are in the s-orbitals (are the electrons in the inner, n=1, shell)

  omega_a_norm=pow(ALPHA_FINE_STRUCTURE_CONSTANT,3)/(kp_cgsUnits_ionization*R0_cgs);
  E_atomic_norm=pow(ALPHA_FINE_STRUCTURE_CONSTANT,4)/(kp_cgsUnits_ionization*R0_cgs);
    
  // 5->6
  n_star=6*sqrt(U_H/U_N5to6_ionization);

  C=(2*l+1)*pow(2, 4.*n_star-4.)/(n_star*tgamma(n_star+l_star+1)*tgamma(n_star-l_star));
  buf1=pow(U_N5to6_ionization/U_H, (6.*n_star-1.)/2.)/pow(fabs(Eionization/E_atomic_norm), 2*n_star-1);
  buf2=2*pow(U_N5to6_ionization/U_H, 1.5)/(3*fabs(Eionization/E_atomic_norm));
  *W56=4*omega_a_norm*C*buf1*exp(-buf2);

  // 6->7
  n_star=7*sqrt(U_H/U_N6to7_ionization);

  C=(2*l+1)*pow(2, 4.*n_star-4.)/(n_star*tgamma(n_star+l_star+1)*tgamma(n_star-l_star));
  buf1=pow(U_N6to7_ionization/U_H, (6.*n_star-1.)/2.)/pow(fabs(Eionization/E_atomic_norm), 2*n_star-1);
  buf2=2*pow(U_N6to7_ionization/U_H, 1.5)/(3*fabs(Eionization/E_atomic_norm));
  *W67=4*omega_a_norm*C*buf1*exp(-buf2);
}

void ionizationPhysicsDetails_NitrogenSimplified(double *W56, double *W67, double zz, double rr)
{
  /* 
     This function sets the values of W56 and/or W67 starting from the amplitude of the ionizing electric field specified in E_laserField_ionization.
     N.B.1. E_laserField_ionization is normalized with E0, so it might be required to denormalize it using E0_amplitudeNormalizationField_ionization (if needed)
  */
  int ii, jj;
  double buf_zz, buf_rr, w00, w01, w10, w11;

  double omega_a_norm, E_atomic_norm, C;
  double n_star, n0_star, l, l_star, Eionization;
  double buf1, buf2;

  if (zz<zmin_ionization || zz>zmax_ionization-dzl_ionization || rr>rmax_ionization-drl)
    {
      *W56=0;
      *W67=0;
      return;
    }

  buf_zz=(zz-zmin_ionization)/dzl_ionization;
  buf_rr=rr/drl_ionization;
  
  ii=(int)buf_zz;
  jj=(int)buf_rr;
  
  buf_zz-=(double)ii;
  buf_rr-=(double)jj;
  
  w00=(1-buf_zz)*(1-buf_rr);
  w01=(1-buf_zz)*buf_rr;
  w10=buf_zz*(1-buf_rr);
  w11=buf_zz*buf_rr;

  Eionization=w00*E_laserField_ionization[ii][jj]+w01*E_laserField_ionization[ii][jj+1]+w10*E_laserField_ionization[ii+1][jj]+w11*E_laserField_ionization[ii+1][jj+1];
  
  n0_star=sqrt(U_H/U_N0to1_ionization);
  l_star=n0_star-1;
  l=0; // electrons #5 and #6 are in the s-orbitals (are the electrons in the inner, n=1, shell)

  omega_a_norm=pow(ALPHA_FINE_STRUCTURE_CONSTANT,3)/(kp_cgsUnits_ionization*R0_cgs);
  E_atomic_norm=pow(ALPHA_FINE_STRUCTURE_CONSTANT,4)/(kp_cgsUnits_ionization*R0_cgs);
    
  // 5->6
  n_star=6*sqrt(U_H/U_N5to6_ionization);

  C=(2*l+1)*pow(2, 4.*n_star-4.)/(n_star*tgamma(n_star+l_star+1)*tgamma(n_star-l_star));
  buf1=pow(U_N5to6_ionization/U_H, (6.*n_star-1.)/2.)/pow(fabs(Eionization/E_atomic_norm), 2*n_star-1);
  buf2=2*pow(U_N5to6_ionization/U_H, 1.5)/(3*fabs(Eionization/E_atomic_norm));
  *W56=4*omega_a_norm*C*buf1*exp(-buf2);
  *W67=0;
}

//==========================

// Krypton 
#define U_Kr0to1_ionization 13.999 // eV
#define U_Kr8to9_ionization 230.39 // eV

void ionizationPhysicsDetails_Krypton(double *W89, double *W910, double zz, double rr)
{
  /* 
     This function sets the values of W89 and/or W910 starting from the amplitude of the ionizing electric field specified in E_laserField_ionization.
     N.B.1. E_laserField_ionization is normalized with E0, so it might be required to denormalize it using E0_amplitudeNormalizationField_ionization (if needed)
  */
  int ii, jj;
  double buf_zz, buf_rr, w00, w01, w10, w11;

  double omega_a_norm, E_atomic_norm, C;
  double n_star, n0_star, l, l_star, Eionization;
  double buf1, buf2;

  if (zz<zmin_ionization || zz>zmax_ionization-dzl_ionization || rr>rmax_ionization-drl)
    {
      *W89=0;
      *W910=0;
      return;
    }

  buf_zz=(zz-zmin_ionization)/dzl_ionization;
  buf_rr=rr/drl_ionization;
  
  ii=(int)buf_zz;
  jj=(int)buf_rr;
  
  buf_zz-=(double)ii;
  buf_rr-=(double)jj;
  
  w00=(1-buf_zz)*(1-buf_rr);
  w01=(1-buf_zz)*buf_rr;
  w10=buf_zz*(1-buf_rr);
  w11=buf_zz*buf_rr;

  Eionization=w00*E_laserField_ionization[ii][jj]+w01*E_laserField_ionization[ii][jj+1]+w10*E_laserField_ionization[ii+1][jj]+w11*E_laserField_ionization[ii+1][jj+1];
  
  n0_star=sqrt(U_H/U_Kr0to1_ionization);
  l_star=n0_star-1;
  l=1; // electron #8 is in the p-orbitals XXXXXXXXXXXXXXXXXXXXXXXXXXXX check

  omega_a_norm=pow(ALPHA_FINE_STRUCTURE_CONSTANT,3)/(kp_cgsUnits_ionization*R0_cgs);
  E_atomic_norm=pow(ALPHA_FINE_STRUCTURE_CONSTANT,4)/(kp_cgsUnits_ionization*R0_cgs);
    
  // 8->9
  n_star=9*sqrt(U_H/U_Kr8to9_ionization);

  C=(2*l+1)*pow(2, 4.*n_star-4.)/(n_star*tgamma(n_star+l_star+1)*tgamma(n_star-l_star));
  buf1=pow(U_Kr8to9_ionization/U_H, (6.*n_star-1.)/2.)/pow(fabs(Eionization/E_atomic_norm), 2*n_star-1);
  buf2=2*pow(U_Kr8to9_ionization/U_H, 1.5)/(3*fabs(Eionization/E_atomic_norm));
  *W89=4*omega_a_norm*C*buf1*exp(-buf2);
  *W910=0;
}

//==========================

// Argon
#define U_Ar0to1_ionization 15.759 // eV
#define U_Ar8to9_ionization 422.44 // eV
#define U_Ar9to10_ionization 478.68 // eV
#define U_Ar10to11_ionization 538.95 // eV

void ionizationPhysicsDetails_Argon8to9_10(double *W8_9, double *W9_10, double zz, double rr)
{
  /* 
     This function sets the values of W8_9 and W9_10 starting from the amplitude of the ionizing electric field specified in E_laserField_ionization.
     N.B.1. E_laserField_ionization is normalized with E0, so it might be required to denormalize it using E0_amplitudeNormalizationField_ionization (if needed)
  */
  int ii, jj;
  double buf_zz, buf_rr, w00, w01, w10, w11;

  double omega_a_norm, E_atomic_norm, C;
  double n_star, n0_star, l, l_star, Eionization;
  double buf1, buf2;

  if (zz<zmin_ionization || zz>zmax_ionization-dzl_ionization || rr>rmax_ionization-drl)
    {
      *W8_9=0;
      *W9_10=0;
      return;
    }

  buf_zz=(zz-zmin_ionization)/dzl_ionization;
  buf_rr=rr/drl_ionization;
  
  ii=(int)buf_zz;
  jj=(int)buf_rr;
  
  buf_zz-=(double)ii;
  buf_rr-=(double)jj;
  
  w00=(1-buf_zz)*(1-buf_rr);
  w01=(1-buf_zz)*buf_rr;
  w10=buf_zz*(1-buf_rr);
  w11=buf_zz*buf_rr;

  Eionization=w00*E_laserField_ionization[ii][jj]+w01*E_laserField_ionization[ii][jj+1]+w10*E_laserField_ionization[ii+1][jj]+w11*E_laserField_ionization[ii+1][jj+1];
  
  n0_star=sqrt(U_H/U_Ar0to1_ionization);
  l_star=n0_star-1;
  l=1; 

  omega_a_norm=pow(ALPHA_FINE_STRUCTURE_CONSTANT,3)/(kp_cgsUnits_ionization*R0_cgs);
  E_atomic_norm=pow(ALPHA_FINE_STRUCTURE_CONSTANT,4)/(kp_cgsUnits_ionization*R0_cgs);
    
  // 8->9
  n_star=9*sqrt(U_H/U_Ar8to9_ionization);

  C=(2*l+1)*pow(2, 4.*n_star-4.)/(n_star*tgamma(n_star+l_star+1)*tgamma(n_star-l_star));
  buf1=pow(U_Ar8to9_ionization/U_H, (6.*n_star-1.)/2.)/pow(fabs(Eionization/E_atomic_norm), 2*n_star-1);
  buf2=2*pow(U_Ar8to9_ionization/U_H, 1.5)/(3*fabs(Eionization/E_atomic_norm));
  *W8_9=4*omega_a_norm*C*buf1*exp(-buf2);

  // 9->10
  n_star=10*sqrt(U_H/U_Ar9to10_ionization);

  C=(2*l+1)*pow(2, 4.*n_star-4.)/(n_star*tgamma(n_star+l_star+1)*tgamma(n_star-l_star));
  buf1=pow(U_Ar9to10_ionization/U_H, (6.*n_star-1.)/2.)/pow(fabs(Eionization/E_atomic_norm), 2*n_star-1);
  buf2=2*pow(U_Ar9to10_ionization/U_H, 1.5)/(3*fabs(Eionization/E_atomic_norm));
  *W9_10=4*omega_a_norm*C*buf1*exp(-buf2);
}

void ionizationPhysicsDetails_Argon9to10_11(double *W9_10, double *W10_11, double zz, double rr)
{
  /* 
     This function sets the values of W9_10 and W10_11 starting from the amplitude of the ionizing electric field specified in E_laserField_ionization.
     N.B.1. E_laserField_ionization is normalized with E0, so it might be required to denormalize it using E0_amplitudeNormalizationField_ionization (if needed)
  */
  int ii, jj;
  double buf_zz, buf_rr, w00, w01, w10, w11;

  double omega_a_norm, E_atomic_norm, C;
  double n_star, n0_star, l, l_star, Eionization;
  double buf1, buf2;

  if (zz<zmin_ionization || zz>zmax_ionization-dzl_ionization || rr>rmax_ionization-drl)
    {
      *W9_10=0;
      *W10_11=0;
      return;
    }

  buf_zz=(zz-zmin_ionization)/dzl_ionization;
  buf_rr=rr/drl_ionization;
  
  ii=(int)buf_zz;
  jj=(int)buf_rr;
  
  buf_zz-=(double)ii;
  buf_rr-=(double)jj;
  
  w00=(1-buf_zz)*(1-buf_rr);
  w01=(1-buf_zz)*buf_rr;
  w10=buf_zz*(1-buf_rr);
  w11=buf_zz*buf_rr;

  Eionization=w00*E_laserField_ionization[ii][jj]+w01*E_laserField_ionization[ii][jj+1]+w10*E_laserField_ionization[ii+1][jj]+w11*E_laserField_ionization[ii+1][jj+1];
  
  n0_star=sqrt(U_H/U_Ar0to1_ionization);
  l_star=n0_star-1;
  l=1; 

  omega_a_norm=pow(ALPHA_FINE_STRUCTURE_CONSTANT,3)/(kp_cgsUnits_ionization*R0_cgs);
  E_atomic_norm=pow(ALPHA_FINE_STRUCTURE_CONSTANT,4)/(kp_cgsUnits_ionization*R0_cgs);
    
  // 9->10
  n_star=10*sqrt(U_H/U_Ar9to10_ionization);

  C=(2*l+1)*pow(2, 4.*n_star-4.)/(n_star*tgamma(n_star+l_star+1)*tgamma(n_star-l_star));
  buf1=pow(U_Ar9to10_ionization/U_H, (6.*n_star-1.)/2.)/pow(fabs(Eionization/E_atomic_norm), 2*n_star-1);
  buf2=2*pow(U_Ar9to10_ionization/U_H, 1.5)/(3*fabs(Eionization/E_atomic_norm));
  *W9_10=4*omega_a_norm*C*buf1*exp(-buf2);

  // 10->11
  n_star=11*sqrt(U_H/U_Ar10to11_ionization);

  C=(2*l+1)*pow(2, 4.*n_star-4.)/(n_star*tgamma(n_star+l_star+1)*tgamma(n_star-l_star));
  buf1=pow(U_Ar10to11_ionization/U_H, (6.*n_star-1.)/2.)/pow(fabs(Eionization/E_atomic_norm), 2*n_star-1);
  buf2=2*pow(U_Ar10to11_ionization/U_H, 1.5)/(3*fabs(Eionization/E_atomic_norm));
  *W10_11=4*omega_a_norm*C*buf1*exp(-buf2);
}
