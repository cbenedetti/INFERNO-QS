#include"bessel.c"

void filterLongitudinallyLaserEnvelope(int Npass)
{
  int i, j, k;
  static double *data_re=NULL, *data_im=NULL;
  
  if (data_re==NULL)
    {
      data_re=(double*)malloc(Nzl*sizeof(double));
      data_im=(double*)malloc(Nzl*sizeof(double));
    }

  for (j=0; j<Nrl; j++)
    {
      for (i=0; i<Nzl; i++)
	{
	  data_re[i]=creal(a[i][j]);
	  data_im[i]=cimag(a[i][j]);
	}
      
      for (k=0; k<Npass; k++)
	binomialFilter(data_re, Nzl, 0.5);
      binomialFilter(data_re, Nzl, 0.5*Npass+1); // compensator

      for (k=0; k<Npass; k++)
	binomialFilter(data_im, Nzl, 0.5);
      binomialFilter(data_im, Nzl, 0.5*Npass+1); // compensator
      
      for (i=0; i<Nzl; i++)
	a[i][j]=data_re[i]+I*data_im[i];
    }
}

double centroidPosition()
{
  int i, j;
  double w, z, zc, norm;
  
  zc=norm=0;
  for (i=0; i<Nz; i++)
    for (j=0; j<Nr; j++)
      {
	z=zmin+i*dz;

	w=a2_half[i][j]*j;	
	norm+=w;
	zc+=z*w;
      }
  zc/=norm;
  
  return zc;
}


double rmsRadius()
{
  int i, j;
  double w, r, r2, norm;
  
  r2=norm=0;
  for (i=0; i<Nz; i++)
    for (j=0; j<Nr; j++)
      {
	r=j*dr;

	w=a2_half[i][j]*j;	
	norm+=w;
	r2+=r*r*w;
      }
  r2/=norm;
  
  return sqrt(r2);
}

double quarticRadius()
{
  int i, j;
  double w, r, r4, norm;
  
  r4=norm=0;
  for (i=0; i<Nz; i++)
    for (j=0; j<Nr; j++)
      {
	r=j*dr;

	w=a2_half[i][j]*j;	
	norm+=w;
	r4+=r*r*r*r*w;
      }
  r4/=norm;
  
  return sqrt(sqrt(r4));
}

double jinc(double x)
{
  if (x<0.01)
    return 1-x*x/8.;
  else
    return 2*bessj(1, x)/x;
}

complex jinc_envelope_atFocus_simplified(double z, double r, double s)
{
  // kpW is the one from the equivalent Gaussian
  return a0*jinc(2.7455*r/kpW)*exp(-pow((z-kpz0)/kpL,2));
}

void absorbingBoundaryCondition(double rmin, double Tdamp_rmin, double Tdamp_rmax)
{
  int i, j, jmin;
  double Tdamp, epsi;

  jmin=rmin/drl;

  for (i=0; i<Nzl; i++)
    for (j=jmin; j<Nrl; j++)
      {
	epsi=(j-jmin)/(double)(Nrl-jmin-1);
	Tdamp=Tdamp_rmin*(1-epsi)+Tdamp_rmax*epsi;

	a[i][j]=a[i][j]*exp(-ds/Tdamp);
      }
}
