void dumpFieldASCII_generalGrid(double **field, double ZMIN__, double DZ__, double DR__, int NZ__, int NR__, char *fileName, double sss)
{
  /*
    must be a "real" field, defined on a NZ x NR grid
  */
  
  int i, j;
  FILE *f;
  char buffer[200];

  //printf("%i %i\n", NZ, NR);
  //printf("%e %e %e\n", ZMIN, ZMIN+(NZ-1)*DZ, DR);
  
  sprintf(buffer, "%s/%s_%.1f", DIRECTORY_OUTPUT, fileName, sss);
  f=fopen(buffer, "w");
  fprintf(f, "%i %i 1\n", NZ__, NR__);
  fprintf(f, "%e 0 %f %f\n", ZMIN__, ZMIN__+(NZ__-1)*DZ__, (NR__-1)*DR__);
  
  for (i=0; i<NZ__; i++)
    for (j=0; j<NR__; j++)
      fprintf(f, "%e %e %e\n", ZMIN__+i*DZ__, j*DR__, field[i][j]);
  fclose(f);
}

void dumpFieldASCII2_generalGrid(double **field, double ZMIN__, double DZ__, double DR__, int NZ__, int NR__, double simmetry, char *fileName, double sss)
{
  /*
    must be a "real" field, defined on a NZ x NR grid
  */
  
  int i, j;
  FILE *f;
  char buffer[200];

  //printf("%i %i\n", NZ, NR);
  //printf("%e %e %e\n", ZMIN, ZMIN+(NZ-1)*DZ, DR);
  
  sprintf(buffer, "%s/%s_%.1f", DIRECTORY_OUTPUT, fileName, sss);
  f=fopen(buffer, "w");
  fprintf(f, "%i %i 1\n", NZ__, 2*NR__-1);
  fprintf(f, "%e %e %e %e\n", ZMIN__, -(NR__-1)*DR__, ZMIN__+(NZ__-1)*DZ__, (NR__-1)*DR__);
  
  for (i=0; i<NZ__; i++)
    {
      for (j=NR__-1; j>0; j--)
	fprintf(f, "%e %e %e\n", ZMIN__+i*DZ__, -j*DR__, simmetry*field[i][j]);

      for (j=0; j<NR__; j++)
	fprintf(f, "%e %e %e\n", ZMIN__+i*DZ__, j*DR__, field[i][j]);
    }
  
  fclose(f);
}

void dumpIonized(char *fileName, double sss)
{ 
  int n;
  FILE *f;
  char buffer[200];
  
  sprintf(buffer, "%s/%s_%.1f", DIRECTORY_OUTPUT, fileName, sss);
  f=fopen(buffer, "w");

  for (n=0; n<Nelectron_ionization; n++)
    fprintf(f, "%e %e %e %e %e %e %e\n", x_ionization[n], y_ionization[n], z_ionization[n], ux_ionization[n], uy_ionization[n], uz_ionization[n], w_ionization[n]); // x,y,z particles
  
  fclose(f);
}

void dumpIonizedStatus(char *fileName, double sss)
{ 
  int n;
  FILE *f;
  char buffer[200];
  
  sprintf(buffer, "%s/%s_%.1f", DIRECTORY_OUTPUT, fileName, sss);
  f=fopen(buffer, "w");

  for (n=0; n<Nelectron_ionization; n++)
    fprintf(f, "%e %e %e %e %e %e %e %e %i %i\n", x_ionization[n], y_ionization[n], z_ionization[n], ux_ionization[n], uy_ionization[n], uz_ionization[n], w_ionization[n], 1., 1, -1);
  
  fclose(f);
}

void dumpCompleteIonizedStatus(char *fileName, double sss)
{ 
  int n;
  FILE *f;
  char buffer[200];
  
  sprintf(buffer, "%s/%s_%.1f", DIRECTORY_OUTPUT, fileName, sss);
  f=fopen(buffer, "w");

  for (n=0; n<Nelectron_ionization; n++)
    fprintf(f, "%e %e %e %e %e %e %e %e %e %e %i\n", x_ionization[n], y_ionization[n], z_ionization[n], ux_ionization[n], uy_ionization[n], uz_ionization[n], x0_ionization[n], y0_ionization[n], z0_ionization[n], w_ionization[n], label_ionization[n]); // x,y,z particles
  
  fclose(f);
}

void dumpCompleteIonizedStatus2(char *fileName, double sss)
{ 
  int n;
  FILE *f;
  char buffer[200];
  
  sprintf(buffer, "%s/%s_%.1f", DIRECTORY_OUTPUT, fileName, sss);
  f=fopen(buffer, "w");

  for (n=0; n<Nelectron_ionization; n++)
    fprintf(f, "%e %e %e %e %e %e %e %e %e %e %e %i\n", x_ionization[n], y_ionization[n], z_ionization[n], ux_ionization[n], uy_ionization[n], uz_ionization[n], x0_ionization[n], y0_ionization[n], z0_ionization[n], zeta0_ionization[n], w_ionization[n], label_ionization[n]); // x,y,z particles
  
  fclose(f);
}
