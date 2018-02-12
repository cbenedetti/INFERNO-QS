void dumpBeamASCII_clean(char *fileName, int indx0, int Np)
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
    {
#ifndef O_NOT_CLEAN
      if (particle_beam_active[indx0_+i]==NO) continue;
#endif
      fprintf(f, "%e %e %e %e %e %e %e %e %i %i\n", xpb[indx0_+i], ypb[indx0_+i], zpb[indx0_+i], uxb[indx0_+i], uyb[indx0_+i], uzb[indx0_+i], q[indx0_+i], me_mb[indx0_+i], particle_beam_active[indx0_+i], self_consistent[indx0_+i]);
    }
  fclose(f);
}
