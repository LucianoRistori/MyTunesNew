// MYTUNES - Luciano Ristori 


# include "MTpkg.c"


int main(){

  // This is a test program for getsample and autocorr

  int i,k;
  int sample, retcode;
  int vmin = 0, vmax = 0;
  double fsample, vmean=0., vrms=0.;
  double *pfourier, *pfourierror;
  double xmax, sigma,seconds,avgvol;;
  int kchar,nchar;
  int kmin=0, kmax=48;
  int mintone, ntones;

  for(i=0; ;++i){
    seconds = i/44100.;
    retcode=getsample(&sample,stdin);
    if(retcode != 0){
      if(retcode==1)fprintf(stderr,"File ends after %d seconds\n", i/44100);
      if(retcode==2)fprintf(stderr,"Read error after %d seconds\n", i/44100);
      fprintf(stderr,"vmin  = %d\n", vmin);
      fprintf(stderr,"vmax  = %d\n", vmax);
      fprintf(stderr,"vmean = %.0f\n", vmean/i);
      fprintf(stderr,"vrms  = %.0f\n", sqrt(vrms/i));
      return 0;
    }
    if(sample < vmin ) vmin=sample;
    if(sample > vmax ) vmax=sample;
    fsample = sample;
    vmean += fsample;
    vrms += fsample*fsample;

    retcode = autofourier(sample,0, &pfourier, &pfourierror);
    avgvol=volume(sample);

    if((i % 4410)==0){
      retcode = autofourier(sample,1, &pfourier, &pfourierror);
      mintone = 0;
      ntones = NFREQS;
      fprintf(stderr,"plotfourier at %.3f s\n",seconds);
      plotfourier(i, seconds, avgvol, pfourier, pfourierror,mintone,ntones);

      if(seconds > 5.) return 0;
   
/*      printf(" %f secs\n",i/44100.);
      xmax =20000000.;
      printf(" %f\n", xmax);
      for(k=kmin;k<kmax;++k){
	sigma=sqrt(*(pfourierror+k));
	nchar =  20+100*(*(pfourier+k)-sigma)/xmax;
	printf(" %6d %+10.0f ",k,*(pfourier+k));
	for(kchar=0;kchar<nchar;++kchar)printf(" ");
	nchar = 20+100*(*(pfourier+k)+sigma)/xmax;
	for(;kchar<=nchar;++kchar)printf("x");
	printf("\n");
	}*/
    }
  }
}
