// MTpkg - Luciano Ristori 
// 
//-----------------------------------------------
// Mytunes project started on November 13, 2005
//
// Mytunes is my first attempt to study the statistical properties of music signals
//
// First goal is to read a .wav file and construct the pseudo-fourier spectrum
// of the auto-correlation function.
//
// Pseudo-fourier is obtained projecting on damped cosine waves. Time constant
// is proportional to period (constant fractional width) and frequencies are
// log distributed (intervals are constant fraction of semitone).
//
// We start by building one utility function to open a wav file and to provide
// one sample at a time till the EOF.
//
// Version 0.0 - Start project - First trials
// 15 Nov 2005 - getsample(int *sample, FILE *infile) is done and tested 
//
// Now attack problem of autocorrelation function.
// 20 Nov 2005 - autocorrelation seems to work.
//
// Now attack fourier
// 23 Nov 2005 - autofourier seems to work
//
// Develop graphic spectrum representation in .ps
//
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "PSpkg.h"
#include "MTpkg.h"

// global constants

// lowest frequency of spectrum
#define MINFREQ (27.5)

// corresponding "note" (0.=C, 1.=C#, 2.=D, 3.=D# ...)
#define MINNOTE (9.) // it's an A

// number of frequencies in spectrum
#define NFREQS (288)

// number of frequency steps in one semitones
#define TONESTEPS (3)


//**********************************************

// Flip four bytes in int - (big endian <-> little endian)

int flip4bytes(int input){
  int result,i;
  result=input & 0xFF;
  for(i=0;i<3;++i){
    input = input >> 8;
    result = (result << 8) | input & 0xFF;
  }
  return result;
}

//**********************************************

// Flip two bytes in short int - (big endian <-> little endian)

short int flip2bytes(short int input){
  short int result;
  result = input & 0xFF;
  input = input >> 8;
  result = (result << 8) | input & 0xFF;
  return result;
}

//**********************************************

// Provides next sample from file
// Assumes file is already open
//
// Return values: 0 = OK
//                1 = EOF
//                2 = Error

int getsample(int *sample, FILE *infile) {
#define MAXDATA 16384
#define RIFF (1179011410)
#define WAVE (1163280727)
#define FMT  (544501094)
#define DATA (1635017060)

    static short int isfirst = 1;
    static short int data[MAXDATA];
    static int datapointer = MAXDATA;
    static size_t nread;

    static long dataStart = 0; // file position of PCM data
    int idata[256];
    int ierr = 0;

    if (isfirst) {
        isfirst = 0;

        // --- Read and validate RIFF/WAVE ---
        nread = fread(idata, 4, 3, infile);
        if (nread != 3) { fprintf(stderr,"File too short\n"); return 2; }
        if (idata[0] != RIFF) { fprintf(stderr,"Not a RIFF file\n"); return 2; }
        if (idata[2] != WAVE) { fprintf(stderr,"WAVE chunk not found\n"); return 2; }

        // --- Search for "data" marker anywhere in first 1 kB ---
        rewind(infile);
        unsigned char buf[1024];
        size_t bytes = fread(buf,1,sizeof(buf),infile);
        int found = 0;
        for (size_t i=0;i+4<bytes;i++){
            if (buf[i]=='d' && buf[i+1]=='a' && buf[i+2]=='t' && buf[i+3]=='a'){
                found = 1;
                dataStart = i + 8; // skip "data" + length
                break;
            }
        }
        if (!found){
            fprintf(stderr,"DATA chunk not found\n");
            return 2;
        }

        // seek to start of data
        fseek(infile, dataStart, SEEK_SET);
    }

    // --- refill data buffer if needed ---
    if (datapointer >= (int)nread - 1) {
        nread = fread(data, 2, MAXDATA, infile);
        if (nread == 0) return feof(infile) ? 1 : 2;
        datapointer = 0;
    }

    // sum left + right (16-bit little-endian stereo)
    *sample = data[datapointer++];
    *sample += data[datapointer++];
    return 0;
}


//*************************************

int autofourier(int sample,int opt, double **pfourier, double **pfourierror){
  // opt = 0: feed next sample
  //     = 1: spew results

#define PI (3.1415926536)

// length of circular buffer for cos
#define NCOSBUFF (4096)

// mask to ensure circularity of cos buffer (to be ANDed to index)
#define MASKCOSBUFF (0xFFF)

// length of buffer for log
#define NLOGBUFF (4096)

// length of circular buffer for raw sound samples
#define NBUFF (65536)

// mask to ensure circularity of data buffer (to be ANDed to index)
#define MASKBUFF (0xFFFF)

// damping constant of cos wave
#define DAMPCONST (100)

// number of damping constants for integration
#define NDAMP (5)

// number of samples for Montecarlo integration
#define NINTEGRATE (4000)

// time constant for running time average
// in milliseconds. Used to calculate alfa and beta

#define TIMECONST (100)

  // constants for running time average
  static double alfa, beta;

  static float minrand;

  // spectrum of frequencies  
  static double omega[NFREQS];
  static double fourier[NFREQS];
  static double fourierror[NFREQS];
  static double fouriersquares[NFREQS];
  static double tinytonestep;

  // circular buffer for holding data
  static int databuffer[NBUFF];
  static int dataindex = 0;

  // circular buffer for cos
  static double cosbuffer[NCOSBUFF];

  // buffer for log
  static double logbuffer[NLOGBUFF];

  // flags first entry time for initialization

  static int firsttime = 1;

  // flags enough data received to start processing

  static int lead = 0;

  // local volatile variables

  int i, k, indtime, ifreq;
  double temp, x, tnorm;

  // Initialize

  if(firsttime){
    firsttime = 0;

    // gives pointers to buffers back to calling programs
    *pfourier=&fourier[0];
    *pfourierror=&fourierror[0];

    tinytonestep = exp(log(2.)/12./TONESTEPS);//relative freq step

    // compute frequencies of spectrum and clear fourier

    omega[0] = MINFREQ*2.*PI;// minimum omega
    fourier[0] = 0.;

    for(i=1;i<NFREQS;++i){
      omega[i] = omega[i-1]*tinytonestep; 
      fourier[i] = 0.;
      fouriersquares[i] = 0.;
    }

    // clear data buffer
    for(i=0;i<NBUFF;++i)databuffer[i]=0;

    // initialize cos buffer
    for(i=0; i<NCOSBUFF; ++i)
      cosbuffer[i] = cos(((double)i)/NCOSBUFF*2.*PI);

    // initialize log buffer
    for(i=0; i<NLOGBUFF; ++i)
      logbuffer[i] = log(((double)(i+1))/NLOGBUFF);
    // compute alfa and beta to obtain the correct
    // time constant for running average

    beta = 1./(44100*TIMECONST/1000.*NINTEGRATE/NFREQS);
    alfa = 1.- beta;

    // minimum value of random throw for time
    minrand = exp(-NDAMP);
  }

  // end initialize

  // normal flow comes here

  if(opt == 1){

    // spew results
 
    for(i=0;i<NFREQS;++i)
      fourierror[i]=(fouriersquares[i] - fourier[i]*fourier[i])*beta;    
    return 0;
  }


  // store next sample in buffer

  ++dataindex;
  dataindex &=MASKBUFF;// ensure circularity

  databuffer[dataindex] = sample;

  // if not enough samples skip processing

  //  if(lead < NBUFF){
  //    ++lead;
  //    return 0;
  //  }

  // if enough samples compute correlations

  for(k=0;k<NINTEGRATE;++k){
    
    //    printf(" %d %d\n", dataindex, k);

    // pick a random time in normalized damped cos wave
    x = (double)rand()/(double)RAND_MAX;
    if(x < minrand) continue;// out of range
    //tnorm = -log(x)*DAMPCONST;// generate exponential distribution
    tnorm = -DAMPCONST*logbuffer[(int)(x*NLOGBUFF-0.005)];
    
    // pick a random frequency in spectrum
    ifreq = rand() % NFREQS;

    // convert normalized time to real time (in samples)
    indtime = 44100.*tnorm/omega[ifreq];
    if(indtime >= NBUFF) continue; // data is lost (too old)

    // fill fourier transforms with autocorrelation * cos wave
    temp = databuffer[dataindex]*databuffer[(dataindex-indtime)&MASKBUFF];
    //temp = temp*cos(tnorm);
    temp = temp*cosbuffer[((int)(tnorm*NCOSBUFF/2./PI)) & MASKCOSBUFF];
    fourier[ifreq] = fourier[ifreq]*alfa + temp*beta;// apply time constant to average
    fouriersquares[ifreq] = fouriersquares[ifreq]*alfa + temp*temp*beta;// apply time constant to average of squares   
    //    printf("%d %d %d %f\n",dataindex, indtime, ifreq, temp);
  }
}

//
// Plot power spectrum with errors on one page starting with mintone
//

void plotfourier(int isample,double seconds, double avgvol, double *fourier, double *fourierror, int mintone, int ntones){

// decay time of maximum in seconds
#define SCALETIMECONSTANT (5.)

// high-left corner coordinates
#define XHILEFT (50)
#define YHILEFT (750)

// length of sides
#define YLEN (700)
#define XLEN (500)

// radius of points
#define POINTRADIUS (3.)

// page counter

  static int ipage = 1;

  static double lastseconds = 0.;
  
  double sigma[NFREQS];
  static double minvalue=0., maxvalue=0., newmaxvalue;
  double scale, xpoint, xlow, xhi, ystep, yleft;
  int itone, inote, ikey;
  char stemp[128];

  // establish range of values to set the scale
  // and compute square root of sigmas2 for all frequencies
  minvalue=-10.;
  newmaxvalue=+10.;
  for(itone=mintone;itone<mintone+ntones;++itone){
    sigma[itone] = sqrt(fourierror[itone]);
    if(newmaxvalue < fourier[itone]+ sigma[itone]) newmaxvalue = fourier[itone]+ sigma[itone];
    if(minvalue > fourier[itone]- sigma[itone]) minvalue = fourier[itone]- sigma[itone];
  }
 
  // set limits based on average volume

  minvalue = 0.;
  //  maxvalue = avgvol*NORMFACTOR;
  if(newmaxvalue > maxvalue) maxvalue = newmaxvalue;
  else {
    maxvalue = newmaxvalue + (maxvalue - newmaxvalue)*exp(-(seconds-lastseconds)/SCALETIMECONSTANT);
  }
  lastseconds=seconds;

  // compute scale factor
  scale = XLEN/(maxvalue-minvalue);

  // compute step in Y
  ystep=YLEN/(float)(ntones+1);

  // draw piano keys
  yleft=10000.;
  for(inote=0; yleft-TONESTEPS*ystep > YHILEFT-YLEN;++inote){
    if(inote >= MINNOTE+(double)mintone/TONESTEPS){
      yleft = YHILEFT-(inote-MINNOTE-(double)mintone/TONESTEPS-0.5)*ystep*TONESTEPS-ystep;
      ikey=inote % 12;
      if(ikey == 0 ) PSline(XHILEFT,yleft,XHILEFT+XLEN,yleft,1.,0.,0.,1.);
      if(ikey == 5) PSline(XHILEFT,yleft,XHILEFT+XLEN,yleft,1.,0.8,0.8,0.8);
      if(ikey == 1 || ikey == 3 || ikey == 6 || ikey == 8 || ikey == 10)
	PSrectangleFill(XHILEFT,yleft,XHILEFT+XLEN,yleft-TONESTEPS*ystep,0.8,0.8,0.8);
    }
  }

 // draw plot frame
  PSline(XHILEFT,YHILEFT,XHILEFT+XLEN,YHILEFT,2.,0.,0.,0.);
  PSline(XHILEFT,YHILEFT,XHILEFT,YHILEFT-YLEN,2.,0.,0.,0.);
  PSline(XHILEFT,YHILEFT-YLEN,XHILEFT+XLEN,YHILEFT-YLEN,2.,0.,0.,0.);
  PSline(XHILEFT+XLEN,YHILEFT,XHILEFT+XLEN,YHILEFT-YLEN,2.,0.,0.,0.);

  // create heading
  sprintf(stemp,"%.3f s - sample %d",(double)isample/44100.,isample);

  // draw heading
  PSsetfont("Times-Roman",12.);
  PSmoveto(XHILEFT,YHILEFT+10.);
  PSshowleft(stemp);

  // create page number
  sprintf(stemp,"page %d",ipage);

  // draw page number
  PSsetfont("Times-Roman",12.);
  PSmoveto(XHILEFT+XLEN,YHILEFT+10.);
  PSshowright(stemp);


  // draw zero line
  //  PSline(XHILEFT-minvalue*scale,YHILEFT,XHILEFT-minvalue*scale,YHILEFT-YLEN,1.,0.,0.,0.);

   // draw all points with error bars
  for(itone=mintone;itone<mintone+ntones;++itone){
    xpoint = (fourier[itone]-minvalue)*scale+XHILEFT;
    xlow = (fourier[itone]-sigma[itone]-minvalue)*scale+XHILEFT;
    if(xlow < XHILEFT) xlow = XHILEFT;
    xhi = (fourier[itone]+sigma[itone]-minvalue)*scale+XHILEFT;
    if(xhi > XHILEFT+XLEN) xhi = XHILEFT+XLEN;
    if(xlow < XHILEFT+XLEN && xhi > XHILEFT)
      PSline(xlow,YHILEFT-(itone-mintone+1)*ystep,xhi,YHILEFT-(itone-mintone+1)*ystep,1.,0.,0.,0.);
    if((xpoint >= XHILEFT) && (xpoint <= XHILEFT+XLEN))
      PScircleFill(xpoint,YHILEFT-(itone-mintone+1)*ystep,POINTRADIUS,0.,0.,0.);
  }

  PSshowpage();
  ++ipage;
  return;
}

//
// compute average volume
//

double volume(int sample){

// time constant for running average in milliseconds
#define VOLTIMECONST (1000.)

  static int firsttime = 1;
  static double alfa, beta, runningaverage=0.;
  int nwrite;

  //
  // initialize
  //
  if(firsttime){ 
    firsttime = 0;
    beta = 1./(44100*VOLTIMECONST/1000.);
    alfa = 1.- beta;
  }
  //
  // end initialize
  //
  runningaverage = runningaverage*alfa + beta*sample*sample;
  return runningaverage;
}


//
// write pseudo-fourier spectrum to stdout
//

int writefourier(int isample, double seconds, double avgvol, double *fourier, double *fourierror){

  static int isfirst = 1;
  static int version = 1;
  static double minfreq = MINFREQ, minnote = MINNOTE;
  static int nfreqs = NFREQS, tonesteps = TONESTEPS;

  size_t nwrite;

  if(isfirst == 1){
    //
    // initialization: executed only first time called
    //
    isfirst = 0;

    // write header
    // file format version
    nwrite=fwrite(&version,sizeof(version),1,stdout);
    if(nwrite != 1){
      fprintf(stderr,"Failed writing header:version\n");
      return 1;
    }
    // minimum frequency
    nwrite=fwrite(&minfreq,sizeof(minfreq),1,stdout);
    if(nwrite != 1){
      fprintf(stderr,"Failed writing header:minfreq\n");
      return 1;
    }
    // minimum note
    nwrite=fwrite(&minnote,sizeof(minnote),1,stdout);
    if(nwrite != 1){
      fprintf(stderr,"Failed writing header:minnote\n");
      return 1;
    }
    // number of frequencies in spectrum
    nwrite=fwrite(&nfreqs,sizeof(nfreqs),1,stdout);
    if(nwrite != 1){
      fprintf(stderr,"Failed writing header:nfreqs\n");
      return 1;
    }
    // number of steps in one semitone
    nwrite=fwrite(&tonesteps,sizeof(tonesteps),1,stdout);
    if(nwrite != 1){
      fprintf(stderr,"Failed writing header:tonesteps\n");
      return 1;
    }
  }

  // body of function
  // executed all times

  nwrite=fwrite(&isample,sizeof(int),1,stdout);
    if(nwrite != 1){
      fprintf(stderr, "Failed writing at %.3f s, nwrite = %zu\n", seconds, nwrite);
      return 1;
    }
  nwrite=fwrite(&seconds,sizeof(double),1,stdout);
    if(nwrite != 1){
      fprintf(stderr, "Failed writing at %.3f s, nwrite = %zu\n", seconds, nwrite);
      return 1;
    }
  nwrite=fwrite(&avgvol,sizeof(double),1,stdout);
    if(nwrite != 1){
       fprintf(stderr, "Failed writing at %.3f s, nwrite = %zu\n", seconds, nwrite);
      return 1;
    }
  nwrite=fwrite(fourier,sizeof(double),nfreqs,stdout);
    if(nwrite != nfreqs){
      fprintf(stderr, "Failed writing at %.3f s, nwrite = %zu\n", seconds, nwrite);
      return 1;
    }
  nwrite=fwrite(fourierror,sizeof(double),nfreqs,stdout);
    if(nwrite != nfreqs){
      fprintf(stderr, "Failed writing at %.3f s, nwrite = %zu\n", seconds, nwrite);
      return 1;
    }
    return 0;
}

// check contents of header of input file

int checkheader(){

  int nread;
  int ierr=0; 
  int correctversion = 1;
  int version, nfreqs, tonesteps;
  double minfreq, minnote;

  // check version number

  nread = fread(&version,sizeof(version),1,stdin);
  if(nread != 1){
    fprintf(stderr,"Failed reading header:version\n");
    return 1;
  }
  if(version != correctversion){
    fprintf(stderr,"File format version is %d. Should be %d\n",
	    version, correctversion);
    ierr=1;
  }

  // check minimum frequency

  nread = fread(&minfreq,sizeof(minfreq),1,stdin);
  if(nread != 1){
    fprintf(stderr,"Failed reading header:minfreq\n");
    return 1;
  }
  if(minfreq != MINFREQ){
    fprintf(stderr,"minfreq is %.2f. Should be %.2f\n",
	    minfreq, MINFREQ);
    ierr=1;
  }

  // check minimum note

  nread = fread(&minnote,sizeof(minnote),1,stdin);
  if(nread != 1){
    fprintf(stderr,"Failed reading header:minnote\n");
    return 1;
  }
  if(minnote != MINNOTE){
    fprintf(stderr,"minnote is %.2f. Should be %.2f\n",
	    minnote, MINNOTE);
    ierr=1;
  }

  // check number of frequencies

  nread = fread(&nfreqs,sizeof(nfreqs),1,stdin);
  if(nread != 1){
    fprintf(stderr,"Failed reading header:nfreqs\n");
    return 1;
  }
  if(nfreqs != NFREQS){
    fprintf(stderr,"nfreqs is %d. Should be %d\n",
	    nfreqs, NFREQS);
    ierr=1;
  }

  // check number of steps in one semitone

  nread = fread(&tonesteps,sizeof(tonesteps),1,stdin);
  if(nread != 1){
    fprintf(stderr,"Failed reading header:tonesteps\n");
    return 1;
  }
  if(tonesteps != TONESTEPS){
    fprintf(stderr,"tonesteps is %d. Should be %d\n",
	    tonesteps, TONESTEPS);
    ierr=1;
  }

  fprintf(stderr,"version %d\n minfreq %.2f\n minnote %.2f\n nfreqs %d\n tonesteps %d\n",
	  version,minfreq,minnote,nfreqs,tonesteps);

  if(ierr != 0) return 1;
  return 0;
}
   
// read one record from input file

int readfourier(int *isample, double *seconds, double *avgvol, 
		double *fourier, double *fourierror){

  int nread;

  // check for EOF
  
  if(feof(stdin)){
    fprintf(stderr," EOF on input\n");
    return 1;
  }

  // start reading

  nread = fread(isample,sizeof(int),1,stdin);

  if(nread != 1){
    fprintf(stderr,"Failed to read input record:isample\n");
    return 2;
  }
  nread = fread(seconds,sizeof(double),1,stdin);
  if(nread != 1){
    fprintf(stderr,"Failed to read input record:seconds\n");
    return 2;
  }
  nread = fread(avgvol,sizeof(double),1,stdin);
  if(nread != 1){
    fprintf(stderr,"Failed to read input record:avgvol\n");
    return 2;
  }


  nread = fread(fourier,sizeof(double),NFREQS,stdin);
  if(nread != NFREQS){
    fprintf(stderr,"Failed to read fourier data: nread = %d. Should be %d\n", nread,NFREQS);
    return 2;
  }
  nread = fread(fourierror,sizeof(double),NFREQS,stdin);
  if(nread != NFREQS){
    fprintf(stderr,"Failed to read fourierror data: nread = %d. Should be %d\n", nread,NFREQS);
    return 2;
  }

  return 0;
}
