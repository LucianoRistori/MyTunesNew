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
#include <string.h>
#include <ctype.h>   // for isspace()


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

int getsample(int *sample, FILE *infile) {
#define MAXDATA 16384
    static short int isfirst = 1;
    static short int data[MAXDATA];
    static int datapointer = MAXDATA;
    static size_t nread = 0;
    static long dataStart = 0; // file position of PCM data

    if (isfirst) {
        isfirst = 0;
        // --- Read and validate RIFF/WAVE header ---
        unsigned char hdr[12];
        if (fread(hdr, 1, 12, infile) != 12) {
            fprintf(stderr, "File too short\n");
            return 2;
        }
        if (memcmp(hdr, "RIFF", 4) != 0) {
            fprintf(stderr, "Not a RIFF file\n");
            return 2;
        }
        if (memcmp(hdr + 8, "WAVE", 4) != 0) {
            fprintf(stderr, "WAVE chunk not found\n");
            return 2;
        }

        // --- Search for "data" marker anywhere in first 1 kB ---
        fseek(infile, 0, SEEK_SET);
        unsigned char buf[1024];
        size_t bytes = fread(buf, 1, sizeof(buf), infile);
        int found = 0;
        for (size_t i = 0; i + 4 < bytes; i++) {
            if (buf[i] == 'd' && buf[i+1] == 'a' && buf[i+2] == 't' && buf[i+3] == 'a') {
                found = 1;
                dataStart = i + 8; // skip "data" + length
                break;
            }
        }
        if (!found) {
            fprintf(stderr, "DATA chunk not found\n");
            return 2;
        }

        fseek(infile, dataStart, SEEK_SET);
        fprintf(stderr, "Reading PCM data starting at byte %ld\n", dataStart);
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

fprintf(stderr, "plotfourier() entered, PSfile=%p\n", (void*)PSfile);/////////////debug


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

// draw piano-style background for this plot area
  paintKeyboardRect(PSfile,
                    XHILEFT,
                    YHILEFT,
                    XLEN,
                    YLEN,
                    ntones,
                    mintone,
                    MINNOTE,   // from MTpkg.h
                    TONESTEPS);


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

  
  //PSline(50, 750, 50+500, 750-700, 1.5, 0.,0.,0.);  // should appear clearly


	// draw all points with error bars
	for (itone = mintone; itone < mintone + ntones; ++itone) {
		xpoint = (fourier[itone] - minvalue) * scale + XHILEFT;
		xlow   = (fourier[itone] - sigma[itone] - minvalue) * scale + XHILEFT;
		if (xlow < XHILEFT) xlow = XHILEFT;
		xhi    = (fourier[itone] + sigma[itone] - minvalue) * scale + XHILEFT;
		if (xhi > XHILEFT + XLEN) xhi = XHILEFT + XLEN;

		// align marker row with piano-style background
		double ybase = YHILEFT
					   - (itone - mintone + 1) * ystep
					   - (SEMITONE_SHIFT * TONESTEPS) * ystep;

		if (xlow < XHILEFT + XLEN && xhi > XHILEFT)
			PSline(xlow, ybase, xhi, ybase, 1., 0., 0., 0.);

		if (xpoint >= XHILEFT && xpoint <= XHILEFT + XLEN)
			PScircleFill(xpoint, ybase, POINTRADIUS, 0., 0., 0.);
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


//============================================================
// writefourier()
// Writes one pseudo-Fourier spectrum record in text format
// Header is written only once at the beginning of the stream
//============================================================

int writefourier(int isample, double seconds, double avgvol,
                 double *fourier, double *fourierror)
{
    static int isfirst = 1;
    if (isfirst) {
        printf("# version 1\n");
        printf("# minfreq %.6f\n", MINFREQ);
        printf("# minnote %.6f\n", MINNOTE);
        printf("# nfreqs %d\n", NFREQS);
        printf("# tonesteps %d\n", TONESTEPS);
        isfirst = 0;
    }

	//fprintf(stderr, "DEBUG1: writing frame %d (NFREQS=%d)\n", isample, NFREQS);


    // --- Record header line ---
    printf("%d %.6f %.6e\n", isample, seconds, avgvol);

    // --- Fourier amplitudes (10 per line) ---
    for (int i = 0; i < NFREQS; ++i) {
        printf(" %.6e", fourier[i]);
        if ((i + 1) % 10 == 0) printf("\n");
    }
    if (NFREQS % 10 != 0) printf("\n");

    // --- Fourier errors (10 per line) ---
    for (int i = 0; i < NFREQS; ++i) {
        printf(" %.6e", fourierror[i]);
        if ((i + 1) % 10 == 0) printf("\n");
    }
    if (NFREQS % 10 != 0) printf("\n");
    
    //fprintf(stderr, "DEBUG2: writing frame %d (NFREQS=%d)\n", isample, NFREQS);


    return 0;
}

// helper function

void writefourier_to_file(FILE *f, int i, double seconds, double avgvol,
                          double *pfourier, double *pfourierror)
{
    FILE *old = stdout;
    stdout = f;
    writefourier(i, seconds, avgvol, pfourier, pfourierror);
    fflush(f);
    stdout = old;
}



// read text header lines (for text-based spectrum.txt)

#include <stdio.h>
#include <string.h>
#include <ctype.h>

int checkheader(FILE *f)
{
    char line[256];
    double val;
    int nfreqs = 0;
    int found_version = 0, found_freq = 0, found_note = 0, found_steps = 0;

    fprintf(stderr, "Reading header:\n");

    while (fgets(line, sizeof(line), f)) {
        // Skip leading spaces
        char *p = line;
        while (isspace((unsigned char)*p)) p++;

        // Skip comment marker #
        if (*p == '#') p++;

        // Skip spaces again after #
        while (isspace((unsigned char)*p)) p++;

        // Stop if line looks like data (starts with digit)
        if (isdigit((unsigned char)*p) || *p == '\0')
            break;

        if (sscanf(p, "version %lf", &val) == 1) {
            fprintf(stderr, "  version = %.0f\n", val);
            found_version = 1;
        } else if (sscanf(p, "minfreq %lf", &val) == 1) {
            fprintf(stderr, "  minfreq = %.3f\n", val);
            found_freq = 1;
        } else if (sscanf(p, "minnote %lf", &val) == 1) {
            fprintf(stderr, "  minnote = %.3f\n", val);
            found_note = 1;
        } else if (sscanf(p, "nfreqs %d", &nfreqs) == 1) {
            fprintf(stderr, "  nfreqs = %d\n", nfreqs);
        } else if (sscanf(p, "tonesteps %lf", &val) == 1) {
            fprintf(stderr, "  tonesteps = %.0f\n", val);
            found_steps = 1;
            break;  // header done
        }
    }

    if (found_version && found_freq && found_note && found_steps) {
        fprintf(stderr, "Header OK.\n");
        return 0;
    } else {
        fprintf(stderr, "Header incomplete or malformed.\n");
        return 1;
    }
}



//============================================================
// readfourier()
// Reads one pseudo-Fourier spectrum record (text format)
// Expects the following structure:
//   # version 1
//   # minfreq ...
//   # minnote ...
//   # nfreqs ...
//   # tonesteps ...
//   isample seconds avgvol
//   <NFREQS Fourier values, 10 per line>
//   <NFREQS Fourier error values, 10 per line>
//============================================================

int readfourier(FILE *f, int *sample, double *seconds, double *avgvol,
                double *pfourier, double *pfourierror)
{
    // Read sample index, time (s), and average volume
    if (fscanf(f, "%d %lf %lf", sample, seconds, avgvol) != 3) {
        if (feof(f)) return 0;  // normal EOF
        return -1;              // read error
    }

    // Read the Fourier amplitudes
    for (int i = 0; i < NFREQS; ++i) {
        if (fscanf(f, "%lf", &pfourier[i]) != 1) {
            if (feof(f)) return 0;
            return -1;
        }
    }

    // Read the Fourier errors
    for (int i = 0; i < NFREQS; ++i) {
        if (fscanf(f, "%lf", &pfourierror[i]) != 1) {
            if (feof(f)) return 0;
            return -1;
        }
    }

    return 1;  // success
}




  
void paintKeyboardRect(FILE *PSfile,
                       double x,
                       double y_top,
                       double width,
                       double height,
                       int ntones,
                       int mintone,
                       double minnote,
                       int tonesteps)
{
    if (!PSfile || ntones <= 0) return;

    // use the SAME spacing as plotfourier()
    double ystep = height / (double)(ntones + 1);

    for (int i = 0; i < ntones; ++i) {
        // this is exactly where plotfourier() draws tone i:
        // y = y_top - (i+1)*ystep
        double y_center = y_top - (i + 1) * ystep;

        // build the semitone number that corresponds to this tone
        int note = (int)(minnote + (double)(mintone + i) / tonesteps);
        int note_in_octave = note % 12;

        int is_white =
            (note_in_octave == 0 || note_in_octave == 2 ||
             note_in_octave == 4 || note_in_octave == 5 ||
             note_in_octave == 7 || note_in_octave == 9 ||
             note_in_octave == 11);

        double gray = is_white ? 0.9 : 0.35;

        // make the band centered on the same line the plot uses
        double half = 0.5 * ystep;
        double y0 = y_center - half;
        double y1 = y_center + half;

        fprintf(PSfile, "%f setgray\n", gray);
        fprintf(PSfile,
                "newpath %.3f %.3f moveto %.3f %.3f lineto %.3f %.3f lineto %.3f %.3f lineto closepath fill\n",
                x, y0,
                x + width, y0,
                x + width, y1,
                x, y1);
    }

    // reset
    fprintf(PSfile, "0 setgray\n");
}
     





