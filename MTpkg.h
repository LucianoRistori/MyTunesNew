//===============================================================
// File: MTpkg.h
// Purpose: Declarations for MTpkg.c functions used by MTconvert
//===============================================================

#ifndef MTPKG_H
#define MTPKG_H

#include <stdio.h>

//===============================================================
// Compile-time constants for Fourier analysis
//===============================================================
// lowest frequency of spectrum
#define MINFREQ (27.5)
// corresponding "note" (0.=C, 1.=C#, 2.=D, 3.=D# ...)
#define MINNOTE (9.) // it's an A
// number of frequencies in spectrum
#define NFREQS (288)
// number of frequency steps in one semitone
#define TONESTEPS (3)

//parameter to tune for perfect alignment of markers to piano keys
#define SEMITONE_SHIFT 0.40   /* positive => move markers DOWN */

// --- byte-swapping helpers ---
int flip4bytes(int input);
short flip2bytes(short input);

// --- main data-handling routines ---
int getsample(int *sample, FILE *infile);
int autofourier(int sample, int opt, double **pfourier, double **pfourierror);

int readfourier(FILE *f, int *sample, double *seconds, double *avgvol,
                double *pfourier, double *pfourierror);


double volume(int sample);
int writefourier(int i, double seconds, double avgvol,
                 double *pfourier, double *pfourierror);
                 
void writefourier_to_file(FILE *f, int i, double seconds, double avgvol,
                          double *pfourier, double *pfourierror);
                 
void plotfourier(int isample, double seconds, double avgvol,
                 double *fourier, double *fourierror,
                 int mintone, int ntones);

int checkheader(FILE *f);


// Paint piano-key background within a plot rectangle
void paintKeyboardRect(FILE *PSfile,
                       double x,
                       double y_top,
                       double width,
                       double height,
                       int ntones,
                       int mintone,
                       double minnote,
                       int tonesteps);


#endif // MTPKG_H
