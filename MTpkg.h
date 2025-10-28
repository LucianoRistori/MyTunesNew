//===============================================================
// File: MTpkg.h
// Purpose: Declarations for MTpkg.c functions used by MTconvert
//===============================================================

#ifndef MTPKG_H
#define MTPKG_H

#include <stdio.h>

// --- byte-swapping helpers ---
int flip4bytes(int input);
short flip2bytes(short input);

// --- main data-handling routines ---
int getsample(int *sample, FILE *infile);
int autofourier(int sample, int opt, double **pfourier, double **pfourierror);
double volume(int sample);
int writefourier(int i, double seconds, double avgvol,
                 double *pfourier, double *pfourierror);
                 
void plotfourier(int isample, double seconds, double avgvol,
                 double *fourier, double *fourierror,
                 int mintone, int ntones);

// --- optional ---
int checkheader(FILE *fp);

#endif // MTPKG_H
