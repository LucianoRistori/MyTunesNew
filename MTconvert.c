//===============================================================
// File: MTconvert.c
// Purpose: Main converter / analyzer for MyTunesNew
//===============================================================
//
// Uses helper functions from:
//   - MTpkg.c   (signal and Fourier analysis)
//   - PSpkg.c   (PostScript plotting utilities)
//
// Build with: make
//===============================================================

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "MTpkg.h"
#include "PSpkg.h"

//---------------------------------------------------------------
// Main program
//---------------------------------------------------------------
int main(void)
{
    int i, retcode;
    int sample;
    double vrms = 0.0;
    double seconds = 0.0;
    double avgvol = 0.0;

    // Fourier analysis results
    double *pfourier = NULL;
    double *pfourierror = NULL;
    
    PSopen("spectrum.ps");      // Open PostScript file for output
	PSsetfont("Times-Roman", 10.0);

	const int writeEvery = 440; // about 100 frames per second at 44.1 kHz
	const int mintone = 0;      // lowest tone index
	const int ntones = 288;     // full spectrum range


    // Read, process, and analyze samples in a loop
    i = 0;
    while ((retcode = getsample(&sample, stdin)) == 0) {

        // Accumulate RMS
        vrms += (double)(sample * sample);
        ++i;

        // Perform Fourier analysis
        retcode = autofourier(sample, 0, &pfourier, &pfourierror);
        if (retcode != 0) {
            fprintf(stderr, "autofourier failed at sample %d\n", i);
            break;
        }

        // Compute volume
        avgvol = volume(sample);

        // Compute elapsed time (example placeholder)
        seconds = (double)i / 44100.0; // assuming 44.1 kHz sampling rate
		// --- Write and plot every N samples ---
		if (i % writeEvery == 0) {
			retcode = writefourier(i, seconds, avgvol, pfourier, pfourierror);
			if (retcode != 0) {
				fprintf(stderr, "writefourier failed at %.3f s\n", seconds);
				break;
			}

		// Also plot the spectrum
			plotfourier(i, seconds, avgvol, pfourier, pfourierror, mintone, ntones);
		}

		// Print RMS every so often (optional)
        if (i % 1000 == 0) {
            double vr = sqrt(vrms / i);
            fprintf(stderr, "vrms = %.0f\n", vr);
        }
    }

    // Final RMS printout
    if (i > 0) {
        double vr = sqrt(vrms / i);
        fprintf(stderr, "Final vrms = %.0f (samples: %d)\n", vr, i);
    }
    
    PSclose();  // close PostScript file


    return 0;
}
