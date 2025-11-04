//===============================================================
// File: MTconvert.c
// Purpose: Main converter / analyzer for MyTunesNew
//
// Usage:
//   MTconvert t1 t2 input.wav spectrum.txt [--ps spectrum.ps]
//
// Examples:
//   MTconvert 0 30 test.wav spectrum.txt
//   MTconvert 10 20 song.wav out.txt --ps out.ps
//===============================================================

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "MTpkg.h"
#include "PSpkg.h"

int main(int argc, char *argv[])
{
    if (argc < 5) {
        fprintf(stderr,
            "Usage: %s t1 t2 input.wav spectrum.txt [--ps spectrum.ps]\n",
            argv[0]);
        return 1;
    }

    //------------------------------------------------------------
    // Parse required arguments
    //------------------------------------------------------------
    double t1 = atof(argv[1]);
    double t2 = atof(argv[2]);
    const char *wavfile = argv[3];
    const char *spectrumfile = argv[4];

    //------------------------------------------------------------
    // Optional --ps argument
    //------------------------------------------------------------
    int have_ps = 0;
    const char *psfile = NULL;

    if (argc >= 7 && strcmp(argv[5], "--ps") == 0) {
        psfile = argv[6];
        have_ps = 1;
    }

    //------------------------------------------------------------
    // Open input and output files
    //------------------------------------------------------------
    FILE *inwav = fopen(wavfile, "rb");
    if (!inwav) {
        perror("Cannot open input WAV file");
        return 1;
    }

    FILE *fspectrum = fopen(spectrumfile, "w");
    if (!fspectrum) {
        perror("Cannot open spectrum output file");
        fclose(inwav);
        return 1;
    }

    if (have_ps) {
        PSopen(psfile);
        fprintf(stderr, "PostScript output: %s\n", psfile);
    }

    //------------------------------------------------------------
    // Processing parameters
    //------------------------------------------------------------
    const double sampleRate = 44100.0;
    int startSample = (int)(t1 * sampleRate);
    int endSample   = (t2 > 0) ? (int)(t2 * sampleRate) : -1;

    fprintf(stderr, "Processing range: %.3f to %.3f s\n", t1, t2);

    //------------------------------------------------------------
    // Main processing loop
    //------------------------------------------------------------
    int i = 0, retcode, sample;
    double vrms = 0.0, seconds = 0.0, avgvol = 0.0;
    double *pfourier = NULL, *pfourierror = NULL;
    const int writeEvery = 441, mintone = 0, ntones = NFREQS;

    while ((retcode = getsample(&sample, inwav)) == 0) {
        ++i;
        if (i < startSample) continue;
        if (endSample > 0 && i > endSample) break;

        vrms += (double)(sample * sample);

        retcode = autofourier(sample, 0, &pfourier, &pfourierror);
        if (retcode != 0) {
            fprintf(stderr, "autofourier failed at sample %d\n", i);
            break;
        }

        avgvol = volume(sample);
        seconds = (double)i / sampleRate;

        if (i % writeEvery == 0) {
            writefourier_to_file(fspectrum, i, seconds, avgvol,
                                 pfourier, pfourierror);
            if (have_ps) {
                plotfourier(i, seconds, avgvol,
                            pfourier, pfourierror,
                            mintone, ntones);
            }
        }

        if (i % 1000 == 0)
            fprintf(stderr, "vrms = %.0f\n", sqrt(vrms / (i - startSample + 1)));
    }

    //------------------------------------------------------------
    // Final report and cleanup
    //------------------------------------------------------------
    fprintf(stderr, "Final vrms = %.0f (samples: %d)\n",
            sqrt(vrms / (i - startSample + 1)), i - startSample + 1);

    if (have_ps) PSclose();
    fclose(fspectrum);
    fclose(inwav);
    fprintf(stderr, "Done.\n");
    return 0;
}
