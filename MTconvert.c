//===============================================================
// File: MTconvert.c
// Purpose: Main converter / analyzer for MyTunesNew
//
// Usage:
//   MTconvert input.wav
//       -> reads input.wav, writes spectrum text to stdout, no PS
//
//   MTconvert input.wav spectrum.txt
//       -> reads input.wav, writes spectrum text to spectrum.txt, no PS
//
//   MTconvert input.wav spectrum.txt spectrum.ps
//       -> reads input.wav, writes spectrum text to spectrum.txt,
//          AND writes per-frame PS plots to spectrum.ps
//
//   (you can also do)
//   MTconvert input.wav - spectrum.ps
//       -> spectrum to stdout, PS to spectrum.ps
//
// Notes:
//   - actual spectrum writing is still done by writefourier(...) in MTpkg.c
//   - actual plotting is still done by plotfourier(...) in MTpkg.c
//===============================================================

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "MTpkg.h"
#include "PSpkg.h"

int main(int argc, char *argv[])
{
    FILE *inwav = stdin;
    const char *spectrum_out = NULL;
    const char *ps_out = NULL;
    int have_ps = 0;

    // ------------------------------------------------------------
    // Parse command line
    // ------------------------------------------------------------
    // argv[1] = input wav
    // argv[2] = spectrum text file (or "-" for stdout)
    // argv[3] = PS output file
    if (argc < 2) {
        fprintf(stderr,
                "Usage: %s input.wav [spectrum.txt|-] [spectrum.ps]\n",
                argv[0]);
        return 1;
    }

    // open input wav
    inwav = fopen(argv[1], "rb");
    if (!inwav) {
        perror("Cannot open input wav");
        return 1;
    }

    // spectrum destination
    if (argc >= 3) {
        spectrum_out = argv[2];
        if (strcmp(spectrum_out, "-") != 0) {
            // redirect stdout to spectrum file
            if (!freopen(spectrum_out, "w", stdout)) {
                perror("Cannot open spectrum output file");
                fclose(inwav);
                return 1;
            }
        }
        // else: leave stdout as is
    }

    // PS destination (optional)
    if (argc >= 4) {
        ps_out = argv[3];
        PSopen(ps_out);
        have_ps = 1;
    } else {
        have_ps = 0; // no PS output
    }

    // ------------------------------------------------------------
    // Processing loop (mostly your original code)
    // ------------------------------------------------------------
    int i = 0;
    int retcode;
    int sample;
    double vrms = 0.0;
    double seconds = 0.0;
    double avgvol = 0.0;
    double *pfourier = NULL;
    double *pfourierror = NULL;

    const int writeEvery = 441;   // every 441 samples (~0.01s @ 44.1kHz)
    const int mintone    = 0;
    const int ntones     = NFREQS;

    while ((retcode = getsample(&sample, inwav)) == 0) {

        // accumulate RMS
        vrms += (double)(sample * sample);
        ++i;

        // feed analyzer
        retcode = autofourier(sample, 0, &pfourier, &pfourierror);
        if (retcode != 0) {
            fprintf(stderr, "autofourier failed at sample %d\n", i);
            break;
        }

        // running volume
        avgvol = volume(sample);

        // time in seconds
        seconds = (double)i / 44100.0;

        // every N samples, dump a frame (text) and optionally make a PS page
        if (i % writeEvery == 0) {
            retcode = writefourier(i, seconds, avgvol, pfourier, pfourierror);
            if (retcode != 0) {
                fprintf(stderr, "writefourier failed at %.3f s\n", seconds);
                break;
            }

            if (have_ps) {
                // this uses global PSfile inside PSpkg.c
                plotfourier(i, seconds, avgvol,
                            pfourier, pfourierror,
                            mintone, ntones);
            }
        }

        // optional monitoring
        if (i % 1000 == 0) {
            double vr = sqrt(vrms / i);
            fprintf(stderr, "vrms = %.0f\n", vr);
        }
    }

    // final RMS
    if (i > 0) {
        double vr = sqrt(vrms / i);
        fprintf(stderr, "Final vrms = %.0f (samples: %d)\n", vr, i);
    }

    // close PS if we opened it
    if (have_ps) {
        PSclose();
    }

    fclose(inwav);
    return 0;
}
