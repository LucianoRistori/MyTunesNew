
//===============================================================
// File: PSpkg.c
// Purpose: PostScript utilities for MTconvert / MTpkg
//===============================================================

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "PSpkg.h"

//---------------------------------------------------------------
// Global PostScript file pointer
//---------------------------------------------------------------
FILE *PSfile = NULL;   // lazily initialized at runtime

//---------------------------------------------------------------
// Initialize PostScript output
//---------------------------------------------------------------
void PSopen(const char *filename)
{
    PSfile = fopen(filename, "w");
    if (!PSfile) {
        fprintf(stderr, "Error: cannot open %s for PostScript output\n", filename);
        exit(1);
    }

    fprintf(PSfile, "%%!PS-Adobe-3.0\n");
    fprintf(PSfile, "%%%%BoundingBox: 0 0 600 800\n");
    fprintf(PSfile, "/Times-Roman findfont 10 scalefont setfont\n");
}

//---------------------------------------------------------------
// Close PostScript output
//---------------------------------------------------------------
void PSclose(void)
{
    if (!PSfile) PSfile = stdout;
    fprintf(PSfile, "showpage\n");
    if (PSfile != stdout)
        fclose(PSfile);
    PSfile = NULL;
}

//---------------------------------------------------------------
// Write simple text at position (x, y)
//---------------------------------------------------------------
void PStext(double x, double y, const char *text)
{
    if (!PSfile) PSfile = stdout;
    fprintf(PSfile, "%.2f %.2f moveto (%s) show\n", x, y, text);
}

//---------------------------------------------------------------
// Set gray level (0=black, 1=white)
//---------------------------------------------------------------
void PSgray(double g)
{
    if (!PSfile) PSfile = stdout;
    fprintf(PSfile, "%.3f setgray\n", g);
}

//---------------------------------------------------------------
// Set line width
//---------------------------------------------------------------
void PSlinewidth(double w)
{
    if (!PSfile) PSfile = stdout;
    fprintf(PSfile, "%.3f setlinewidth\n", w);
}

//---------------------------------------------------------------
// Draw a colored line with given width
//---------------------------------------------------------------
void PSline(double x1, double y1, double x2, double y2,
            double width, double r, double g, double b)
{
    if (!PSfile) PSfile = stdout;
    fprintf(PSfile, "%.3f setlinewidth\n", width);
    fprintf(PSfile, "%.3f %.3f %.3f setrgbcolor\n", r, g, b);
    fprintf(PSfile, "%.2f %.2f moveto %.2f %.2f lineto stroke\n",
            x1, y1, x2, y2);
    fprintf(PSfile, "0 0 0 setrgbcolor\n");
}

//---------------------------------------------------------------
// Draw a filled rectangle with RGB color
//---------------------------------------------------------------
void PSrectangleFill(double x1, double y1, double x2, double y2,
                     double r, double g, double b)
{
    if (!PSfile) PSfile = stdout;
    fprintf(PSfile, "newpath %.2f %.2f moveto ", x1, y1);
    fprintf(PSfile, "%.2f %.2f lineto ", x2, y1);
    fprintf(PSfile, "%.2f %.2f lineto ", x2, y2);
    fprintf(PSfile, "%.2f %.2f lineto closepath\n", x1, y2);
    fprintf(PSfile, "%.3f %.3f %.3f setrgbcolor fill\n", r, g, b);
    fprintf(PSfile, "0 0 0 setrgbcolor\n");
}

//---------------------------------------------------------------
// Set font and size
//---------------------------------------------------------------
void PSsetfont(const char *fontname, double size)
{
    if (!PSfile) PSfile = stdout;
    fprintf(PSfile, "/%s findfont %.2f scalefont setfont\n", fontname, size);
}

//---------------------------------------------------------------
// Move current point
//---------------------------------------------------------------
void PSmoveto(double x, double y)
{
    if (!PSfile) PSfile = stdout;
    fprintf(PSfile, "%.2f %.2f moveto\n", x, y);
}

//---------------------------------------------------------------
// Show left-justified text
//---------------------------------------------------------------
void PSshowleft(const char *text)
{
    if (!PSfile) PSfile = stdout;
    fprintf(PSfile, "(%s) show\n", text);
}

//---------------------------------------------------------------
// Show right-justified text
//---------------------------------------------------------------
void PSshowright(const char *text)
{
    if (!PSfile) PSfile = stdout;
    fprintf(PSfile, "(%s) dup stringwidth pop neg 0 rmoveto show\n", text);
}

//---------------------------------------------------------------
// Draw circle outline
//---------------------------------------------------------------
void PScircle(double x, double y, double r)
{
    if (!PSfile) PSfile = stdout;
    fprintf(PSfile, "newpath %.2f %.2f %.2f 0 360 arc stroke\n", x, y, r);
}

//---------------------------------------------------------------
// Draw filled circle
//---------------------------------------------------------------
void PScircleFill(double x, double y, double r,
                  double rcol, double gcol, double bcol)
{
    if (!PSfile) PSfile = stdout;
    fprintf(PSfile, "%.3f %.3f %.3f setrgbcolor\n", rcol, gcol, bcol);
    fprintf(PSfile, "newpath %.2f %.2f %.2f 0 360 arc fill\n", x, y, r);
    fprintf(PSfile, "0 0 0 setrgbcolor\n");
}

//---------------------------------------------------------------
// Issue PostScript "showpage" command
//---------------------------------------------------------------
void PSshowpage(void)
{
    if (!PSfile) PSfile = stdout;
    fprintf(PSfile, "showpage\n");
}
