//===============================================================
// File: PSpkg.h
// Purpose: PostScript utilities — header file
//===============================================================

#ifndef PSPKG_H
#define PSPKG_H

#include <stdio.h>

//---------------------------------------------------------------
// Global PostScript file handle (declared in PSpkg.c)
//---------------------------------------------------------------
extern FILE *PSfile;

//---------------------------------------------------------------
// Function prototypes
//---------------------------------------------------------------
void PSopen(const char *filename);
void PSclose(void);
void PStext(double x, double y, const char *text);
void PSgray(double g);
void PSlinewidth(double w);
void PSline(double x1, double y1, double x2, double y2,
            double width, double r, double g, double b);
void PSrectangleFill(double x1, double y1, double x2, double y2,
                     double r, double g, double b);
void PSsetfont(const char *fontname, double size);
void PSmoveto(double x, double y);
void PSshowleft(const char *text);
void PSshowright(const char *text);
void PScircle(double x, double y, double r);
void PScircleFill(double x, double y, double r,
                  double rcol, double gcol, double bcol);
void PSshowpage(void);

#endif // PSPKG_H
