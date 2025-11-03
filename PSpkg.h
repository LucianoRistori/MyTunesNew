#ifndef PSPKG_H
#define PSPKG_H

#include <stdio.h>

/* Global PostScript output file */
extern FILE *PSfile;

/* File management */
void PSopen(const char *filename);
void PSclose(void);


/* Basic PostScript primitives */
void PSnewpath(void);
void PSrmoveto(float x, float y);
void PSmoveto(float x, float y);
void PSrlineto(float x, float y);
void PSlineto(float x, float y);
void PSclosepath(void);
void PSsetlinewidth(float width);
void PSsetdash(float n1, float n2, float n3, float n4);
void PSstroke(void);
void PSgsave(void);
void PSgrestore(void);
void PSsetgrey(float level);
void PSsetrgbcolor(float r, float g, float b);
void PSfill(void);
void PSarc(float x, float y, float r, float a1, float a2);
void PSarcn(float x, float y, float r, float a1, float a2);
void PSsetfont(char *fontname, float fontsize);
void PSshow(char *text);
void PSshowleft(char *text);
void PSshowright(char *text);
void PSshowcenter(char *text);
void PSshowpage(void);
void PStranslate(float x, float y);
void PSrotate(float degrees);

/* Higher-level drawing functions */
void PSrectangleBorder(float xLowerLeft, float yLowerLeft,
                       float xUpperRight, float yUpperRight,
                       float lineWidth,
                       float lineRed, float lineGreen, float lineBlue);

void PSrectangleFill(float xLowerLeft, float yLowerLeft,
                     float xUpperRight, float yUpperRight,
                     float fillRed, float fillGreen, float fillBlue);

void PScircleFill(float xCenter, float yCenter, float radius,
                  float fillRed, float fillGreen, float fillBlue);

void PScircleBorder(float xCenter, float yCenter, float radius,
                    float lineWidth,
                    float lineRed, float lineGreen, float lineBlue);

void PSline(float xLowerLeft, float yLowerLeft,
            float xUpperRight, float yUpperRight,
            float lineWidth,
            float lineRed, float lineGreen, float lineBlue);

/* Image drawing (grayscale) */
void PSimage(float xlowleft, float ylowleft,
             float xlen, float ylen,
             int nx, int ny, char *picture);

#endif /* PSPKG_H */
