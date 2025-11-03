
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/* file for all .ps output */
FILE * PSfile = NULL;

/* Optional wrappers for backward compatibility */
void PSopen(const char *filename)
{
    PSfile = fopen(filename, "w");
    if (!PSfile) {
        fprintf(stderr, "Error: cannot open %s for PostScript output\n", filename);
        exit(1);
    }
    fprintf(PSfile, "%%!PS-Adobe-3.0\n");
    fprintf(PSfile, "%%%%BoundingBox: 0 0 800 600\n");
}

void PSclose(void)
{
    if (!PSfile) return;
    fprintf(PSfile, "showpage\n");
    fclose(PSfile);
    PSfile = NULL;
}




/* BASIC FUNCTIONS */ 


void PSnewpath(){
fprintf(PSfile,"newpath\n");
}
void PSrmoveto(float x, float y){
fprintf(PSfile,"%10.3f %10.3f rmoveto\n",x,y);
}
void PSmoveto(float x, float y){
fprintf(PSfile,"%10.3f %10.3f moveto\n",x,y);
}
void PSrlineto(float x, float y){
fprintf(PSfile,"%10.3f %10.3f rlineto\n",x,y);
}
void PSlineto(float x, float y){
fprintf(PSfile,"%10.3f %10.3f lineto\n",x,y);
}
void PSclosepath(){
fprintf(PSfile,"closepath\n");
}
void PSsetlinewidth(float width){
fprintf(PSfile,"%10.3f setlinewidth\n",width);
}
void PSsetdash(float n1, float n2, float n3, float n4){
fprintf(PSfile,"[%10.3f %10.3f %10.3f %10.3f] 0 setdash\n",n1,n2,n3,n4);
}
void PSstroke(){
fprintf(PSfile,"stroke\n");
}
void PSgsave(){
fprintf(PSfile,"gsave\n");
}
void PSgrestore(){
fprintf(PSfile,"grestore\n");
}
void PSsetgrey(float level){
fprintf(PSfile,"%10.3f setgrey\n",level);
}
void PSsetrgbcolor(float r, float g, float b){
fprintf(PSfile,"%10.3f %10.3f %10.3f setrgbcolor\n",r,g,b);
}
void PSfill(){
fprintf(PSfile,"fill\n");
}
void PSarc(float x, float y, float r, float a1, float a2){
fprintf(PSfile,"%10.3f %10.3f %10.3f %10.3f %10.3f arc\n",x,y,r,a1,a2);
}
void PSarcn(float x, float y, float r, float a1, float a2){
fprintf(PSfile,"%10.3f %10.3f %10.3f %10.3f %10.3f arcn\n",x,y,r,a1,a2);
}
void PSsetfont(char* fontname, float fontsize){
fprintf(PSfile,"/%s findfont\n",fontname);
fprintf(PSfile,"%10.3f scalefont setfont\n",fontsize);
}
void PSshow(char* text){
fprintf(PSfile,"(%s) show\n",text);
}
void PSshowleft(char* text){
fprintf(PSfile,"(%s) show\n",text);
}
void PSshowright(char* text){
fprintf(PSfile,"(%s) dup stringwidth pop -1 div 0 rmoveto show\n",text);
}
void PSshowcenter(char* text){
fprintf(PSfile,"(%s) dup stringwidth pop -2 div 0 rmoveto show\n",text);
}
void PSshowpage(){
fprintf(PSfile,"showpage\n");
}
void PStranslate(float x, float y){
fprintf(PSfile,"%10.3f %10.3f translate\n",x,y);
}
void PSrotate(float degrees){
fprintf(PSfile,"%10.3f rotate\n",degrees);
}

/* END OF BASIC FUNCTIONS */


/* HIGHER LEVEL FUNCTIONS */
/* Al these functions generate PS output through calls to
   basic functions ONLY */

/* define standard colors for six barrels */

#define COL0 {0.,0.,0.}
#define COL1 {1.,0.,0.}
#define COL2 {0.,.5,0.}
#define COL3 {0.,0.,1.}
#define COL4 {0.,.7,.7}
#define COL5 {.5,0.,.5}


/* Draws the border of a rectangle */

void PSrectangleBorder(float xLowerLeft, float yLowerLeft,
                 float xUpperRight, float yUpperRight,
                 float lineWidth,
                 float lineRed, float lineGreen, float lineBlue){
  PSgsave();

  PSnewpath();
  PSsetlinewidth(lineWidth);
  PSmoveto(xLowerLeft, yLowerLeft);
  PSlineto(xUpperRight,yLowerLeft);
  PSlineto(xUpperRight,yUpperRight);
  PSlineto(xLowerLeft,yUpperRight);
  PSclosepath();
  PSsetrgbcolor(lineRed, lineGreen, lineBlue);
  PSstroke();

  PSgrestore();
}

/* Draws a color filled rectangle (no border) */

void PSrectangleFill(float xLowerLeft, float yLowerLeft,
                 float xUpperRight, float yUpperRight,
                 float fillRed, float fillGreen, float fillBlue){
  PSgsave();

  PSnewpath();
  PSmoveto(xLowerLeft, yLowerLeft);
  PSlineto(xUpperRight,yLowerLeft);
  PSlineto(xUpperRight,yUpperRight);
  PSlineto(xLowerLeft,yUpperRight);
  PSclosepath();
  PSsetrgbcolor(fillRed, fillGreen, fillBlue);
  PSfill();

  PSgrestore();
}

/* Draws a color filled circle (no border) */

void PScircleFill(float xCenter, float yCenter, float radius,
                  float fillRed, float fillGreen, float fillBlue){
  PSgsave();

  PSnewpath();
  PSarc(xCenter, yCenter, radius,0.,360.);
  PSsetrgbcolor(fillRed, fillGreen, fillBlue);
  PSfill();

  PSgrestore();
}

/* Draws the border of a circle */

void PScircleBorder(float xCenter, float yCenter, float radius,
                    float lineWidth,
		    float lineRed, float lineGreen, float lineBlue){
  PSgsave();

  PSnewpath();
  PSsetlinewidth(lineWidth);
  PSarc(xCenter, yCenter, radius,0.,360.);
  PSsetrgbcolor(lineRed, lineGreen, lineBlue);
  PSstroke();

  PSgrestore();
}

/* Draws a line of specified width and color */

void PSline(float xLowerLeft, float yLowerLeft,
	    float xUpperRight, float yUpperRight,
	    float lineWidth,
	    float lineRed, float lineGreen, float lineBlue){
  PSgsave();

  PSnewpath();
  PSsetlinewidth(lineWidth);
  PSmoveto(xLowerLeft, yLowerLeft);
  PSlineto(xUpperRight,yUpperRight);
  PSsetrgbcolor(lineRed, lineGreen, lineBlue);
  PSstroke();

  PSgrestore();
}

// Paints a grey image from a bitmap (one byte per pixel)

void PSimage(float xlowleft, float ylowleft, float xlen, float ylen, int nx, int ny, char *picture){

  int ind,j;

  PSgsave();

  fprintf(PSfile,"/picstr 256 string def\n");
  fprintf(PSfile,"%.0f %.0f translate\n", xlowleft, ylowleft);
  fprintf(PSfile,"%.0f %.0f scale\n", xlen, ylen);

  fprintf(PSfile,"%d %d 8\n", nx, ny);
  fprintf(PSfile,"{%d 0 0 %d 0 0}\n", nx, ny);
  fprintf(PSfile,"{currentfile picstr readhexstring pop} image\n");

  for (ind=0;ind< nx*ny;){
    for(j=0;j<256;++j){
      if(ind < nx*ny)
	fprintf(PSfile," %02X",(*(picture+ind))&0xFF);
      else
	fprintf(PSfile," 00");
      if(ind%16 == 15) printf("\n");
      ++ind;
    }
  }
  fprintf(PSfile,"\n\n");

  PSgrestore();
}
