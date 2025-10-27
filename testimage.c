//
// 
//

# include "MTpkg.c"

int main(){

  // This is a test program for postscript image

  char pic[256];
  int i;

  for(i=0;i<256;++i)pic[i]=i;
  PSimage(100.,100.,100.,400.,8,32,pic);

  PSshowpage();

  PSimage(100.,100.,200.,350.,16,16,pic);

  PSshowpage();
}
