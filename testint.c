#include <stdio.h>

int main(){

  size_t nread, nbytes, nobj;
  char bdata[16];
  int idata[16];
  short int sdata[16];

  printf("char %d %d\n", &bdata[0], &bdata[1]);
  printf("short int %d %d\n", &sdata[0], &sdata[1]);
  printf("int %d %d\n", &idata[0], &idata[1]);

}
