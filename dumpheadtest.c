#include <stdio.h>

int main(){

  size_t nread, nbytes, nobj;
  short int data[32];
  int i,k,byte;

  nread=fread(&data,2,32,stdin);
  printf("nread=%d\n",nread);
  for(i=0;i<nread;i+=4){
    for(k=0;k<4 && (i+k)<nread;++k){
      printf(" %02X",data[i+k]);
    }
    printf("\n");
  }
}
