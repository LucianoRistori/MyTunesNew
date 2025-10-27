#include <stdio.h>

int main(){

  size_t nread, nbytes, nobj;
  char data[64];
  int i,k,byte;

  nread=fread(&data,1,64,stdin);
  printf("nread=%d\n",nread);
  for(i=0;i<nread;i+=8){
    for(k=0;k<8 && (i+k)<nread;++k){
      byte=data[i+k];
      byte=byte & 0xFF;
      printf(" %02X",byte);
    }
    printf("  |  ");
    for(k=0;k<8 && (i+k)<nread;++k){
      printf("%2c",data[i+k]);
    }
    printf("\n");
  }
}
