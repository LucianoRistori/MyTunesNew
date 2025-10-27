#include <stdio.h>

int flip4bytes(int input){
  int result,i;
  result=input & 0xFF;
  for(i=0;i<3;++i){
    input = input >> 8;
    result = (result << 8) | input & 0xFF;
  }
  return result;
}

int main(){

  size_t nread, nbytes, nobj;
  char bdata[64];
  int i,k,byte,idata[16];

  nread=fread(&idata,4,16,stdin);
  printf("nread=%d\n",nread);

  for(i=0;i<16;++i)idata[i]=flip4bytes(idata[i]);

  printf("RIFF=%d\n",idata[0]);
  printf("File length=%d\n",idata[1]);
  printf("WAVE=%d\n",idata[2]);
  printf("fmt=%d\n",idata[3]);
  printf("nbytes=%d\n",idata[4]);
  printf("Code=%d\n",idata[5]&0xFF);
  printf("Channels=%d\n",idata[5]>>16);
  printf("Sample rate=%d\n",idata[6]);
  printf("Bytes/sec=%d\n",idata[7]);
  printf("Block align=%d\n",idata[8]&0xFF);
  printf("Bits/sample=%d\n",idata[8]>>16);
  printf("DATA=%d\n",idata[9]);
  printf("Data length=%d\n",idata[10]);

}
