#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(void){
    int version, nfreqs, tonesteps;
    double minfreq, minnote;

    if(fread(&version,sizeof(int),1,stdin)!=1) return 1;
    if(fread(&minfreq,sizeof(double),1,stdin)!=1) return 1;
    if(fread(&minnote,sizeof(double),1,stdin)!=1) return 1;
    if(fread(&nfreqs,sizeof(int),1,stdin)!=1) return 1;
    if(fread(&tonesteps,sizeof(int),1,stdin)!=1) return 1;

    int isample; double seconds, avgvol;
    double *f = malloc(sizeof(double)*nfreqs);
    double *e = malloc(sizeof(double)*nfreqs);

    size_t recs=0;
    while (1){
        if(fread(&isample,sizeof(int),1,stdin)!=1) break;
        if(fread(&seconds,sizeof(double),1,stdin)!=1) break;
        if(fread(&avgvol,sizeof(double),1,stdin)!=1) break;
        if(fread(f,sizeof(double),nfreqs,stdin)!=(size_t)nfreqs) break;
        if(fread(e,sizeof(double),nfreqs,stdin)!=(size_t)nfreqs) break;
        recs++;
    }

    double step = exp(log(2.0)/(12.0*tonesteps));
    const int TOP=10;
    int idx[TOP]; double val[TOP];
    for(int i=0;i<TOP;i++){ idx[i]=-1; val[i]=-1e300; }

    for(int i=0;i<nfreqs;i++){
        double v=f[i];
        for(int k=0;k<TOP;k++){
            if(v>val[k]){
                for(int j=TOP-1;j>k;j--){ val[j]=val[j-1]; idx[j]=idx[j-1]; }
                val[k]=v; idx[k]=i; break;
            }
        }
    }

    printf("Records: %zu  last sample=%d  time=%.4fs\n",recs,isample,seconds);
    printf("Top %d bins:\n",TOP);
    for(int k=0;k<TOP && idx[k]>=0;k++){
        double freq=minfreq*pow(step,idx[k]);
        printf("  #%d  %.2f Hz  value=%.6g\n",k+1,freq,val[k]);
    }

    free(f); free(e);
    return 0;
}
