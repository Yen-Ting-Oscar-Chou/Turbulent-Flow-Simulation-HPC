#include  <stdlib.h>
#include   <stdio.h>
#define N 100000

typedef struct myvec{
    size_t len;
    double *data;
} myvec_t;

void init(myvec_t *s);

int main(){
   myvec_t s;

   s.data = (double *)calloc(N,sizeof(double));
   s.len  = N;

   #pragma omp target parallel for defaultmap(default: pointer) map(s, s.data[0:s.len])
   for(int i=0; i<s.len; i++) {
    s.data[i]=3*i+1;
   } 

   printf("s.data[%d]=%lf\n",N-1,s.data[N-1]);  //s.data[99]=99.000000
}
