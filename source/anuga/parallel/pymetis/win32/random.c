#include <stdlib.h>

double drand48(){
  return (((double)rand())/RAND_MAX);
}

void srand48(long seed){
  srand(seed);
}
