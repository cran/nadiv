#include "ga.h"

extern "C"{  

void ga(
	int *dam,
	int *sire,
	int *generation,
	int *n
){         

  int     k, kdam, ksire;

  for(k = 0; k < n[0]; k++){
     kdam = dam[k];
     ksire = sire[k];
     if((kdam != -999) & (ksire != -999)){
        generation[k] = max(generation[kdam], generation[ksire]) + 1;
     }
     else{
        if((kdam != -999)){
           generation[k] = generation[kdam] + 1;
        }
        if((ksire != -999)){
           generation[k] = generation[ksire] + 1;
        }
     }
  }
        
}
}
