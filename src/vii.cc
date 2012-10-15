#include "vii.h"

extern "C"{  

void vii(
        int *N,       
	int *dam,
	double *sex,
	int *nN,
	int *nonfound,
	int *iTP,              
	int *pTP,	         
	double *xTP,
 	int *nzmaxTP,
	double *v,
	double *f
){         

  int     k, j, l, m;
  double  sumTcol;
  cs *T;

  T = cs_spalloc(N[0], N[0], nzmaxTP[0], true, false);  

         for (k = 0 ; k < nzmaxTP[0] ; k++){
           T->i[k] = iTP[k];
           T->x[k] = xTP[k];
         }
         for (k = 0 ; k <= N[0]; k++){
           T->p[k] = pTP[k];
         }

  for(j=0; j < nN[0]; j++){  //replace k with j
      k=nonfound[j];
      sumTcol=0.0;      

      if(dam[k] != N[0]){
         v[k] = 0.25*(1-f[dam[k]]);
      }
         else{
            v[k] = 0.25*(((sex[k] - 0.5) * 4) + 1.0 - f[dam[k]]);
         }
 
      for(m = T->p[k]; m < T->p[k+1]; m++){
          l = T->i[m];
          if(l <= k){
              sumTcol += T->x[m] * T->x[m] * v[l];
          }
     }

     f[k] = sumTcol - sex[k];
   }
            
  cs_spfree(T);
}
}
