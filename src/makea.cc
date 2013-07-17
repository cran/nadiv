#include "ddiag.h"

extern "C"{  

void makea(
        int *dam,       
        int *sire,         
        double *f,     
        double *dii,    
        int *N,
	int *iAP,
	int *pAP,
	double *xAP,
	int *nzmaxAP
){         

  int     i, j, k, cnt, icnt, sj, dj;
  double  ai;
  double  *AN = new double[2*N[0]];
  double  *li = new double[N[0]];
  cs *T, *D, *L, *tL, *A;

  for(i=0; i<N[0]; i++){
     li[i]=0.0;               // set l to zero
  }
  for(i=0; i<N[0]; i++){
     AN[i]=-1;               // set AN to zero
  }



  T = cs_spalloc(N[0], N[0], nzmaxAP[0], true, false);  

  D = cs_spalloc(N[0], N[0], N[0], true, false);  
        for (i = 0 ; i < N[0] ; i++){
           D->i[i] = i;
           D->x[i] = 1.0;
           D->p[i] = i;
         }
         D->p[N[0]] = N[0];  


  icnt = 0;

  for(i=0; i<N[0]; i++){  // iterate through each row of l 
  T->p[i] = icnt;
  li[i] = 1.0;                   // set l_ii to one
    ai=0.0;                        // set a_ii to zero

    dii[i] = 0.5-0.25*(f[dam[i]]+f[sire[i]]);
    D->x[i] = sqrt(dii[i]);

    j=i;
    cnt=0;

    while(j>=0){

      sj=sire[j];
      dj=dam[j];

      if(sj!= N[0]){
        AN[cnt] = sj;
        li[sj] += 0.5*li[j];
        cnt++;
      }

      if(dj!= N[0]){ 
        AN[cnt] = dj;
        li[dj] += 0.5*li[j];
        cnt++;
      }

      ai += li[j]*li[j]*dii[j];
      j=-1;

      for(k=0; k<cnt; k++){   // find eldest individual
       if(AN[k]>j){
         j = AN[k];
       }
      }
      for(k=0; k<cnt; k++){   // delete duplicates
        if(AN[k]==j){
          AN[k] -= N[0]; 
        }
      }
    }  // end of while

    f[i] = ai-1.0;
    
    for(k=0; k<i; k++){
      if(li[k] != 0.0){
        T->i[icnt] = k;
        T->x[icnt] = li[k];
        icnt++;
        li[k]  = 0.0;            
      }
    }
    T->i[icnt] = i;
    T->x[icnt] = li[i];
    li[i] = 0;
    icnt++;
  } 
  T->p[N[0]] = icnt;


  L = cs_multiply(D, T);
  tL = cs_transpose(L, true);
  A = cs_multiply(tL, L);



  for (i = 0 ; i < A->nzmax; i++){
    iAP[i] = A->i[i];
    xAP[i] = A->x[i];
  }
  for (i = 0 ; i <= A->n; i++){
    pAP[i] = A->p[i];
  }
  nzmaxAP[0] = A->nzmax;

//  nzmaxAP[0] = ((A->nzmax - A->n) / 2) + A->n;
//  for(i = 0; i <= nzmaxAP[0]; i++){
    


  free(AN);
  free(li);
  free(T);
  free(D);
  free(L);
  free(tL);
  free(A);

}
}
