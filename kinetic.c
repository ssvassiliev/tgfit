#include "tgfit.h"

#define Exc_Flow_matrix Exc_Flow
//#define Exc_Flow_propagator  Exc_Flow

void Exc_Flow_propagator(float* buf, int mode)
{
// Propagate using Euler method
/*---- Constraint function is different between propagator and solution !!!! */
  int i, j,k, p;
  double v[fit.n];
  char TRANS = 'T';
  int LDA = fit.n;
  int INCX = 1, INCY = 1;
  double BETA = 1.0;
  double dt;
  int steps = 10;
  double M[fit.n][fit.n], mtmp;

  int npoints=mode*ctrl.NData;
  
  for(i=1; i <= fit.n; i++){
    for(j = 1; j <= fit.n; j++)
      M[i-1][j-1] = model.TransferMatrix[j][i]; 
  }

  for(i=0;i<fit.n;i++){  /* Fill diagonal */
    for(j=0;j<fit.n;j++)
      if(i != j)
  	M[i][i] += M[j][i];
    M[i][i] *= -1;
  }

  for(i=0; i < fit.n; i++){
    for(j = 0; j < fit.n; j++)
      printf("%10.5f",M[i][j]);printf("\n");
	}
printf("\n");
   
  for(i = 1; i <= npoints; i++) {  /* RESET TRACE MATRIX */
    buf[i]=0.0;
    for(j=1; j<=fit.n; j++)
      model.trace[j][i]=0.0;
  }

  for (i = 1; i <= fit.n; i++) {
    model.trace[i][1] = v[i-1] = fit.Bound[i];
    buf[1] += fit.SAS[i]*model.trace[i][1];
  }
    
  dt = Set[1].TCal / (double)steps;
  
  for (i = 2; i <= npoints; i++) {
      for (p = 0; p < steps; p++)            
	dgemv_ (&TRANS, &fit.n, &fit.n, &dt, &M, &LDA, &v, &INCX, &BETA, &v, &INCY);
      for (j = 1; j <= fit.n; j++) { 
	  model.trace[j][i] = v[j-1];
	  buf[i] += fit.SAS[j]*model.trace[j][i];
	}
    }
}

 

void Exc_Flow_matrix(float* buf, int mode) /*** COMPUTE MODEL DECAY CURVES ***/
     /* Solves system of rate equations, computes fluorescence decay and
	stores it in vector [buf]
	mode is a multiplier for the number of points to compute [mode*ctrl.NData]
        valid values are 1-4
	User supply model parameters in struct MODEL */
/*---- Constraint function is different betrween propagator and solution !!!! */
{
  extern void sgeev_();
  extern void sgetrf_();
  extern void sgetri_();
  float TmpMatrix[fit.n][fit.n],diag=1.;
  char jobvl = 'N', jobvr = 'V';
  int info, lwork = 4*fit.n;
  float vl[fit.n][fit.n], vr[fit.n][fit.n],
    wi[fit.n], wr[fit.n],
    work[lwork], ipiv[fit.n];
  float S[fit.n+1][fit.n+1],A[fit.n+1],Rho_0[fit.n+1],u,du;  

  int  i,j,k,l,m, npoints;

  npoints=mode*ctrl.NData;

  for(i=1;i<=npoints;i++){  /* RESET TRACE MATRIX */
    buf[i]=0.0;
    for(j=1;j<=fit.n;j++)
      model.trace[j][i]=0.0;
  } 
 
 for(i=1; i <= fit.n; i++){
    for(j = 1; j <= fit.n; j++)
      {
      S[i][j]=0.0;
      TmpMatrix[i-1][j-1] = model.TransferMatrix[i][j];
      if((i != j) && (model.TransferMatrix[i][j] != .0))diag=0;
      } 
  }

  for(i=0;i<fit.n;i++){  /* Fill diagonal */
    for(j=0;j<fit.n;j++)
      if(i != j)
  	TmpMatrix[i][i] += TmpMatrix[i][j];
    TmpMatrix[i][i] *= -1;
  }

 for(i=0; i < fit.n; i++)
    for(j = 0; j < fit.n; j++)
      TmpMatrix[i][j] *= -1;

 /*  for(i=0; i < fit.n; i++){ */
/*     for(j = 0; j < fit.n; j++) */
/*       printf("%10.5f", TmpMatrix[i][j]);printf("\n"); */
/* 	} */
/* printf("\n"); */


  /* IF MATRIX DIAGONAL  */
  if(diag){  
    for(i=1;i<=npoints;i++)
      for(j=1;j<=fit.n;j++){
	model.trace[j][i]+=exp(-model.TransferMatrix[j][j]*model.XTime[i]);
	if(fit.SAS[j])buf[i]+=fit.SAS[j]*model.trace[j][i];
      }
return;
  }

  /* MATRIX FACTORIZATION  */
  sgeev_(&jobvl, &jobvr, &fit.n, &TmpMatrix, &fit.n, &wr, &wi, &vl, &fit.n, &vr, &fit.n, &work, &lwork, &info);
  if(info)nrerror("Eigensystem solution failed");   /* FIND EIGENVALUES AND RIGHT EIGENVECTORS */
  for(i=0;i<fit.n;i++)  /* MAKE A COPY OF RIGHT EIGENVECTORS */
    for(j=0;j<fit.n;j++)vl[j][i]=vr[j][i]; 
  sgetrf_(&fit.n, &fit.n, &vl, &fit.n, &ipiv, &info);  /* COMPUTE LU FACTORIZATION OF THE MATRIX VR */
  if(info)nrerror("Matrix inversion failed (1)");  
  sgetri_(&fit.n, &vl, &fit.n, &ipiv, &work, &lwork, &info); /* COMPUTE INVERSE MATRIX: RESULT IN VL */  
  if(info)nrerror("Matrix inversion failed (2)");
 

/* CALCULATE EXCITATION PROBABILITY FLOW */

  for(i=0;i<fit.n;i++){ 
    model.EigenvalRe[i+1]=wr[i];
    for(j=0;j<fit.n;j++){
      model.EigenvecMatrix[i+1][j+1]=vr[i][j]; 
      model.EigenvecMatrixInverse[i+1][j+1]=vl[i][j]; 
    }
  }

  for(k=0;k<=fit.n;k++){Rho_0[k]= A[k]=0.0;}

    for(k=1;k<=fit.n;k++){
      u = 0;
      for(j=1;j<=fit.n;j++)
	u += model.EigenvecMatrixInverse[j][k]*fit.Bound[j];
      S[k][i]=0;
      for(i=1;i<=fit.n;i++)       
	S[k][i]+= model.EigenvecMatrix[k][i]*u;
    }

    for(i=1;i<=fit.n;i++){
      du = exp(-model.EigenvalRe[i]*Set[1].TCal);
      for(k=1;k<=fit.n;k++){
	u = S[i][k];
	model.trace[k][1] += u;
	for(l=2;l<=npoints;l++){
	  u *= du;
	  model.trace[k][l] += u;
	}
      }
    }  

    for(m=1;m<=fit.n;m++)
      for(l=1;l<=npoints;l++)
	if(fit.SAS[m])buf[l]+=fit.SAS[m]*model.trace[m][l];
}


