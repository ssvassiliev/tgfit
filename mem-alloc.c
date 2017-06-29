
#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>


void nrerror(char* error_text)
{
  fprintf(stderr," ");
  fprintf(stderr,"\n\nTGFIT run-time error: %s\n\n",error_text);
  exit(0);
}

double *dvector(int nl, int nh)
{
  double *v;
  v=NULL;
  v=(double *)malloc((unsigned long) (nh-nl+1)*sizeof(double));
  if (v==NULL) nrerror("allocation failure in vector()");
  return v-nl;
}

float *vector(int nl, int nh)
{
  float *v;
  v=NULL;
  v=(float *)malloc((unsigned long) (nh-nl+1)*sizeof(double));
  if (v==NULL) nrerror("allocation failure in vector()");
  return v-nl;
}

int *ivector(int nl, int nh)
{
  int *v;
  v=NULL;
  v=(int *)malloc((unsigned long) (nh-nl+1)*sizeof(int));
  if (v==NULL) nrerror("allocation failure in ivector()");
  return v-nl;
}

float **matrix(int nrl,int nrh,int ncl,int nch)
{
  int i;
  float **m;
  m=NULL;
  m=(float **) malloc((unsigned long) (nrh-nrl+1)*sizeof(float*));
  if (m==NULL) nrerror("allocation failure 1 in matrix()");
  m -= nrl;
  
  for(i=nrl;i<=nrh;i++) {
    m[i]=(float *) malloc((unsigned long) (nch-ncl+1)*sizeof(float));
    if (m[i]==NULL) nrerror("allocation failure 2 in matrix()");
    m[i] -= ncl;
  }
  return m;
}

void free_vector(float* v,int nl,int nh){free( (v+nl));}
void free_dvector(double* v,int nl,int nh){free( (v+nl));}
void free_ivector(int* v, int nl, int nh){ free( (v+nl));}

void free_matrix(float** m,int nrl,int nrh,int ncl,int nch)
{
  int i;
  for(i=nrh;i>=nrl;i--) free( (m[i]+ncl));
  free( (m+nrl));
}






















