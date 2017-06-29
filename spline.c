#include "tgfit.h"

void spline(float* x, float* y,int n, float* y2) /*** COMPUTE CUBIC SPLINE ***/
{
  int i,k;
  float p,qn,sig,un;
  y2[1]=ndiff.u[1]=0.0;
  for (i=2;i<=n-1;i++) {
    sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p=sig*y2[i-1]+2.0;
    y2[i]=(sig-1.0)/p;
    ndiff.u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    ndiff.u[i]=(6.0*ndiff.u[i]/(x[i+1]-x[i-1])-sig*ndiff.u[i-1])/p;
  }
  qn=un=0.0;
  y2[n]=(un-qn*ndiff.u[n-1])/(qn*y2[n-1]+1.0);
  for (k=n-1;k>=1;k--)
    y2[k]=y2[k]*y2[k+1]+ndiff.u[k];
}

void splint(float *xa, float* ya, float* y2a, int n, float x, float* y)
{
	int klo,khi,k;
	float h,b,a;

	klo=1;
	khi=n;
	while (khi-klo > 1) {
		k=(khi+klo) >> 1;
		if (xa[k] > x) khi=k;
		else klo=k;
	}
	h=xa[khi]-xa[klo];
	if (h == 0.0) nrerror("Bad XA input to routine SPLINT");
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;
	*y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}
