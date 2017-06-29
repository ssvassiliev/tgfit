#include "tgfit.h"


#define CON 1.4
#define CON2 (CON+CON)
#define BIG 1.0e30
#define SAFE 2.0


void Derivs() /*** COMPUTE DERIVATIVES ***/
  /* Computes derivatives of model function with respect to all fitting
     parameters at ctrl.NData points and stores them  in matrix
     Set[nDatasets].Derivs[nParams][nDatapoints].
     The value h is as an estimated initial stepsize; it need not be small,
     but rather should be an increment in x over which func changes sustantially.
     An estimate of the error in the derivative is returned as err */
{
  float errt, err, fac, h, hh, ans, TParam;
  int i,j,k,pr,l,ii;
  ans=0.0;


  for(l=1;l<=fit.nDatasets;l++)
    {
      Expand_Params2Fit(Set[l].Params,Set[l].SAS_scale_factor,Set[l].data_max);
      // Skip Bond and SAS; SAS derivatives are calculated in mrqcof() 
      for(pr=2*fit.n+1;pr<=fit.nParams;pr++)
	if(Set[l].FitList[pr]&&pr!=fit.nParams-1) // if not background perform function evaluations
	  {
	    TParam=Set[l].Params[pr];
	    h=TParam/ctrl.denom_param;
	    if(h==0)nrerror("h must be nonzero in dfridr");
	    hh=h;
	    
	    Set[l].Params[pr]=TParam+hh;
	    Expand_Params2Fit(Set[l].Params,Set[l].SAS_scale_factor,Set[l].data_max);
	    if(pr==fit.nParams) // Shift
	      {
		Exc_Flow(model.Fluor,1);
		if(ctrl.convlv)
		  convolve(l,ndiff.ymod_p[1]);
		else for(ii=1;ii<=ctrl.NData;ii++)
		  ndiff.ymod_p[1][ii]=model.Fluor[ii];
	      }
	    else Exc_Flow(ndiff.ymod_p[1],1);  // rates
	    
	    Set[l].Params[pr]=TParam-hh;
	    Expand_Params2Fit(Set[l].Params,Set[l].SAS_scale_factor,Set[l].data_max);
	    if(pr==fit.nParams) // Shift
	      {
		Exc_Flow(model.Fluor,1);
		if(ctrl.convlv)
		  convolve(l,ndiff.ymod_n[1]);
		else for(ii=1;ii<=ctrl.NData;ii++)
		  ndiff.ymod_n[1][ii]=model.Fluor[ii];
	      }
	    else Exc_Flow(ndiff.ymod_n[1],1);  // rates
	    
	    for(i=2;i<=ctrl.NTAB;i++)
	      {
		hh/=CON;
		Set[l].Params[pr]=TParam+hh;
		Expand_Params2Fit(Set[l].Params,Set[l].SAS_scale_factor,Set[l].data_max);
		if(pr==fit.nParams){ // Shift
		  Exc_Flow(model.Fluor,1);
		  if(ctrl.convlv)
		    convolve(l,ndiff.ymod_p[i]);
		  else for(ii=1;ii<=ctrl.NData;ii++)
		    ndiff.ymod_p[i][ii]=model.Fluor[ii];
		}
		else Exc_Flow(ndiff.ymod_p[i],1);  //  rates
		
		Set[l].Params[pr]=TParam-hh;
		Expand_Params2Fit(Set[l].Params,Set[l].SAS_scale_factor,Set[l].data_max);
		if(pr==fit.nParams) // Shift
		  {
		    Exc_Flow(model.Fluor,1);
		    if(ctrl.convlv)
		      convolve(l,ndiff.ymod_n[i]);
		    else for(ii=1;ii<=ctrl.NData;ii++)
		      ndiff.ymod_n[i][ii]=model.Fluor[ii];
		  }
		else Exc_Flow(ndiff.ymod_n[i],1); // rates
	      }
	    Set[l].Params[pr]=TParam;
	    Expand_Params2Fit(Set[l].Params,Set[l].SAS_scale_factor,Set[l].data_max);
	    // Calculate derivatives 
	    for(k=1;k<=ctrl.NData;k++)
	      {
		hh=h;
		ndiff.a[1][1]=(ndiff.ymod_p[1][k]-ndiff.ymod_n[1][k])/(2.0*hh);
		err=BIG;
		for (i=2;i<=ctrl.NTAB;i++)
		  {
		    hh /= CON;
		    ndiff.a[1][i]=(ndiff.ymod_p[i][k]-ndiff.ymod_n[i][k])/(2.0*hh);
		    fac=CON2;
		    for (j=2;j<=i;j++){
		      ndiff.a[j][i]=(ndiff.a[j-1][i]*fac - ndiff.a[j-1][i-1])/(fac-1.0);
		      fac=CON2*fac;
		      errt=FMAX(fabs(ndiff.a[j][i]-ndiff.a[j-1][i]),
				fabs(ndiff.a[j][i]-ndiff.a[j-1][i-1]));
		      if (errt <= err)
			{
			  err=errt;
			  ans=ndiff.a[j][i];
			}
		    }
		    if (fabs(ndiff.a[i][i]-ndiff.a[i-1][i-1]) >= SAFE*err) break;
		  }
		if(pr==fit.nParams)Set[l].Derivs[Set[l].LnkList[pr]][k]=ans; // Derivative of the shift is ready 
		else model.Fluor[k]=ans; // Derivatives of the rates need convolution
	      }

	    if(pr!=fit.nParams)
	      {
		if(ctrl.convlv)
		  convolve(l,Set[l].Derivs[Set[l].LnkList[pr]]);   
		else 
		  for(ii=1;ii<=ctrl.NData;ii++)
		    Set[l].Derivs[Set[l].LnkList[pr]][ii]=model.Fluor[ii];
	      }
	  }
    }
}



void Derivs_diag() /*** COMPUTE DERIVATIVES ***/
  /* Computes derivatives of the model function with respect to the Shift
     and Rates at ctrl.NData points and stores them  in matrix
     Set[nDatasets].Derivs[nParams][nDatapoints].
     The value h is as an estimated initial stepsize; it need not be small,
     but rather should be an increment in x over which func changes sustantially.
     An estimate of the error in the derivative is returned as err */
{
  float errt, err, fac, h, hh, TParam; 
    double A0, dA, ans;
  int i,j,k,pr,l,ii,jj,n;
  ans=0.0;
  pr=fit.nParams;

  for(l=1;l<=fit.nDatasets;l++)
    {  
      if(Set[l].FitList[pr]) // Numerical Derivatives of the  shift
	{ 
	  Expand_Params2Fit(Set[l].Params,Set[l].SAS_scale_factor,Set[l].data_max);
	  TParam=Set[l].Params[pr];// save current parameters
	  h=TParam/ctrl.denom_param;
	  if(h==0)nrerror("h must be nonzero in dfridr");
	  hh=h;	  
	  Set[l].Params[pr]=TParam+hh; // + delta shift 
	  Expand_Params2Fit(Set[l].Params,Set[l].SAS_scale_factor,Set[l].data_max);
	  Exc_Flow(model.Fluor,1);
	  if(ctrl.convlv)
	    convolve(l,ndiff.ymod_p[1]);
	  else for(ii=1;ii<=ctrl.NData;ii++)
	    ndiff.ymod_p[1][ii]=model.Fluor[ii];
	  
	  Set[l].Params[pr]=TParam-hh; // - delta shift
	  Expand_Params2Fit(Set[l].Params,Set[l].SAS_scale_factor,Set[l].data_max);      
	  Exc_Flow(model.Fluor,1);
	  if(ctrl.convlv)
	    convolve(l,ndiff.ymod_n[1]);
	  else for(ii=1;ii<=ctrl.NData;ii++)
	    ndiff.ymod_n[1][ii]=model.Fluor[ii];
	 
	  for(i=2;i<=ctrl.NTAB;i++) // fill Neville tableau
	    {
	      hh/=CON;
	      Set[l].Params[pr]=TParam+hh;
	      Expand_Params2Fit(Set[l].Params,Set[l].SAS_scale_factor,Set[l].data_max);
	      Exc_Flow(model.Fluor,1);
	      if(ctrl.convlv)
		convolve(l,ndiff.ymod_p[i]);
	      else for(ii=1;ii<=ctrl.NData;ii++)
		ndiff.ymod_p[i][ii]=model.Fluor[ii];
	      
	      Set[l].Params[pr]=TParam-hh;
	      Expand_Params2Fit(Set[l].Params,Set[l].SAS_scale_factor,Set[l].data_max);
	      Exc_Flow(model.Fluor,1);
	      if(ctrl.convlv)
		convolve(l,ndiff.ymod_n[i]);
	      else for(ii=1;ii<=ctrl.NData;ii++)
		ndiff.ymod_n[i][ii]=model.Fluor[ii];  
	    }
	  Set[l].Params[pr]=TParam; // restore parameters
	  Expand_Params2Fit(Set[l].Params,Set[l].SAS_scale_factor,Set[l].data_max);
	  
	  for(k=1;k<=ctrl.NData;k++)   // Calculate derivatives 
	    {
	      hh=h;
	      ndiff.a[1][1]=(ndiff.ymod_p[1][k]-ndiff.ymod_n[1][k])/(2.0*hh);
	      err=BIG;
	      for (i=2;i<=ctrl.NTAB;i++)
		{
		  hh /= CON;
		  ndiff.a[1][i]=(ndiff.ymod_p[i][k]-ndiff.ymod_n[i][k])/(2.0*hh);
		  fac=CON2;
		  for (j=2;j<=i;j++){
		    ndiff.a[j][i]=(ndiff.a[j-1][i]*fac - ndiff.a[j-1][i-1])/(fac-1.0);
		    fac=CON2*fac;
		    errt=FMAX(fabs(ndiff.a[j][i]-ndiff.a[j-1][i]),
			      fabs(ndiff.a[j][i]-ndiff.a[j-1][i-1]));
		    if (errt <= err)
		      {
			err=errt;
			ans=ndiff.a[j][i];
		      }
		  }
		  if (fabs(ndiff.a[i][i]-ndiff.a[i-1][i-1]) >= SAFE*err) break;
		}
	      Set[l].Derivs[Set[l].LnkList[pr]][k]=ans; // Derivatives of the shift are ready 
	    }
	}
      
      for(ii=1;ii<=fit.n;ii++) // Analytical Derivatives of the rates = -t*A*exp(-kt)
	{      
	  if(Set[l].FitList[fit.n*(ii+1)+ii])
	    {
	      A0=Set[l].Params[fit.n+ii]*Set[l].Params[ii];
	      dA=exp(-Set[l].Params[fit.n*(ii+1)+ii]*Set[1].TCal);
	      ans=A0;
	      model.Fluor[1]=0.0;
	      if(!ctrl.convlv)
		Set[l].Derivs[Set[l].LnkList[fit.n*(ii+1)+ii]][1]=model.Fluor[1];
	      for(jj=2;jj<=ctrl.NData;jj++)
		{
		  ans*=dA;
		  model.Fluor[jj]=-model.XTime[jj]*ans;  
		  if(!ctrl.convlv){
		    Set[l].Derivs[Set[l].LnkList[fit.n*(ii+1)+ii]][jj]=model.Fluor[jj];
		  }
		}
	      if(ctrl.convlv)
		convolve(l,Set[l].Derivs[Set[l].LnkList[fit.n*(ii+1)+ii]]); 
	      
	    }
	}
    }	 
}


#undef CON 
#undef CON2 
#undef BIG 
#undef SAFE 
