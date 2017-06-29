#include "tgfit.h"
#include <sys/timeb.h>

void mrqcof(float* a, float** alpha, float* beta)/*** COMPILE FITING MATRIX AND EVALUATE CHISQUARE ***/
  /* Given the global parameter vector *a this function compiles (symmetric)
     alpha[nGFitVarParams][nGFitVarParams] matrix, beta[nGFitVarParams] vector
     and evaluates global chisquare */
{
  int k,j,i,l,m,ds;
  float wt,dy;
  struct  timeb start, end;
  double elapsed_t;

  ftime(&start); 
  Expand_Lnk2Params(a);

  for(ds=1;ds<=fit.nDatasets;ds++)
    {
      Expand_Params2Fit(Set[ds].Params, Set[ds].SAS_scale_factor,Set[ds].data_max);
      Exc_Flow(model.Fluor,1);
      
      if(ctrl.convlv)
	convolve(ds,Set[ds].ymod);
      else for(i=1;i<=ctrl.NData;i++)
	Set[ds].ymod[i]=model.Fluor[i];
      
      for(m=fit.n+1;m<=2*fit.n;m++)  /* Derivatives for SAS */
	{
	  if(Set[ds].FitList[m])
	    for(i=1;i<=ctrl.NData;i++)
	      model.Fluor[i]=model.trace[m-fit.n][i];    
	  
	  if(ctrl.convlv)
	    convolve(ds,Set[ds].Derivs[Set[ds].LnkList[m]]);
	  else 
	    for(i=1;i<=ctrl.NData;i++)
	      Set[ds].Derivs[Set[ds].LnkList[m]][i]=model.Fluor[i];
	}
      
  
      if(Set[ds].FitList[fit.nParams-1])  /*  Derivatives for Background = const */
	for(i=1;i<=ctrl.NData;i++)
	  Set[ds].Derivs[Set[ds].LnkList[fit.nParams-1]][i]=1;
    }

    if(ctrl.diag)  // Derivatives for rates and shift 
      Derivs_diag();
    else
    Derivs();
  for (j=1;j<=fit.nGFitParams;j++)
    { 
      for (k=1;k<=j;k++) alpha[j][k]=0.0;
      beta[j]=0.0;
    }
  fit.chisq=0.0;
  for(ds=1;ds<=fit.nDatasets;ds++)
    {
      for(i=Set[ds].Start;i<=Set[ds].End;i++)
	{
	  if((Set[ds].data[i]>0) && ctrl.weight)
	    fit.sig[i]=1/Set[ds].data[i];
	  else fit.sig[i]=1;
	  dy=Set[ds].data[i]-Set[ds].ymod[i];
	  for (j=0,l=1;l<=fit.nGFitParams;l++) 
	    {
	      if(fit.GFitList[l])
		{
		  wt=Set[ds].Derivs[l][i]*fit.sig[i];
		  for (j++,k=1,m=1;m<=l;m++)
		    if(fit.GFitList[m])
		      alpha[j][k++] += wt*Set[ds].Derivs[m][i];
		  beta[j] += dy*wt;
		}
	    }
	  fit.chisq+=dy*dy*fit.sig[i];
	}
    }
  for (j=2;j<=fit.nGFitParams;j++)
    for (k=1;k<j;k++) alpha[k][j]=alpha[j][k];


  ftime(&end); 
  elapsed_t=(1000.0*(end.time-start.time)+(end.millitm-start.millitm))/1000;
  if(ctrl.verb==2)
    fprintf(stdout, "\n Iterration time: %.3f sec\n", elapsed_t);
}



void mrqmin(void) /*** LEVENBERG-MARQUARDT MINIMIZATION ***/
{

  int k,j,l;
  fit.delta_vect=0.0;
  /* First Call */
  if(fit.lamda < 0.0){
    fit.lamda=0.001;
    mrqcof(fit.GFitVect,fit.alpha,fit.beta);
    fit.old_chisq=fit.chisq;
    for(j=1;j<=fit.nGFitParams;j++)fit.atry[j]=fit.GFitVect[j];
  }
  /* Second Call */

  for(j=1;j<=fit.nGFitVarParams;j++){
    for(k=1;k<=fit.nGFitVarParams;k++) fit.covar[j][k]=fit.alpha[j][k];
    fit.covar[j][j]=fit.alpha[j][j]*(1.0+fit.lamda);
    fit.oneda[j]=fit.beta[j];
  }
  ludcmp(fit.covar,fit.nGFitVarParams,model.gindx,&model.d);
  lubksb(fit.covar,fit.nGFitVarParams,model.gindx,fit.oneda);
  for(j=1;j<=fit.nGFitVarParams;j++)
    fit.delta_params[j]=fit.oneda[j];
  /* Last Call */
  if(fit.lamda == 0.0) {
    /* covsrt(covar,ma,lista,mfit); */
    return;
  }
  
  for(j=0,l=1;l<=fit.nGFitParams;l++)
    if(fit.GFitList[l]){
      j++;
      //     if(fabs(fit.delta_params[j]/fit.GFitVect[l]) > 0.2)
      //   fit.delta_params[j]=fit.delta_params[j]*0.2*fabs(fit.delta_params[j]/fit.GFitVect[l]);
      fit.atry[l]=fit.GFitVect[l]+fit.delta_params[j];
      fit.delta_vect+=sqrt(fit.delta_params[j]*fit.delta_params[j]/(fit.GFitVect[l]*fit.GFitVect[l]));
    }
 
  Expand_Lnk2Params(fit.atry);
  Check_Constraints();
  Parse_Links(fit.atry,0);
  Chisquare(fit.atry);
  
  if(fit.chisq < fit.old_chisq) 
    {
      mrqcof(fit.atry,fit.covar,fit.delta_params);
      fit.lamda *= 0.1;
      fit.delta_chisq=(fit.old_chisq-fit.chisq)/(fit.Glb_NData-fit.nGFitVarParams);
      fit.old_chisq=fit.chisq;
      for(j=1;j<=fit.nGFitVarParams;j++) {
	for(k=1;k<=fit.nGFitVarParams;k++)fit.alpha[j][k]=fit.covar[j][k];
	fit.beta[j]=fit.delta_params[j];
      }
      for(l=1;l<=fit.nGFitParams;l++)
	if(fit.GFitList[l])fit.GFitVect[l]=fit.atry[l];
    }
  else 
    {
      fit.lamda *= 10;
      fit.chisq=fit.old_chisq;
    }
  //  printf("%i parameters, delta =%.15e\n\n",fit.nGFitVarParams,fit.delta_vect);
}

