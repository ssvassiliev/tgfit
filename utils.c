#include "tgfit.h"


void Alloc_Memory(void)
{
int i;
  model.Fluor=vector(1,8*ctrl.NData);
  model.TmpMatrix=matrix(1,fit.n,1,fit.n);
  model.EigenvecMatrix=matrix(1,fit.n,1,fit.n);
  model.EigenvecMatrixInverse=matrix(1,fit.n,1,fit.n);
  model.gindx=ivector(1,fit.nGFitParams);

  model.EigenvalRe=vector(1,fit.n);
  model.trace=matrix(1,fit.n,1,4*ctrl.NData);
  model.DAS=matrix(1,fit.n,1,fit.nDatasets);
  model.XTime=vector(1,4*ctrl.NData);
  ndiff.ymod_p=matrix(1,ctrl.NTAB,1,ctrl.NData);
  ndiff.ymod_n=matrix(1,ctrl.NTAB,1,ctrl.NData);
  ndiff.a=matrix(1,ctrl.NTAB,1,ctrl.NTAB);
  ndiff.u=vector(1,ctrl.NData);
  fit.oneda=vector(1,fit.nGFitParams);
  fit.atry=vector(1,fit.nGFitParams);
  fit.delta_params=vector(1,fit.nGFitParams);
  fit.alpha=matrix(1,fit.nGFitParams,1,fit.nGFitParams);
  fit.covar=matrix(1,fit.nGFitParams,1,fit.nGFitParams);
  fit.Fluor=vector(1,4*ctrl.NData);
  fit.beta=vector(1,fit.nGFitParams);
  fit.sig=vector(1,ctrl.NData);
  fit.corr=vector(1,ctrl.NData);

  for(i=1;i<=fit.nDatasets;i++){
    Set[i].ymod=vector(1,8*ctrl.NData);
    Set[i].Derivs=matrix(1,fit.nGFitParams,1,ctrl.NData);
    Set[i].auto_corr=vector(1,(Set[i].End-Set[i].Start)/2);
    Set[i].resid=vector(1,ctrl.NData);
    Set[i].irf_spl=vector(1,8*ctrl.NData);
    Set[i].spline=vector(1,8*ctrl.NData);
  } 
}
 



void Parse_Links(float *fit_vect, int first_call)
     /*** PARSE GLOBAL LINKS AND COMPILE GLOBAL FIT PARAMETER VECTOR ***/
{
  int i,k,nlp=0;
  
  for(k=1;k<=fit.nDatasets;k++)
    {
      for(i=1;i<=fit.nParams;i++)
	if(Set[k].LnkList[i]>nlp){
	  fit_vect[nlp+1]=Set[k].Params[i];
	  if(first_call)fit.GFitList[nlp+1]=Set[k].FitList[i];
	  nlp++;
	}
    }
  if(first_call){
    fit.nGFitParams=nlp;
    fit.nGFitVarParams=nlp;
    for(i=1;i<=nlp;i++)
      if(!fit.GFitList[i])fit.nGFitVarParams--;
  }
  //  printf("FITTING PARAMETERS = %i\n",fit.nGFitVarParams);
}


void Expand_Lnk2Params(float *a)
{
  /*** EXPAND GLOBAL FIT PARAMETER VECTOR TO ARRAY DATASET[].PARAMETERS[] ***/
  int i,k;
  for(k=1;k<=fit.nDatasets;k++)
    for(i=1;i<=fit.nParams;i++)
      if(Set[k].FitList[i])Set[k].Params[i]=a[Set[k].LnkList[i]];
}

void Simulate(void)
{
  int i,ds;                   
  FILE *fp;
  printf("Saving  data:");
  for(ds=1;ds<=fit.nDatasets;ds++){
    Set[ds].data_max=0;
    /* fp = fopen(Set[ds].label, "wt"); */
    /* if(fp==NULL)nrerror("simulated file could not be opened for writing"); */
    /* printf(" %i",ds); */
    /* if(fp==NULL)nrerror("Error on open of data file"); */
    Expand_Params2Fit(Set[ds].Params,Set[ds].SAS_scale_factor,Set[ds].data_max);
    Exc_Flow(model.Fluor,1);
    
    if(ctrl.convlv)
      convolve(ds,Set[ds].ymod);
    else for(i=1;i<=ctrl.NData;i++)
	   Set[ds].ymod[i]=model.Fluor[i]*10000;
    
    /* fprintf(fp,"%f %i %i\n",Set[ds].TCal,Set[ds].Start,Set[ds].End); */
    for(i=1;i<=ctrl.NData;i++){
      if(Set[ds].ymod[i]<0)Set[ds].ymod[i]=0.001;
      if(ctrl.noise)  //  Add noise to simulated curve
    	Set[ds].ymod[i]+=sqrt(Set[ds].ymod[i])*(drand48()-0.5)*ctrl.noise_factor;
      if(Set[ds].ymod[i]<0)
    	Set[ds].ymod[i]=0.001;
    /*   fprintf(fp,"%f %f\n",Set[ds].irf[i],Set[ds].ymod[i]); */
      Set[ds].data_max=FMAX(Set[ds].data_max,Set[ds].ymod[i]);
    }
    /* fclose(fp); */
  }
printf("  >>\n");
}


void Check_Constraints(void)
{
  int i,j,k,ds;
  int *bcount, *scount, *mcount, *dmcount;
  
  bcount=calloc(1,fit.nDatasets*sizeof(int));
  scount=calloc(1,fit.nDatasets*sizeof(int));
  mcount=calloc(1,fit.nDatasets*sizeof(int));
  dmcount=calloc(1,fit.nDatasets*sizeof(int));

  for(ds=1;ds<=fit.nDatasets;ds++){
    for(i=1;i<=fit.n;i++)  /* check bond */
      if(Set[ds].Params[i]<0){
	Set[ds].Params[i]=-Set[ds].Params[i]/5;
	bcount[ds-1]++; 
      }
    
    if(ctrl.SAS_positive)
      {
	for(;i<=2*fit.n;i++)  /* check SAS */
	  if(Set[ds].Params[i]<0){
	    Set[ds].Params[i]=-Set[ds].Params[i]/5;
	    scount[ds-1]++;
	  }
      }
    for(;i<=2*fit.n;i++);    

    for(k=1;k<=fit.n;k++) /* check off-diagonal matrix elements*/
      for(j=1;j<=fit.n;j++,i++)
	if(k!=j&&Set[ds].Params[i]<0){
	  Set[ds].Params[i]=-Set[ds].Params[i]/5;
	  mcount[ds-1]++; 	
	}
	  
    for(i=1;i<=fit.n;i++) /* check diagonal matrix elements*/
      if(Set[ds].Params[(i+1)*fit.n+i]<0){
	Set[ds].Params[(i+1)*fit.n+i]=-Set[ds].Params[(i+1)*fit.n+i];
	dmcount[ds-1]++;	  
      }
  } 

  if(ctrl.verb==2)
    {
      printf("\nConstrained BND ");
      for(i=0;i<fit.nDatasets;i++)
	printf("%i ",bcount[i]);
      printf("\n");
      printf("Constrained SAS ");
      for(i=0;i<fit.nDatasets;i++)
	printf("%i ",scount[i]);
      printf("\n");
      printf("Constrained ME  ");
      for(i=0;i<fit.nDatasets;i++)
	printf("%i ",mcount[i]);
      printf("\n");
      printf("Constrained DME ");
      for(i=0;i<fit.nDatasets;i++)
	printf("%i ",dmcount[i]);
      printf("\n");
    }
  free(bcount); free(scount); free(mcount); free(dmcount);
}



void convolve(int ds, float* buf) /*** CONVOLUTION ***/
{
  int i;
  /* Convolve model.Fluor[] with shifted irf from Set[ds],
     result returned in array buf[] */

  for(i=1;i<=ctrl.NData;i++){
    if((i-1-fit.shift)>=0.0)
      splint(model.XTime,Set[ds].irf,Set[ds].spline,ctrl.NData,
	     (i-1-fit.shift)*Set[ds].TCal,&Set[ds].irf_spl[i]);
    else Set[ds].irf_spl[i]=Set[ds].irf[i];
    if(Set[ds].irf_spl[i]<0.0)Set[ds].irf_spl[i]=0.0;
  }
 
  convlv(model.Fluor,ctrl.NData,Set[ds].irf_spl,fit.Fluor);
#pragma omp parallel for
  for(i=1;i<=ctrl.NData;i++)buf[i]=fit.Fluor[i]/Set[ds].irf_int+fit.bkgnd;
}

void Expand_Params2Fit(float *params, float *sas_scaling, float d_max)
{
  int i,j,k;
  /* EXPAND LIST OF PARAMETERS TO STRUCT FIT */
  for(i=1;i<=fit.n;i++)fit.Bound[i]=params[i];
  // SAS of each dataset are scaled by data maximum 
  for(j=1;j<=fit.n;j++,i++)fit.SAS[j]=sas_scaling[j]*params[i]*d_max/2.4e4;

  for(k=1;k<=fit.n;k++) /* Fill transfer matrix */
    for(j=1;j<=fit.n;j++,i++)
      model.TransferMatrix[k][j]=params[i];
  
  fit.bkgnd=params[fit.nParams-1];
  fit.shift=params[fit.nParams];
}



void Chisquare(float* a) /*** COMPUTE CHISQUARE ***/
  /* Given global parameter vector *a this function evaluates global
     chisquare and displays chisquares of each dataset */
{
  int i,ds;
  float dy,tmp_chisq;
  
  Expand_Lnk2Params(a);
  if(ctrl.verb) 
    fprintf(stdout,"* "); 
  fit.chisq=0.0;
  tmp_chisq=0.0;

  for(ds=1;ds<=fit.nDatasets;ds++)
    {
      Set[ds].res_max=-1e8;
      Set[ds].res_min=1e8;
    Expand_Params2Fit(Set[ds].Params,Set[ds].SAS_scale_factor,Set[ds].data_max);
    Exc_Flow(model.Fluor,1);
    
    if(ctrl.convlv)
      convolve(ds,Set[ds].ymod);
    else for(i=1;i<=ctrl.NData;i++)
      Set[ds].ymod[i]=model.Fluor[i];
    
    for(i=Set[ds].Start;i<=Set[ds].End;i++){
      if(Set[ds].data[i]>0&&ctrl.weight)fit.sig[i]=1/Set[ds].data[i];
      else fit.sig[i]=1;
      dy=Set[ds].data[i]-Set[ds].ymod[i];
      fit.chisq+=dy*dy*fit.sig[i];
      Set[ds].resid[i]=dy*sqrt(fit.sig[i]);
      Set[ds].res_max=FMAX(Set[ds].res_max,Set[ds].resid[i]);
      Set[ds].res_min=FMIN(Set[ds].res_min,Set[ds].resid[i]);
      tmp_chisq+=dy*dy*fit.sig[i];
    }
    Set[ds].chisq=tmp_chisq/(Set[ds].End-Set[ds].Start-Set[ds].nvarParams);
    if(ctrl.verb) 
      fprintf(stdout,"%.3f ",Set[ds].chisq);
    tmp_chisq=0.0;
  } 
  if(ctrl.verb) 
    fprintf(stdout,"*  %.3f [%.3f] L=%.0e",
	    fit.chisq/(fit.Glb_NData-fit.nGFitVarParams),
	    fit.old_chisq/(fit.Glb_NData-fit.nGFitVarParams), fit.lamda);
  fflush(stdout);
}





