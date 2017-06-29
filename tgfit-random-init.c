#include "tgfit.h"
#include <time.h>
#include </usr/include/gsl/gsl_rng.h>
#include </usr/include/gsl/gsl_randist.h>
#include </usr/include/gsl/gsl_sort_float.h>

#define Exc_Flow_matrix Exc_Flow
//#define Exc_Flow_propagator  Exc_Flow


float** Load_Start_Params(int);
void Randomize(gsl_rng *);

int main( int argc, char *argv[])
{
  FILE *fp;
  double secs;
  time_t t_start, t_end;
  int i, j, h, k, l, seed, trials;
  float best_chisq;
  char label_ans[100];
  gsl_rng *r;

  size_t *permut;

  // Initialize RNG
  r = gsl_rng_alloc (gsl_rng_mt19937);
  seed = time (NULL) * getpid ();
  gsl_rng_set (r, seed);

  float **Saved_Fit_Parameters;
  float *Saved_Fit_Chisquares;
 
  /*---------- Default parameters ------------*/
  int ntrials = 10;    // Number of cycles
  int nrand = 800;     // Number of randomizations in 1 cycle
  int nbest = 3;       // Number the best for second minimization
  int iterr_1 = 100;   // Number of initial fitting iterrations
  int iterr_2 = 300;   // Number of refinement fitting iterrations 
 
  if (argc<2){
    printf(
	   "\x1B[31mSYNOPSIS\x1B[0m\n" 
	   "        tgfit-rand [INPUT.ans FILE] [OUTPUT.ans FILE - optional]\n\n" 
	   "\x1B[31mFILES\x1B[0m\n"
	   "       ~/.tgfit/.tgfitrc    - fit configuration file\n"
	   "       ~/.tgfit/?           - initial guesses \n"
	   "       ./[your file].ans    - model setup file\n"
	   "       ./spec.dat           - fluorescence spectrum\n\n");
    return(0);
  }

  printf("rand [answerfile] [ntrials]\n");
  fprintf(stdout,"** Initial command was: rand %s %s\n",argv[1],argv[2]);  
  /*--------------------- setup fitting ----------------*/
  Load_Tbl();
  Load_Ans(argv[1]);
  Load_Data();
  Parse_Links(fit.GFitVect,1);
  strcpy(fit.ans_in,argv[1]);
  ntrials=atoi(argv[2]);
  Alloc_Memory();
  permut=malloc(ntrials*sizeof(size_t));
  Saved_Fit_Parameters=matrix(1,ntrials,1,fit.nGFitParams);
  Saved_Fit_Chisquares=malloc(ntrials*sizeof(float));
  
  ctrl.verb=0; 
  // Compute IRF spline
  for(i=ctrl.NData+1,h=1;i<=2*ctrl.NData;i++,h++)
    {
      model.XTime[h]=(h-1)*Set[1].TCal;
      model.Fluor[i]=0.0;
      for(j=1;j<=fit.nDatasets;j++)
	{
	  Set[j].irf[i]=0.0; 
	  Set[j].irf_spl[i]=0.0;
	}
    }
  for(j=1;j<=fit.nDatasets;j++)
    {
      for(i=1;i<=ctrl.NData;i++)
	Set[j].data_max=FMAX(Set[j].data_max,Set[j].data[i]);
      spline(model.XTime,Set[j].irf,ctrl.NData,Set[j].spline);
    }

  t_start=time(NULL);
  
  printf("** Data peak counts:\n" );
  for(j=1,k=0;j<=fit.nDatasets;j++)
    {
      printf("%7i",(int)Set[j].data_max);
      if(k++==4)
	{
	  printf("\n");
	  k=0;
	}
    }
 
  for(i=1;i<=fit.nDatasets;i++)
    Expand_Params2Fit(Set[i].Params,Set[i].SAS_scale_factor,Set[i].data_max);
  Parse_Links(fit.GFitVect,0); 
    
  Chisquare(fit.GFitVect);
  printf("\n Initial chisq = %10.3f \n",			\
	 fit.chisq/(fit.Glb_NData-fit.nGFitVarParams));
  /*------------------------------------------------------*/

  /*------------------- MAIN LOOP----------------------- */
  for(trials=0;trials<ntrials;trials++)
    {
      best_chisq=1e30;
      /*--------------------------------------------------------*/
      printf("Random initialization %5i", trials);
      for(i=0;i<nrand;i++)
	{
	  Randomize(r);
	  Parse_Links(fit.atry,0);
	  Chisquare(fit.atry);
	  if(fit.chisq<best_chisq)
	    {
	      best_chisq=fit.chisq;
	      for(l=1;l<=fit.nGFitParams;l++)
		if(fit.GFitList[l])fit.GFitVect[l]=fit.atry[l];
	      //	      printf("Global   chisq=%f \n", fit.chisq/(fit.Glb_NData-fit.nGFitVarParams));  
	      Saved_Fit_Chisquares[trials]=fit.chisq/(fit.Glb_NData-fit.nGFitVarParams);
	    }
	}
      
      printf(" Chisq = %10.3f",Saved_Fit_Chisquares[trials]);  
      /*--------------------------------------------------------*/
      printf(" Fitting ");
      ctrl.verb=0;
      fit.lamda=-1;
      for(h=1;h<=iterr_1;h++)
	{
	  mrqmin();
	  Saved_Fit_Chisquares[trials]=fit.chisq/(fit.Glb_NData-fit.nGFitVarParams); 
	  if((fit.delta_chisq<ctrl.tolerance)&&(fit.delta_chisq!=0))
	    break;
	  if(fit.lamda>ctrl.max_lamda)
	    break;
	  if(fit.delta_vect==0.0)
	    break;
	}
      
      // Save fit 
      for(l=1;l<=fit.nGFitParams;l++)
      	Saved_Fit_Parameters[trials+1][l]=fit.GFitVect[l];
      sprintf(label_ans,"rrr_%d.ans",trials);
      Save_Ans(label_ans);
      printf("Chisq = %10.3f\n",Saved_Fit_Chisquares[trials]);
    }
  /*--------------------------------------------------------*/
 
 // Sort results 
  gsl_sort_float_index(permut,Saved_Fit_Chisquares,1,ntrials);

  printf("Refining best %i trials:\n", nbest);
  for(j=0;j<nbest;j++)
    printf("%3i %10.3f\n", j, Saved_Fit_Chisquares[permut[j]]);

  for(i=0;i<nbest;i++)
    {  
      Chisquare(Saved_Fit_Parameters[1+permut[i]]);
      /*--------------------------------------------------------*/
      printf("Fitting: %3i ",i);
      printf("Chisq = %10.3f: ",fit.chisq/(fit.Glb_NData-fit.nGFitVarParams));
      ctrl.verb=0;
      fit.lamda=-1;
      for(l=1;l<=fit.nGFitParams;l++)
	if(fit.GFitList[l])fit.GFitVect[l]=Saved_Fit_Parameters[1+permut[i]][l];
      for(h=1;h<=iterr_2;h++)
	{
	  mrqmin();
	  //	  Saved_Fit_Chisquares[trials]=fit.chisq/(fit.Glb_NData-fit.nGFitVarParams); 
	  if((fit.delta_chisq<ctrl.tolerance)&&(fit.delta_chisq!=0))
	    break;
	  if(fit.lamda>ctrl.max_lamda)
	    break;
	  if(fit.delta_vect==0.0)
	    break;
	}
      printf(" Final Chisq = %10.3f \n",			\
	     fit.chisq/(fit.Glb_NData-fit.nGFitVarParams));
      sprintf(label_ans,"rdres%d.ans",i);
      Save_Ans(label_ans);
    }
 /*------------------- END OF MAIN LOOP----------------------- */
   
  t_end=time(NULL);   
  secs=difftime(t_end,t_start);
  fprintf(stdout,"\nJob CPU time:  %g seconds\n", secs);
  return(0);
}       


void Randomize(gsl_rng *r)
{
  int i,j,k,ds;
  float *rnd_bond;
  float *SAS_center, *SAS_fwhm, *SAS_ampl;  
 
  /*---------- Default parameters ------------*/
  /* Initial populations are currently not fitted */
  float max_BOUND[6]={0.5,0.5,0.5,0.5,0.5,0.5};
  // SAS are asumed to be gaussians
  float max_DIAG_MAT=20.0;
  float max_OFF_DIAG_MAT=40.0;
  float min_MAT=0.1;
  float min_SAS_fwhm=10;
  float max_SAS_fwhm=40;
  float min_SAS_ampl=0.1;
  float max_SAS_ampl=0.7;
  /*----------------------------------------------*/
  int nc=fit.n;
  float wmax=Set[fit.nDatasets].wlg;
  float wmin=Set[1].wlg;
  
  SAS_center=malloc(nc*sizeof(float));
  SAS_fwhm=malloc(nc*sizeof(float));
  SAS_ampl=malloc(nc*sizeof(float));

  
  for(i=0;i<nc;i++) 
    SAS_ampl[i] = min_SAS_ampl + gsl_rng_uniform(r)*(max_SAS_ampl - min_SAS_ampl);
  for(i=0;i<nc;i++)
    {
      SAS_center[i] = wmin + gsl_rng_uniform(r)*(wmax-wmin);
      SAS_fwhm[i] = min_SAS_fwhm + gsl_rng_uniform(r)*(max_SAS_fwhm - min_SAS_fwhm);
    }
 
  // SAS_center[0]=720;
  //  SAS_fwhm[0]=50.0;
   
 // Randomize parameters in the first dataset
  ds=1;
  /* Randomize bond */
  for(i=1;i<=nc;i++) 
    if(Set[ds].FitList[i]){
      Set[ds].Params[i]=gsl_rng_uniform (r)*max_BOUND[i-1]+0.0001;
    }
  /* Randomize SAS */
  
  for(i=nc+1,j=0;i<=2*nc;i++,j++)
    if(Set[ds].FitList[i]){
      // Gaussian
      Set[ds].Params[i]=SAS_ampl[j]/(exp((Set[ds].wlg-SAS_center[j])*(Set[ds].wlg-SAS_center[j])/(SAS_fwhm[j]*SAS_fwhm[j])));
      // Uniform
      // Set[ds].Params[i]=SAS_ampl[j];
      // Random
      // Set[ds].Params[i] = min_SAS_ampl + gsl_rng_uniform(r)*(max_SAS_ampl - min_SAS_ampl);
    } 

  /* Normalization of SAS does not work in Bond=0! */
 /*  float sumSAS=0; */
 /*  for(i=nc+1,j=0;i<=2*nc;i++,j++) */
 /*    if(Set[ds].FitList[i]) */
 /*      { */
 /* 	sumSAS+=Set[ds].Params[i]*Set[ds].Params[i-nc]; */

 /*      } */
 
 
 /* for(i=nc+1,j=0;i<=2*nc;i++,j++) */
 /*    if(Set[ds].FitList[i]) */
 /*      { */
 /* 	Set[ds].Params[i]=Set[ds].Params[i]/sumSAS; */
 /*      } */


  /* Randomize off-diagonal matrix elements*/
  for(k=1;k<=nc;k++) 
    for(j=1;j<=nc;j++,i++)
      if(k !=j && Set[ds].FitList[i]){
	Set[ds].Params[i]=gsl_rng_uniform (r)*max_OFF_DIAG_MAT + min_MAT;	
      }	 
  /* Randomize diagonal matrix elements*/
  for(i=1;i<=nc;i++) 
    if(Set[ds].FitList[(i+1)*nc+i]){
      Set[ds].Params[(i+1)*nc+i]=gsl_rng_uniform (r)*max_DIAG_MAT + min_MAT;	  
    }   
  
  // Copy shared randomized values from the first dataset to the rest 
  for(ds=2;ds<=fit.nDatasets;ds++)
    {
      /* Randomize bond */
      for(i=1;i<=nc;i++) 
	if(Set[ds].FitList[i]){
	  Set[ds].Params[i]=Set[1].Params[i];
	}
      /* Randomize SAS */
      for(i=nc+1,j=0;i<=2*nc;i++,j++)
	if(Set[ds].FitList[i]){
	  // Gaussian
	  Set[ds].Params[i]=SAS_ampl[j]/(exp((Set[ds].wlg-SAS_center[j])*(Set[ds].wlg-SAS_center[j])/(SAS_fwhm[j]*SAS_fwhm[j])));
	  // Uniform
	  //  Set[ds].Params[i]=SAS_ampl[j];
	  // Random
	  // Set[ds].Params[i] = min_SAS_ampl + gsl_rng_uniform(r)*(max_SAS_ampl - min_SAS_ampl);
	}  

 /* Normalization of SAS does not work in Bond=0! */
      /* sumSAS=0; */
      /* for(i=nc+1,j=0;i<=2*nc;i++,j++) */
      /* 	if(Set[ds].FitList[i]) */
      /* 	  { */
      /* 	    sumSAS+=Set[ds].Params[i]*Set[ds].Params[i-nc]; */
      /* 	  } */

      /* //  sumSAS*=Set[ds].data_max/2e4; */

      /* for(i=nc+1,j=0;i<=2*nc;i++,j++) */
      /* 	if(Set[ds].FitList[i]) */
      /* 	  { */
      /* 	    Set[ds].Params[i]=Set[ds].Params[i]/sumSAS; */
      /* 	  } */


      /* Randomize off-diagonal matrix elements*/
      for(k=1;k<=nc;k++) 
	for(j=1;j<=nc;j++,i++)
	  if(k !=j && Set[ds].FitList[i]){
	    Set[ds].Params[i]=Set[1].Params[i];	
	  }	 
      /* Randomize diagonal matrix elements*/
      for(i=1;i<=nc;i++) 
	if(Set[ds].FitList[(i+1)*nc+i]){
	  Set[ds].Params[(i+1)*nc+i]=Set[1].Params[(i+1)*nc+i];	  
	}     
    }
  free(SAS_center);
  free(SAS_fwhm);
  free(SAS_ampl);
}



























