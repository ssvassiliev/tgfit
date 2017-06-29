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

                  
  int iterr_1 = 500;    // Number of fitting iterrations
 

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


  fprintf(stdout,"** Initial command was: tgfit-rand %s %s\n",argv[1],argv[2]);  
  /*--------------------- setup fitting ----------------*/
  Load_Tbl();
  Load_Ans(argv[1]);
  Load_Data();
  Parse_Links(fit.GFitVect,1);
  strcpy(fit.ans_in,argv[1]);
  // if(argc==3)strcpy(fit.ans_out,argv[2]);
  if(argc==2)strcpy(fit.ans_out,argv[1]);
  Alloc_Memory();
  
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
  printf("\n Initial chisq = %10.3f \n",		\
	 fit.chisq/(fit.Glb_NData-fit.nGFitVarParams));
  
  /*------------------- MAIN LOOP----------------------- */
  
  fit.lamda=-1;
  for(h=1;h<=iterr_1;h++)
    {
      if(h%20==0)
	{ctrl.verb=1;printf("\n%i/%i",h,iterr_1);}
      mrqmin(); 
      ctrl.verb=0;
      if((fit.delta_chisq<ctrl.tolerance)&&(fit.delta_chisq!=0))
	break;
      if(fit.lamda>ctrl.max_lamda)
	break;
      if(fit.delta_vect==0.0)
	break;
    }
  
   printf("\n Final chisq = %10.3f \n",		\
	 fit.chisq/(fit.Glb_NData-fit.nGFitVarParams));
   // Save_Ans(label_ans);
  
   sprintf(label_ans,"block_%i.asc\n", atoi(argv[2]));
  Grow_Ans(fit.ans_in,label_ans,Set[fit.nDatasets].wlg+5);
 /*------------------- END OF MAIN LOOP----------------------- */
     
  t_end=time(NULL);   
  secs=difftime(t_end,t_start);
  fprintf(stdout,"\nJob CPU time:  %g seconds\n", secs);
  return(0);
}       















