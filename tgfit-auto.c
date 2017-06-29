#include "tgfit.h"
#include <time.h>

#define MAX_NSTART_VALUES 100 
#define Exc_Flow_matrix Exc_Flow
//#define Exc_Flow_propagator  Exc_Flow

int N_S, N_K1, N_K2, N_K3, N_K4, N_K5, N_K6;
float S[MAX_NSTART_VALUES], K1[MAX_NSTART_VALUES], K2[MAX_NSTART_VALUES], 
  K3[MAX_NSTART_VALUES], K4[MAX_NSTART_VALUES], K5[MAX_NSTART_VALUES], K6[MAX_NSTART_VALUES];

float** Load_Start_Params(int);


int main( int argc, char *argv[])
{
  double secs;
  time_t t_start, t_end;
  float **SAS; //Matrix of initial SAS
  int i,j,h,k;
  float best_chisq=1e10;
  int ii,jj,k1,k2,k3,k4,k5,k6,s,count,n_trials,cp;

  N_K1 = N_K2 = N_K3 = N_K4 = N_K5 = N_K6 = N_S = 1;

  count=0.0;cp=0;

  fprintf(stdout,"\n"
"\x1B[34m                  :->>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<-:\x1B[0m\n"
"\x1B[34m                  :->      TARGET GLOBAL DATA FIT      <-:\x1B[0m\n"
"\x1B[34m                  :->          VERSION 2.2.06          <-:\x1B[0m\n"
"\x1B[34m                  :->    (C) 2000-2006 S.Vasil'ev      <-:\x1B[0m\n"
"\x1B[34m                  :->         Brock University         <-:\x1B[0m\n"
"\x1B[34m                  :->>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<-:\x1B[0m\n\n");

  if (argc<2){
    printf(
	   "\x1B[31mSYNOPSIS\x1B[0m\n" 
	   "        tgfit-auto [INPUT.ans FILE] [OUTPUT.ans FILE - optional]\n\n" 
	   "\x1B[31mFILES\x1B[0m\n"
	   "       ~/.tgfit/.tgfitrc    - fit configuration file\n"
	   "       ~/.tgfit/auto-n.rc   - initial guesses \n"
	   "       ./[your file].ans    - model setup file\n"
	   "       ./spec.dat           - fluorescence spectrum\n\n");
    return(0);
  }


  fprintf(stdout,"** Initial command was: tgfit-auto %s %s\n",argv[1],argv[2]);  
  Load_Tbl();
  Load_Ans(argv[1]);
  Load_Data();
  Parse_Links(fit.GFitVect,1);
  strcpy(fit.ans_in,argv[1]);
  if(argc==3)strcpy(fit.ans_out,argv[2]);
  if(argc==2)strcpy(fit.ans_out,argv[1]);
  Alloc_Memory();
  SAS=Load_Start_Params(fit.n);
  ctrl.verb=0; 
  n_trials=N_K1*N_K2*N_K3*N_K4*N_K5*N_K6*N_S;
 

  fprintf(stdout,"** The number of steps in this run will be %i\n\n",n_trials);
  fflush(stdout);
  
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
  printf("\n** I will scale decay amplitudes using peak count values\n"
	 "** Use \x1B[31m SAS_scale_factor \x1B[0m for functional links: SAS[i] = k * SAS[j]  \x1B[0m\n\n");

  for(k1=0;k1<N_K1;k1++)
    for(k2=0;k2<N_K2;k2++)
      for(k3=0;k3<N_K3;k3++)
	for(k4=0;k4<N_K4;k4++)
	  for(k5=0;k5<N_K5;k5++)
	    for(k6=0;k6<N_K6;k6++)
	      for(s=0;s<N_S;s++)
		{
		  for(jj=1;jj<=fit.nDatasets;jj++)
		    {
		      for(ii=1+fit.n;ii<=2*fit.n;ii++)      
			Set[jj].Params[ii]=SAS[ii-fit.n][jj];
		      Set[jj].Params[2*fit.n+1]=K1[k1];
		      Set[jj].Params[3*fit.n+2]=K2[k2];
		      if(fit.n>=3)
			Set[jj].Params[4*fit.n+3]=K3[k3];
		      if(fit.n>=4)
			Set[jj].Params[5*fit.n+4]=K4[k4];
		      if(fit.n>=5)
			Set[jj].Params[6*fit.n+5]=K5[k5];
		      if(fit.n>=6)
			Set[jj].Params[7*fit.n+6]=K6[k6];
		      Set[jj].Params[fit.n*(fit.n+2)+1]=0.1;
		      Set[jj].Params[fit.n*(fit.n+2)+2]=S[s];
		    }
		  count++;
		  for(i=1;i<=fit.nDatasets;i++)
		    Expand_Params2Fit(Set[i].Params,Set[i].SAS_scale_factor,Set[i].data_max);
		  Parse_Links(fit.GFitVect,0);     
		  fit.lamda=-1;
		  for(h=1;h<=ctrl.max_iterrations;h++)
		    {
		      mrqmin();    
		      if((fit.delta_chisq<ctrl.tolerance)&&(fit.delta_chisq!=0))
			break;
		      if(fit.lamda>ctrl.max_lamda)
			break;
		      if(fit.delta_vect==0.0)
			break;
		    }
		  if(best_chisq>fit.chisq)
		    {
		      Save_Ans(fit.ans_out);					
		      best_chisq=fit.chisq;
		    }
		  cp++;
		  if(cp>=10)
		    {
		      fprintf(stdout,"step %i out of %i; current Chi^2 = %f, best Chi^2 = %f\n",
			     count,n_trials,fit.chisq/(fit.Glb_NData-fit.nGFitVarParams),
			     best_chisq/(fit.Glb_NData-fit.nGFitVarParams));
		      cp=0;
		      fflush(stdout);
		    }		
 		} 
  
  t_end=time(NULL);  
    
  secs=difftime(t_end,t_start);

  fprintf(stdout,"\nJob CPU time:  %g seconds\n", secs);
  return(1);
}       




float** Load_Start_Params(int n)
{
  int i,j;
  FILE *fp;
  char g_str[80];
  float **SAS;

  SAS=matrix(1,fit.n,1,fit.nDatasets);

  strcpy(g_str,getenv("HOME"));
  if(n==2)
    strcat(g_str,"/.tgfit/auto-2.rc");
  if(n==3)
    strcat(g_str,"/.tgfit/auto-3.rc");
  if(n==4)
    strcat(g_str,"/.tgfit/auto-4.rc");
  if(n==5)
    strcat(g_str,"/.tgfit/auto-5.rc");
  if(n==6)
    strcat(g_str,"/.tgfit/auto-6.rc");
  fp=fopen(g_str,"rt");
  if(fp==NULL)nrerror("Starting parameters not found");
  while(fgetc(fp)!='<');  /* SAS */
  for(j=1;j<=fit.nDatasets;j++)
    for(i=1;i<=fit.n;i++)
      SAS[i][j]=1/(float)n; // Default values
  printf("Initial SAS values:\n");
  for(j=1;j<=fit.nDatasets;j++) // Load SAS values
    {
      for(i=1;i<=fit.n;i++)
	{
	  fscanf(fp, "%f", &SAS[i][j]);
	  if(!ctrl.convlv)
	    SAS[i][j]*=10000;
	  printf("  %f ",SAS[i][j]);
	}
      printf("\n");
    }

  while(fgetc(fp)!='<');  /* K1 */
  fscanf(fp, "%i", &N_K1);
  while(fgetc(fp)!='<');
  for(i=0;i<N_K1;i++)
    fscanf(fp, "%f", &K1[i]);
  while(fgetc(fp)!='<');  /* K2 */
  fscanf(fp, "%i", &N_K2);
  while(fgetc(fp)!='<');
  for(i=0;i<N_K2;i++)
    fscanf(fp, "%f", &K2[i]);
  if(fit.n>=3)
    {
      while(fgetc(fp)!='<');  /* K3 */
      fscanf(fp, "%i", &N_K3);
      while(fgetc(fp)!='<');
      for(i=0;i<N_K3;i++)
	fscanf(fp, "%f", &K3[i]);
    }
  if(fit.n>=4)
    {
      while(fgetc(fp)!='<');  /* K4 */
      fscanf(fp, "%i", &N_K4);
      while(fgetc(fp)!='<');
      for(i=0;i<N_K4;i++)
	fscanf(fp, "%f", &K4[i]);
    }
  if(fit.n>=5)
    {
      while(fgetc(fp)!='<');  /* K5 */
      fscanf(fp, "%i", &N_K5);
      while(fgetc(fp)!='<');
      for(i=0;i<N_K5;i++)
	fscanf(fp, "%f", &K5[i]);
    }
  if(fit.n>=6)
    {
      while(fgetc(fp)!='<');  /* K5 */
      fscanf(fp, "%i", &N_K6);
      while(fgetc(fp)!='<');
      for(i=0;i<N_K6;i++)
	fscanf(fp, "%f", &K6[i]);
    }
  while(fgetc(fp)!='<');  /* SHIFT */
  fscanf(fp, "%i", &N_S);
  while(fgetc(fp)!='<');
  for(i=0;i<N_S;i++)
    fscanf(fp, "%f", &S[i]);
  fclose(fp);
  return SAS;
}


























