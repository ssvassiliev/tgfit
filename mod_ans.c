#include "tgfit.h"

void Save_Mod_Ans (char *);

int main( int argc, char *argv[])
{
  int i,j;
  char g_str[80];
  strcpy(g_str,getenv("HOME"));
  strcat(g_str,"/.tgfit/input.tgf");

  printf("\nKeywords:\n"
	"NCompartments, NDatasets\n"
	"[Fit]InitialState, [Fit]SAS, [Fit]Background\n"
	"[Fit]Shift, [Fit]RateMatrix\n\n");
 
  Load_Ans(argv[1]);

  if(!read_input(g_str))
    {
      Save_Mod_Ans("out.ans");
      printf("** modified answer file saved as [out.ans]\n");
    }
  else 
    printf("** E ** answer file was not modified\n");
  
}

int read_input(char *filename)
{

  FILE *fp;
  char line_buf[200];
  int i,j,k,l,n,tmp, ok1, ok2;

  fp=fopen(filename, "rt");
  if(fp==NULL)
    {printf("** E ** input mod file not found\n");return 1;}

  ok1=ok2=0;
  while(1)
    {
      if(fgets(line_buf,82,fp)==NULL)
	break;
	
      /*-------------------------------------------------------------*/
      n=strlen("NCompartments");
      if(!strncmp(line_buf,"NCompartments",n))
	{
	  sscanf(&line_buf[n], "%i", &tmp);
	  printf("Ncompartments=%i\n",tmp);
	  ok1=1;
	  if(tmp != fit.n)
	    {printf("** E ** number of compartments does not match\n");return 1;}   
	}
      /*-------------------------------------------------------------*/
      n=strlen("NDatasets");
      if(!strncmp(line_buf,"NDatasets",n))
	{
	  sscanf(&line_buf[n], "%i", &tmp);
	  printf("NDatasets=%i\n",tmp);
	  ok2=1;
	  if(tmp > fit.nDatasets)
	    {printf("** E ** too many datasets\n");return 1;}   
	  else
	    fit.nDatasets=tmp;
	}
      /*-------------------------------------------------------------*/
      n=strlen("InitialState");
      if(!strncmp(line_buf,"InitialState",n))
	{
	  printf("** replacing initial states\n");
	  for(l=1;l<=fit.nDatasets;l++)
	    {
	      for(i=1;i<=fit.n;i++)
		fscanf(fp, "%f", &Set[l].Params[i]);
	    }
	}
      /*-------------------------------------------------------------*/
      n=strlen("FitInitialState");
      if(!strncmp(line_buf,"FitInitialState",n))
	{
	  printf("** replacing fit initial state flags\n");
	  for(l=1;l<=fit.nDatasets;l++)
	    {
	      for(i=1;i<=fit.n;i++)
		fscanf(fp, "%i", &Set[l].FitList[i]);
	    }
	}
      /*-------------------------------------------------------------*/
      n=strlen("SAS");
      if(!strncmp(line_buf,"SAS",n))
	{
	  printf("** replacing SAS\n");
	  for(l=1;l<=fit.nDatasets;l++)
	    {
	      for(i=1+fit.n;i<=2*fit.n;i++)
		{fscanf(fp, "%f", &Set[l].Params[i]); Set[l].Params[i]*=1e4;}
	    }
	}
     /*-------------------------------------------------------------*/
      n=strlen("FitSAS");
      if(!strncmp(line_buf,"FitSAS",n))
	{
	  printf("** replacing fit SAS flags\n");
	  for(l=1;l<=fit.nDatasets;l++)
	    {
	      for(i=1+fit.n;i<=2*fit.n;i++)
		{fscanf(fp, "%i", &Set[l].FitList[i]);}
	    }
	}
      /*-------------------------------------------------------------*/
      n=strlen("RateMatrix");
      if(!strncmp(line_buf,"RateMatrix",n))
	{
	  printf("** replacing RateMatrix\n");
	  i=2*fit.n+1;
	  for(j=1;j<=fit.n*fit.n;j++)
	    {fscanf(fp, "%f", &Set[1].Params[i]);i++;}

	  for(l=2;l<=fit.nDatasets;l++)
	    {
	      i=2*fit.n+1;
	      for(j=1;j<=fit.n*fit.n;j++)
		{Set[l].Params[i]=Set[1].Params[i];i++;}
	    }
	}
     /*-------------------------------------------------------------*/
      n=strlen("FitRateMatrix");
      if(!strncmp(line_buf,"FitRateMatrix",n))
	{
	  printf("** replacing FitRateMatrix flags\n");
	  i=2*fit.n+1;
	  for(j=1;j<=fit.n*fit.n;j++)
	    {fscanf(fp, "%i", &Set[1].FitList[i]);i++;}

	  for(l=2;l<=fit.nDatasets;l++)
	    {
	      i=2*fit.n+1;
	      for(j=1;j<=fit.n*fit.n;j++)
		{Set[l].FitList[i]=Set[1].FitList[i];i++;}
	    }
	}
  /*-------------------------------------------------------------*/
      n=strlen("FitBackground");
      if(!strncmp(line_buf,"FitBackground",n))
	{
	  printf("** replacing fit background flags\n");
	  for(l=1;l<=fit.nDatasets;l++)
	    fscanf(fp, "%i", &Set[l].FitList[fit.nParams-1]);
	}
 /*-------------------------------------------------------------*/
      n=strlen("Background");
      if(!strncmp(line_buf,"Background",n))
	{
	  printf("** replacing background\n");
	  for(l=1;l<=fit.nDatasets;l++)
	    fscanf(fp, "%f", &Set[l].Params[fit.nParams-1]);
	}
 /*-------------------------------------------------------------*/
      n=strlen("FitShift");
      if(!strncmp(line_buf,"FitShift",n))
	{
	  printf("** replacing fit shift flags\n");
	  for(l=1;l<=fit.nDatasets;l++)
	    fscanf(fp, "%i", &Set[l].FitList[fit.nParams]);
	}
 /*-------------------------------------------------------------*/
      n=strlen("Shift");
      if(!strncmp(line_buf,"Shift",n))
	{
	  printf("** replacing shift\n");
	  for(l=1;l<=fit.nDatasets;l++)
	    fscanf(fp, "%f", &Set[l].Params[fit.nParams]);
	}

   }
  if(ok1&&ok2){printf("Consistency check OK\n");return 0;}
  else {printf("Keywords: NCompartments and NDatasets are required\n");return 1;}
}


void Save_Mod_Ans (char *label_ans) /*** SAVE ANSWER FILE ***/
{
  FILE *fp;
  int i,ii,j,k,l;
  float DW;
  fp = fopen(label_ans, "wt");
  if(fp==NULL)nrerror("Answer could not be opened for writing");
 
  fprintf(fp,"Number of datasets: <%i\n",fit.nDatasets);
  fprintf(fp,"Number of compartments <%i\n\n",fit.n);
  
  
  for(l=1;l<=fit.nDatasets;l++)
    {
      Expand_Params2Fit(Set[l].Params,Set[l].SAS_scale_factor,Set[l].data_max);
      fprintf(fp,"Datafile: <%s\n",Set[l].label);
      fprintf(fp,"Wavelength: <%.2f\n",Set[l].wlg);
      fprintf(fp,"fit start: <%i ",Set[l].Start);
      fprintf(fp,"fit end: <%i\n",Set[l].End);
      for(i=1;i<=19*fit.n;i++)fprintf(fp,"-");fprintf(fp,"\n");
      fprintf(fp,"Bound\n<");
      for(i=1;i<=fit.n;i++)
	fprintf(fp,"%12.9f %3i %1i ",Set[l].Params[i],Set[l].LnkList[i],Set[l].FitList[i]);
      fprintf(fp,"\nSAS\n<");
      for(;i<=2*fit.n;i++){
	if(!ctrl.convlv)
	  Set[l].Params[i]/=10000;
	fprintf(fp,"%12.9f %3i %1i ",Set[l].Params[i],Set[l].LnkList[i],Set[l].FitList[i]);
      }
      fprintf(fp,"\nSAS SCALE FACTORS\n<");
      for(ii=1;ii<=fit.n;ii++)
	fprintf(fp," %.9f       ",Set[l].SAS_scale_factor[ii]);

      fprintf(fp,"\nTransfer Matrix\n");
      for(k=1;k<=fit.n;k++){
	fprintf(fp,"<");
	for(j=1;j<=fit.n;j++){
	  fprintf(fp,"%12.9f %3i %1i ", model.TransferMatrix[k][j],\
		  Set[l].LnkList[i],Set[l].FitList[i]);
	  i++;
	}
	fprintf(fp,"\n");
      }
      fprintf(fp,"Diagonal\n<");
      for(i=fit.n*2+1;i<=fit.nParams-1;i+=fit.n+1) /* Diagonal */
	fprintf(fp, "%12.9f %3i %1i ", Set[l].Params[i],Set[l].LnkList[i],Set[l].FitList[i]);
      fprintf(fp,"\nBackground\n<");
      fprintf(fp, "% 12.8f %3i %1i\n", Set[l].Params[fit.nParams-1],\
	      Set[l].LnkList[fit.nParams-1],Set[l].FitList[fit.nParams-1]);
      fprintf(fp,"Shift\n<");
      fprintf(fp, "% 12.9f %3i %1i\n", Set[l].Params[fit.nParams],\
	      Set[l].LnkList[fit.nParams],Set[l].FitList[fit.nParams]);
      for(i=1;i<=19*fit.n;i++)fprintf(fp,"-");
      fprintf(fp,"\n\n");
    }
  fclose(fp);
}

