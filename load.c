#include "tgfit.h"

void Load_Ans(char* label_ans) /*** LOAD ANSWER FILE ***/
{
  FILE *fp;
  int i,ii,j,k,l,n_nonzero_elem;
  fit.Glb_NData=0;
  fp=fopen(label_ans,"rt");
  if(fp==NULL)nrerror("Answer file not found");
  while(fgetc(fp)!='<');
  fscanf(fp, "%i", &fit.nDatasets);
  while(fgetc(fp)!='<');
  fscanf(fp, "%i", &fit.n);
  fit.Bound=vector(1,fit.n);
  fit.nParams=fit.n*(fit.n+2)+2;
  fit.SAS=vector(1,fit.n);
  model.TransferMatrix=matrix(1,fit.n,1,fit.n);
  fit.GFitVect=vector(1,fit.nParams*fit.nDatasets);
  fit.GFitList=ivector(1,fit.nParams*fit.nDatasets);
  fit.shift=fit.nParams;
  
  n_nonzero_elem=0;

  for(l=1;l<=fit.nDatasets;l++){
    Set[l].Params=vector(1,fit.nParams);
    Set[l].SAS_scale_factor=vector(1,fit.n);
    Set[l].LnkList=ivector(1,fit.nParams);
    Set[l].FitList=ivector(1,fit.nParams);

    for(i=1;i<=fit.n;i++)   /* SAS scale factors */
      Set[l].SAS_scale_factor[i]=1.0;


    while(fgetc(fp)!='<');
    Set[l].label[strlen(fgets(Set[l].label,90,fp))-1]=0;
    while(fgetc(fp)!='<');  /* Wavelength */
    fscanf(fp, "%f", &Set[l].wlg);
    while(fgetc(fp)!='<');  /* Start */
    fscanf(fp, "%i", &Set[l].Start);
    while(fgetc(fp)!='<');  /* End */
    fscanf(fp, "%i", &Set[l].End);
    fit.Glb_NData+=Set[l].End-Set[l].Start;
    while(fgetc(fp)!='<');
    for(i=1;i<=fit.n;i++)   /* Bound */
      fscanf(fp, "%f %i %i", &Set[l].Params[i],&Set[l].LnkList[i],&Set[l].FitList[i]);
    while(fgetc(fp)!='<');
    for(;i<=fit.n*2;i++){    /* SAS */
      fscanf(fp, "%f %i %i", &Set[l].Params[i],&Set[l].LnkList[i],&Set[l].FitList[i]);
        if(!ctrl.convlv)
         Set[l].Params[i]*=10000;
    }

    while(fgetc(fp)!='<');
    for(ii=1;ii<=fit.n;ii++)    /* SAS scale factors*/
      fscanf(fp, "%f", &Set[l].SAS_scale_factor[ii]);


    for(k=1;k<=fit.n;k++){  /* Transfer Matrix */
      while(fgetc(fp)!='<');
      for(j=1;j<=fit.n;j++){
	fscanf(fp, "%f %i %i", &Set[l].Params[i],&Set[l].LnkList[i],&Set[l].FitList[i]);
	if(Set[l].Params[i]!=0.0)
	  n_nonzero_elem++;
	i++;
      }
    }
    while(fgetc(fp)!='<');
    for(i=fit.n*2+1;i<=fit.nParams-1;i+=fit.n+1) /* Diagonal */
      fscanf(fp, "%f %i %i", &Set[l].Params[i],&Set[l].LnkList[i],&Set[l].FitList[i]);
    while(fgetc(fp)!='<');  /* Background */
    fscanf(fp, "%f %i %i", &Set[l].Params[fit.nParams-1],&Set[l].LnkList[fit.nParams-1],\
	   &Set[l].FitList[fit.nParams-1]);
    while(fgetc(fp)!='<');  /* Shift */
    fscanf(fp, "%f %i %i", &Set[l].Params[fit.nParams],&Set[l].LnkList[fit.nParams],\
	   &Set[l].FitList[fit.nParams]);
  }
  for(i=1;i<=fit.nDatasets;i++) {
    Set[i].nvarParams=0;
    for(j=1;j<=fit.nParams;j++)
      if(Set[i].FitList[j])Set[i].nvarParams++;
  }
  if(n_nonzero_elem <= fit.n*fit.nDatasets)
    {
      ctrl.diag=1;
      printf("** Transfer Matrix is diagonal,\n** I will use analytic rate derivatives (fast)\n");
    }
  else
    {
      ctrl.diag=0;
      printf("** Transfer Matrix is not diagonal,\n** I will use numerical rate derivatives (slow)\n");
    }
  fclose(fp);
}


void Load_Data() /*** LOAD DATA FILES ***/
{
  FILE *fp;
  int i,j;
  char linebuf[200];

  //Find number of datapoints 
  fp = fopen(Set[1].label, "rt");
  if(fp==NULL)nrerror("Error on open of data file");
  fgets(linebuf,80,fp);
  for(i=0;;i++)
    {
      if(!fgets(linebuf,80,fp)) // Read line from datafile
	break; 
    }
  ctrl.NData=i;


  printf("** Reading data: ");
  for(j=1;j<=fit.nDatasets;j++)
    {
      Set[j].data_max=0.0; Set[j].data_min=1e6;
      Set[j].res_max=-100; Set[j].res_min=1e6;
      Set[j].irf_int=0.0;
      fp = fopen(Set[j].label, "rt");
      printf("%i ",j);
      if(fp==NULL)nrerror("Error on open of data file");
      fscanf(fp, "%f", &Set[j].TCal);
      Set[j].irf=vector(1,8*ctrl.NData);
      Set[j].data=vector(1,ctrl.NData);
      fscanf(fp, "%f",  &Set[j].irf[1]);
      fscanf(fp, "%f",  &Set[j].data[1]);
      for(i=1;i<=ctrl.NData;i++){
	fscanf(fp, "%f",  &Set[j].irf[i]);
	Set[j].irf_int+=Set[j].irf[i];
	fscanf(fp, "%f",  &Set[j].data[i]);
      }
      Set[j].irf_int/=21295.2;
      fclose(fp);
    }
  printf(" >>\n");
}





void Load_Tbl(void) /*** LOAD CONTROL TABLE FILE <tgfit.tbl> ***/
{
  FILE *fp;
  int i;
  char g_str[80];
 strcpy(g_str,getenv("HOME"));
 strcat(g_str,"/.tgfit/.tgfitrc");
  fp=fopen(g_str,"rt");
  if(fp==NULL)nrerror("Table file not found");
  while(fgetc(fp)!='<');
  fscanf(fp, "%i", &ctrl.max_iterrations);
  while(fgetc(fp)!='<');
  fscanf(fp, "%f", &ctrl.tolerance);
  while(fgetc(fp)!='<');
  fscanf(fp, "%i", &ctrl.NData);
  while(fgetc(fp)!='<');
  fscanf(fp, "%i", &ctrl.NTAB);
  while(fgetc(fp)!='<');
  fscanf(fp, "%i", &ctrl.denom_param);
  while(fgetc(fp)!='<');
  fscanf(fp, "%e", &ctrl.max_lamda);
  while(fgetc(fp)!='<');
  fscanf(fp, "%i", &ctrl.SAS_positive);
  while(fgetc(fp)!='<');
  fscanf(fp, "%i", &ctrl.weight);
  while(fgetc(fp)!='<');
  fscanf(fp, "%i", &ctrl.correct);
  while(fgetc(fp)!='<');
  fscanf(fp, "%i", &ctrl.simulate);
  while(fgetc(fp)!='<');
  fscanf(fp, "%i", &ctrl.noise);
  while(fgetc(fp)!='<');
  fscanf(fp, "%f", &ctrl.noise_factor);
  while(fgetc(fp)!='<');
  fscanf(fp, "%i", &ctrl.convlv);
  while(fgetc(fp)!='<');
  fscanf(fp, "%i", &ctrl.resid);
  while(fgetc(fp)!='<');
  fscanf(fp, "%i", &ctrl.autocorr);
  while(fgetc(fp)!='<');
  fscanf(fp, "%i", &ctrl.fit);
  while(fgetc(fp)!='<');
  fscanf(fp, "%i", &ctrl.DAS);
  while(fgetc(fp)!='<');
  for(i=1;i<=fit.n;i++)
    fscanf(fp, "%i", &ctrl.save_DAS[i]);
  while(fgetc(fp)!='<');
  fscanf(fp, "%i", &ctrl.trace);
  while(fgetc(fp)!='<');
  fscanf(fp, "%i", &ctrl.trace_number);
  while(fgetc(fp)!='<');
  fscanf(fp, "%i", &ctrl.derv);
  while(fgetc(fp)!='<');
  fscanf(fp, "%i", &ctrl.derv_number);
  fclose(fp);
}
