#include "tgfit.h"

void Save_residuals() /*** SAVE RESIDUALS ***/
{
  int i,j;
  FILE *fp;
  fp=fopen("residuals.dat","wt");
  if(fp==NULL)nrerror("residuals file could not be opened for writing");  
  for(i=1;i<=fit.nDatasets;i++){
    Expand_Params2Fit(Set[i].Params,Set[i].SAS_scale_factor,Set[i].data_max);
    Exc_Flow(model.Fluor,1);
    if(ctrl.convlv)
      convolve(i,Set[i].ymod);
    else
      for(j=1;j<=ctrl.NData;j++)
	Set[i].ymod[j]=model.Fluor[j];
  }
  
  for(j=1;j<=ctrl.NData;j++){
    /*fprintf(fp,"%f ",(j-1)*Set[1].TCal); x-scale*/
    for(i=1;i<=fit.nDatasets;i++)
      if(j>=Set[i].Start && j<=Set[i].End)
	fprintf(fp,"%f ",(Set[i].ymod[j]-Set[i].data[j])/sqrt(Set[i].data[j]+1)/*+(i-1)*4 - offset */);
      else fprintf(fp,"%i ",0 /*(i-1)*4* - offset*/);
    fprintf(fp,"\n");
  }
  fclose(fp);
}


void Save_autocorr(void) /*** SAVE AUTOCORRELATION ***/
{
  FILE *fp;
  int i,m1,j,npts,ds;
  float sum,denom;
  fp=fopen("autocorr.dat","wt");
    if(fp==NULL)nrerror("autocorr file could not be opened for writing");
  for(ds=1;ds<=fit.nDatasets;ds++) {
    sum = .0;  denom = .0; npts = Set[ds].End-Set[ds].Start;
    
    if(ctrl.weight)
      for(i=Set[ds].Start;i<=Set[ds].End;i++) {
	Set[ds].resid[i]=(Set[ds].ymod[i]-Set[ds].data[i])/sqrt(Set[ds].data[i]+1);
	denom+=Set[ds].resid[i]*Set[ds].resid[i];
	if (i>Set[ds].Start)sum+=(Set[ds].resid[i]-Set[ds].resid[i-1])*\
			      (Set[ds].resid[i]-Set[ds].resid[i-1]);
      }
    
    else
      for(i=Set[ds].Start;i<=Set[ds].End;i++) {
	Set[ds].resid[i]=(Set[ds].ymod[i]-Set[ds].data[i]);
	denom+=Set[ds].resid[i]*Set[ds].resid[i];
	if (i>Set[ds].Start)sum+=(Set[ds].resid[i]-Set[ds].resid[i-1])*
			      (Set[ds].resid[i]-Set[ds].resid[i-1]);
      }
   
    for(i=Set[ds].Start;i<=Set[ds].End;i+=2){
      if (Set[ds].resid[i]+Set[ds].resid[i+1]<fit.resid_min)
	fit.resid_min=Set[ds].resid[i]+Set[ds].resid[i+1];
      if (Set[ds].resid[i]+Set[ds].resid[i+1]>fit.resid_max)
	fit.resid_max=Set[ds].resid[i]+Set[ds].resid[i+1];
    }
    
    Set[ds].DW=sum/denom;
    denom/=npts;
    
    for(i=1;i<=npts/2+1;i++) {
      sum=0.0;
      m1=npts-i;
      for(j=Set[ds].Start;j<m1;j++) sum+=Set[ds].resid[j]*Set[ds].resid[j+i];
      sum/=m1;
      Set[ds].auto_corr[i]=sum/denom;
      if (Set[ds].auto_corr[i]<fit.auto_corr_min)
	fit.auto_corr_min=Set[ds].auto_corr[i];
      if (Set[ds].auto_corr[i]>fit.auto_corr_max)
	fit.auto_corr_max=Set[ds].auto_corr[i];
    }
  }
  
  for(i=1;i<=ctrl.NData/2;i++){
    for(ds=1;ds<=fit.nDatasets;ds++)
      if(i<=(Set[ds].End-Set[ds].Start)/2)
	fprintf(fp,"%f ", Set[ds].auto_corr[i]);
      else fprintf(fp,"%f ",0.0);
    fprintf(fp,"\n");
  }
  fclose(fp);
}


void Save_fit(void) /*** SAVE FIT ***/
{
  int i,j;
  FILE *fp;
  fp=fopen("fit.dat","wt");
  if(fp==NULL)nrerror("trace file could not be opened for writing");
  for(i=1;i<=fit.nDatasets;i++){
    Expand_Params2Fit(Set[i].Params,Set[i].SAS_scale_factor,Set[i].data_max);
    Exc_Flow(model.Fluor,1);
    if(ctrl.convlv)
      convolve(i,Set[i].ymod);
    else
      for(j=1;j<=ctrl.NData;j++)
	Set[i].ymod[j]=model.Fluor[j];
  }
  
  for(j=1;j<=ctrl.NData;j++){
    fprintf(fp,"%f ",j*Set[1].TCal);
    for(i=1;i<=fit.nDatasets;i++){
      if(j>=Set[i].Start && j<=Set[i].End)
	fprintf(fp,"%f %f ",Set[i].data[j],Set[i].ymod[j]);
      else fprintf(fp," 1.0 1.0 ");
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
}

void Save_Model(void) /*** SAVE DECONVOLVED DECAYS ***/
{
  int i,j;
  FILE *fp,*fp1;
  int mode=4;
  float *yield,*spec,*peak,pmax;

  spec=vector(1,fit.nDatasets);
  peak=vector(1,fit.nDatasets);
  yield=vector(1,fit.nDatasets);

  fp1=fopen("spec.dat","rt");
  if(fp1==NULL)nrerror("cannot open file spec.dat");
  for(i=1;i<=fit.nDatasets;i++)
    fscanf(fp1, "%f", &spec[i]);
    

  fp=fopen("model.dat","wt");
  if(fp==NULL)nrerror("trace file could not be opened for writing");
  int h;
  for(i=1;i<=4*ctrl.NData;i++)
    model.XTime[i]=(i-1)*Set[1].TCal;
 
  for(i=1;i<=fit.nDatasets;i++){
    yield[i]=0;
    peak[i]=0;
    Expand_Params2Fit(Set[i].Params,Set[i].SAS_scale_factor,Set[i].data_max);
    Exc_Flow(model.Fluor,mode);
    for(j=1;j<=mode*ctrl.NData;j++)
      {
	if(model.Fluor[j]>peak[i])
	  peak[i]=model.Fluor[j];
	yield[i]+=model.Fluor[j];
	Set[i].ymod[j]=model.Fluor[j];
      }
  }
  
  for(j=1;j<=ctrl.NData;j++){
    //   fprintf(fp,"%f ",model.XTime[j]);
    for(i=1;i<=fit.nDatasets;i++){
      fprintf(fp,"%f ",Set[i].ymod[j]/peak[i]);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
  free_vector(yield,1,fit.nDatasets);
  free_vector(spec,1,fit.nDatasets);
  free_vector(peak,1,fit.nDatasets);
}


void Save_DAS(void) /*** SAVE DAS ***/
{
  int i,j,k;
  FILE *fp, *fp1;
  float *spec, *yield, *ampl, **lifetimes, spec_max, ampl_max;
  
  yield=vector(1,fit.nDatasets);
  spec=vector(1,fit.nDatasets);
  ampl=vector(1,fit.nDatasets);
  lifetimes=matrix(1,fit.n,1,fit.nDatasets);
  ampl_max=0;

  fp=fopen("DAS.dat","wt");
  if(fp==NULL)nrerror("DAS.dat file could not be opened for writing");
  fp1=fopen("spec.dat","rt");
  if(fp1==NULL)nrerror("cannot open file spec.dat");
  for(i=1;i<=fit.nDatasets;i++)
    fscanf(fp1, "%f %f", &Set[i].wlg, &spec[i]);
  fclose(fp1);
   

  float  max_yield, max_spec;
  max_yield=max_spec=0.0f;
  for(j=1;j<=fit.nDatasets;j++)
    yield[j]=0.0f; 
  for(j=1;j<=fit.nDatasets;j++){
    Expand_Params2Fit(Set[j].Params,Set[j].SAS_scale_factor,Set[j].data_max); 
    Exc_Flow(model.Fluor,4);
    for(i=1;i<=fit.n;i++){
      lifetimes[i][j]=model.TransferMatrix[i][i];
    // Unscale DAS 
      model.DAS[i][j]=fit.SAS[i]*2.4e4/Set[j].data_max;
    }
    // Compute emission yield of each fluorescence decay 
    for(i=1;i<4*ctrl.NData;i++){
      yield[j] =  yield[j] + model.Fluor[i];
    }
    if(yield[j]>max_yield)max_yield=yield[j];
    if(spec[j]>max_spec)max_spec=spec[j];
  }
  
  // We will correct SAS using fluorescence yield computed by
  // numerical integration of the decay kinetics and measured fluorescence yield
  printf("%10s%10s%10s\n","Spec","Yield","Fact");
  for(j=1;j<=fit.nDatasets;j++){
    printf("%10.3e %10.3e %10.3e\n", spec[j], yield[j], spec[j]/(max_spec*yield[j]));
    for(i=1;i<=fit.n;i++){
      model.DAS[i][j] = model.DAS[i][j]*spec[j]/(max_spec*yield[j]);
    }
  }


  //  fprintf(fp,"Chisquare in this run: %f\n", fit.chisq/(fit.Glb_NData-fit.nGFitVarParams)); /* save global chisquare */
  fprintf(fp,"W");
  for(j=1;j<=fit.n;j++)
    fprintf(fp,"            T%i            A%i",j,j);
  fprintf(fp,"\n");  
  

  for(i=1;i<=fit.nDatasets;i++){
    fprintf(fp," %5.1f ", Set[i].wlg);  /* save wavelength */
    for(j=1;j<=fit.n;j++)
      if(ctrl.save_DAS[j]) fprintf(fp, "%14e%14e", 1/lifetimes[j][i],model.DAS[j][i]);    /* save SAS */
    fprintf(fp,"\n");
  }
  
  fclose(fp);
 
  free_vector(yield,1,fit.nDatasets);
  free_vector(spec,1,fit.nDatasets);
  free_vector(ampl,1,fit.nDatasets);
  free_matrix(lifetimes,1,fit.n,1,fit.nDatasets);

}


void Save_trace(int n)
{
  FILE *fp;
  int i,j;
  fp=fopen("trace.dat","wt");
  if(fp==NULL)nrerror("trace file could not be opened for writing");
  Expand_Params2Fit(Set[n].Params,Set[n].SAS_scale_factor,Set[n].data_max);
  Exc_Flow(model.Fluor,1);
  
  for(i=1;i<=ctrl.NData;i++){fprintf(fp, "%f ",(i-1)*Set[n].TCal);
  for(j=1;j<=fit.n;j++)fprintf(fp,"%f ",model.trace[j][i]);
  fprintf(fp,"\n ");}
  fclose(fp);
}

void Save_derivs(int n)
{
  FILE *fp;
  int i,j;
  fp=fopen("derivs.dat","wt");
  if(fp==NULL)nrerror("derivs file could not be opened for writing");
  Expand_Params2Fit(Set[n].Params,Set[n].SAS_scale_factor,Set[n].data_max);
  Exc_Flow(model.Fluor,1);
  
  for(i=1;i<=ctrl.NData;i++){
    fprintf(fp, "%f ",(i-1)*Set[n].TCal);
    for(j=1;j<=fit.nParams;j++)
      if(Set[n].FitList[j])
	fprintf(fp,"%f ",Set[n].Derivs[Set[n].LnkList[j]][i]);
    fprintf(fp,"\n");
  }
  fclose(fp);
}


void Save_Ans (char *label_ans) /*** SAVE ANSWER FILE ***/
{
  FILE *fp;
  int i,ii,j,k,l;
  float DW;
  fp = fopen(label_ans, "wt");
  if(fp==NULL)nrerror("Answer could not be opened for writing");
  fprintf(fp,"Target global fit>>\nInput answer file: %s\nOutput answer file: %s\n",
	  fit.ans_in, fit.ans_out);
  DW=0.0;

  Chisquare(fit.GFitVect);

  for(i=1;i<=fit.nDatasets;i++){
    DW+=Set[i].DW;
    fprintf(fp,"data [%i] chisq=%f DW=%f \n", i, Set[i].chisq, Set[i].DW);
  }

  fprintf(fp,"Global   chisq=%f DW=%f \n",\
	  fit.chisq/(fit.Glb_NData-fit.nGFitVarParams),(float)(DW/fit.nDatasets));

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

void Grow_Ans (char *label_ans, char *label_datafile, float wavelength) /*** SAVE GROWN ANSWER FILE ***/
{
  FILE *fp;
  int i,ii,j,k,l;
  float DW;
  fp = fopen(label_ans, "wt");
  if(fp==NULL)nrerror("Answer could not be opened for writing");
  fprintf(fp,"Target global fit>>\nInput answer file: %s\nOutput answer file: %s\n",
	  fit.ans_in, fit.ans_out);
  DW=0.0;

  Chisquare(fit.GFitVect);

  for(i=1;i<=fit.nDatasets;i++){
    DW+=Set[i].DW;
    fprintf(fp,"data [%i] chisq=%f DW=%f \n", i, Set[i].chisq, Set[i].DW);
  }

  fprintf(fp,"Global   chisq=%f DW=%f \n",\
	  fit.chisq/(fit.Glb_NData-fit.nGFitVarParams),(float)(DW/fit.nDatasets));

  fprintf(fp,"Number of datasets: <%i\n",fit.nDatasets+1);
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
  l=fit.nDatasets;
  int last_link=Set[l].LnkList[fit.nParams]+1;
  // INCREMENT DATASET
  fprintf(fp,"Datafile: <%s\n",label_datafile);
  fprintf(fp,"Wavelength: <%.2f\n",wavelength);
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
    fprintf(fp,"%12.9f %3i %1i ",Set[l].Params[i],last_link++,Set[l].FitList[i]);
  }
  fprintf(fp,"\nSAS SCALE FACTORS\n<");
  for(ii=1;ii<=fit.n;ii++)
    fprintf(fp," %.9f       ",Set[l].SAS_scale_factor[ii]);
  
  fprintf(fp,"\nTransfer Matrix\n");
  for(k=1;k<=fit.n;k++){
    fprintf(fp,"<");
    for(j=1;j<=fit.n;j++){
      fprintf(fp,"%12.9f %3i %1i ", model.TransferMatrix[k][j],	\
	      Set[l].LnkList[i],Set[l].FitList[i]);
      i++;
    }
    fprintf(fp,"\n");
  }
  fprintf(fp,"Diagonal\n<");
  for(i=fit.n*2+1;i<=fit.nParams-1;i+=fit.n+1) /* Diagonal */
    fprintf(fp, "%12.9f %3i %1i ", Set[l].Params[i],Set[l].LnkList[i],Set[l].FitList[i]);
  fprintf(fp,"\nBackground\n<");
  fprintf(fp, "% 12.8f %3i %1i\n", Set[l].Params[fit.nParams-1],	\
	  last_link++,Set[l].FitList[fit.nParams-1]);
  fprintf(fp,"Shift\n<");
  fprintf(fp, "% 12.9f %3i %1i\n", Set[l].Params[fit.nParams],		\
	  last_link++,Set[l].FitList[fit.nParams]);
  for(i=1;i<=19*fit.n;i++)fprintf(fp,"-");
  fprintf(fp,"\n\n");


  fclose(fp);
}

