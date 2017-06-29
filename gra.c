#include "tgfit.h"
#include "plcdemos.h"
#include <time.h>

#define Exc_Flow_matrix Exc_Flow
//#define Exc_Flow_propagator  Exc_Flow

void Display_data(void) /*** DISPLAY DATA ***/
{
  int i,j,npts;
  float yi, yl,A,Sp,Ts,Ymax,Ymin;
  static PLFLT mtime[4096], ampl_m[4096], ampl_d[4096];

  
  Ts=0.2*ctrl.NData/4096; // Time shift
  Ts=0.0;
  Ymax=5.0; // Y scale Maximum 
  Sp=0.3*Ymax; // Separation 
  
  A=pow(10,(Ymax-Sp*(fit.nDatasets-1)));  // Data Max
  Ymin=log10(A*Set[1].data_min/Set[1].data_max); // Data Min
  //  Ymin=.0; 

  // pladv(0);
  plvpor(0.12, 0.95, 0.12, 0.9); 
  npts=Set[1].End-Set[1].Start;
  plwind(-Ts*fit.nDatasets, Set[1].TCal*ctrl.NData , Ymin, Ymax);

  plcol0(1);pllsty(1);
  plbox("bcnst", 0, 0, "bcvnstl", 0, 0);    
 
  for(i=1,yl=0,yi=Sp;i<=fit.nDatasets;i++,yl+=yi)
    {  
      npts=Set[i].End-Set[i].Start;
      for(j=1;j<=npts;j++)
	{
	  mtime[j-1]=model.XTime[Set[i].Start+j]-Ts*i;
	  if(Set[i].ymod[Set[i].Start+j]>Set[i].data_min)
	    ampl_m[j-1]=yl+log10(A*Set[i].ymod[Set[i].Start+j]/Set[i].data_max);
	  else
	    ampl_m[j-1]=yl+log10(A*(Set[i].data_min-0.4)/Set[i].data_max);
	  if(Set[i].data[Set[i].Start+j]>0.0)
	    ampl_d[j-1]=yl+log10(A*Set[i].data[Set[i].Start+j]/Set[i].data_max);
	  else
	    ampl_d[j-1]=yl+log10(A*(Set[i].data_min-0.4)/Set[i].data_max);
	}
      plcol0(1);
      plline(npts, mtime, ampl_m); // Plot model
      plcol0(i+1);
      plline(npts, mtime, ampl_d); // Plot data
    } 
  plcol0(2);
  pllab("Time, ns", "Amplitude, a.u.","Data and Fit");
  printf("\x1B[31m Right click on the graphics window to continue \x1B[0m\n");
  return;

}


void Display_trace(void) /*** DISPLAY DATA ***/
{
  int i,j,npts;
  float Ymax,Yinc;
  char  text[3];
  static PLFLT mtime[4096], ampl_m[4096];

  Yinc=1.0/fit.n;
  Ymax=0;
  for(i=1;i<=fit.n;i++)
    for(j=1;j<=ctrl.NData;j++)
      if(model.trace[i][j]>Ymax)Ymax=model.trace[i][j];
  
  // pladv(0); 
  plvpor(0.12, 0.95, 0.12, 0.9); 
  npts=ctrl.NData;
  plwind(0, Set[1].TCal*npts, 0, 1.05*Ymax);
  plcol0(1);pllsty(1);
  plbox("bcnst", 0, 0, "bcnst", 0, 0);    
  plschr(0,0.6);

  for(i=1;i<=fit.n;i++)
    {  
      for(j=1;j<=ctrl.NData;j++)
	{
	  mtime[j-1]=(j-1)*Set[1].TCal;
	  ampl_m[j-1] = model.trace[i][j];
	}
      sprintf( text, "%d", i );
      plcol0(i+1);
      plline(npts, mtime, ampl_m); // Plot model
     
      plmtex("rv", 2.0, Yinc*i - 0.5/(fit.n+1), 0, text );
    } 

  plcol0(2);
  plschr(0,1.0);
  pllab("Time, ns", "Amplitude, a.u.","Model kinetics");
  printf("\x1B[31m Right click on the graphics window to continue \x1B[0m\n");
  return;
}




void Display_residuals(void)
{
  int i,j,npts,n;
  float  Sp,Ymax,Ymin,Yinc;
  static PLFLT  ampl_m[4096], mtime[4096] ;
  static char Chi[50], filename[200], date[25];
  time_t tm1;

  tm1=time(NULL); 
  sprintf(date,"%s",&ctime(&tm1)[4]);
  date[20]='\0';
  sprintf(filename,"%s/%s",getenv("PWD"),fit.ans_in); 

  n=fit.nDatasets;
  Sp=1.0; // Separation
  Ymax=fit.nDatasets+Set[n].res_max/(Set[n].res_max-Set[n].res_min); // Y scale Maximum
  Ymin=1+Set[1].res_min/(Set[1].res_max-Set[1].res_min);
  Yinc=1.0/fit.nDatasets;

  // pladv(0);
  plvpor(0.12, 0.85, 0.12, 0.9); 
  npts=Set[1].End-Set[1].Start;
  plwind(0.0, Set[1].TCal*(ctrl.NData+2) , Ymin, Ymax);

  plcol0(1);pllsty(1);
  plbox("bcnst", 0.0, 0, "bcvnst", 0.0, 0);    
 
  for(i=1;i<=fit.nDatasets;i++)
    {
      sprintf(Chi,"#gx#u2#d=%.3f",Set[i].chisq);  
      npts=Set[i].End-Set[i].Start;
      for(j=1;j<=npts;j++){
	ampl_m[j-1]=i+Set[i].resid[j]/(Set[i].res_max-Set[i].res_min);
	mtime[j-1]=model.XTime[Set[i].Start+j];
      }
      plcol0(i+1);pllsty(1);
      plline(npts, mtime, ampl_m); // Plot model
      plcol0(1);
      pljoin(0,i,mtime[npts],i);
      plcol0(i+1);
      plschr(0,0.6);
      plmtex("rv", 3.0, Yinc*i - 0.5/(fit.nDatasets+1), 0, Chi );
    }  
  plcol0(2);
  plschr(0,0.5);
  plmtex("t", 1.5, 0.5, 0.5, filename);  
  plmtex("b", 7.0, 1.03, 0.5, date);
  plschr(0,1.0);
  pllab("Time, ns", "Amplitude, a.u.","Weighted Residuals");
  printf("\x1B[31m Right click on the graphics window to continue \x1B[0m\n");
  return;
}


void Display_DAS(void) /*** DISPLAY DAS ***/
{  
  time_t tm1;
  FILE *fp1;
  int  i, j, k, base, npts=100;
  static char leg_str_begin[50], filename[200], Chi[50], date[25];
  static PLFLT plp_wlg[max_datasets],plp_DAS[10][max_datasets],plp_wlg_spl[4096],plp_DAS_spl[10][4096];

  float x, 
    y, 
    Rho_0,
    ampl_min = 1000000.0, ampl_max =-100000.0,
    *wlg,
    *yield_n,
    **spl,
    **das_spl,
    **comp_yield,
    wlg_spl[npts];
  
  double *yield;
  float *spec, max_yield, max_spec;
  spec=vector(1,fit.nDatasets);

  base=849; // Font offset for plplot

 
  fp1=fopen("spec.dat","rt");
  if(fp1==NULL)nrerror("cannot open file spec.dat");
  for(i=1;i<=fit.nDatasets;i++)
     fscanf(fp1, "%f %f", &Set[i].wlg, &spec[i]);

  tm1=time(NULL); 
  sprintf(date,"%s",&ctime(&tm1)[4]);
  date[20]='\0';
  sprintf(filename,"%s/%s",getenv("PWD"),fit.ans_in); 
  sprintf(Chi,"#gx#u2#d=%.3f",fit.chisq/(fit.Glb_NData-fit.nGFitVarParams));
  
  wlg=vector(1,fit.nDatasets);
  yield=dvector(1,fit.nDatasets);
  yield_n=vector(1,fit.nDatasets);
  spl=matrix(1,fit.n,1,fit.nDatasets);
  das_spl=matrix(1,fit.n,1,npts);
  comp_yield=matrix(1,fit.n,1,fit.nDatasets);
  
  
  for(i=1;i<=fit.nDatasets;i++)
    for(j=1;j<=fit.n;j++)
      comp_yield[j][i]=.0;
  
  // Analytic normalization. Will not work for sequential models, 
  // will ignore rize components 
  /* printf("\n  W     Norm. yield\n");     */
  /* for(i=1;i<=fit.nDatasets;i++) */
  /*   { */
  /*     wlg[i]=Set[i].wlg; /\* ARRAY OF WAVELENGTH *\/ */
  /*     yield_n[i]=.0; */
  /*     yield[i]=.0; */
         
  /*    Rho_0=0.0; */
  /*     for(j=1;j<=fit.n;j++) */
  /* 	if(model.DAS[j][i] > 0.0) */
  /* 	  Rho_0+=model.DAS[j][i]; */
      
  /*     for(j=1;j<=fit.n;j++){ */
  /* 	if(model.DAS[j][i] > 0.0) { */
  /* 	  yield[i]+=model.DAS[j][i]/model.TransferMatrix[j][j]; /\* FLUORESCENCE YIELD *\/ */
  /* 	  yield_n[i]=yield[i]/Rho_0; /\* FLUORESCENCE YIELD *\/ */
  /* 	  comp_yield[j][i]=model.DAS[j][i]/(model.TransferMatrix[j][j]*Rho_0); */
	  
  /* 	} */
  /*     } */
  /*     printf("%5.0f%14g\n", wlg[i],yield_n[i]);  */
  /*   } */
  
  /* printf("\n             Normalized yield of components:\n  "); */
  /* for(j=1;j<=fit.n;j++) */
  /*   printf("%12i",j);    printf("\n");   */
  
  /* for(i=1;i<=fit.nDatasets;i++){ */
  /*   printf("%5.0f", wlg[i]); */
  /*   for(j=1;j<=fit.n;j++) */
  /*     printf("%12f",comp_yield[j][i]); */
  /*   printf("\n");   */
   /* } */

  // Numeric normalization 
  max_yield=max_spec=0.0f;
  for(j=1;j<=fit.nDatasets;j++){
    yield[j]=0.0f; 
    wlg[j]=Set[j].wlg; /* ARRAY OF WAVELENGTH */}
  for(j=1;j<=fit.nDatasets;j++){
    Expand_Params2Fit(Set[j].Params,Set[j].SAS_scale_factor,Set[j].data_max); 
    Exc_Flow(model.Fluor,4);
    // Prepare DAS matrix for display
    for(i=1;i<=fit.n;i++)
    // Unscale DAS 
      model.DAS[i][j]=fit.SAS[i]*2.4e4/Set[j].data_max;
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

  for(i=1;i<=fit.nDatasets;i++)
    {
      plp_wlg[i]=Set[i].wlg;
      for(j=1;j<=fit.n;j++){
	plp_DAS[j][i]=model.DAS[j][i];
	}
    }
 

  if(fit.nDatasets>2)
    {
      for(i=1;i<=fit.n;i++)
	if(ctrl.save_DAS[i]){ /* INTERPOLATE DAS USING CUBIC SPLINE */ 
	  spline(wlg,model.DAS[i],fit.nDatasets,spl[i]); 
	  for(j=0;j<npts;j++){ 
	    wlg_spl[j]=wlg[1]+(wlg[fit.nDatasets]-wlg[1])*j/npts;
	    splint(wlg,model.DAS[i],spl[i],fit.nDatasets,wlg_spl[j],&das_spl[i][j+1]);
	    plp_DAS_spl[i][j]=das_spl[i][j+1];
	    plp_wlg_spl[j]=wlg_spl[j];
	    if(das_spl[i][j+1]>ampl_max)ampl_max=das_spl[i][j+1];  /* FIND Y-SCALE */
	    if(das_spl[i][j+1]<ampl_min)ampl_min=das_spl[i][j+1];  
	  }  
	}
    }
  else    
    for(i=1;i<=fit.n;i++)
      for(j=1;j<=fit.nDatasets;j++)
	if(ctrl.save_DAS[i]){
	  if(model.DAS[i][j]>ampl_max)ampl_max=model.DAS[i][j];  /* FIND Y-SCALE */
	  if(model.DAS[i][j]<ampl_min)ampl_min=model.DAS[i][j];   
	}
  
  // pladv(0);
  plvpor(0.12, 0.72, 0.1, 0.9); 
  plwind(Set[1].wlg-1, Set[fit.nDatasets].wlg+1 , ampl_min-0.04*ampl_max, 1.04*ampl_max);
  
  plcol0(1);pllsty(1);
  plbox("bcnst", 0.0, 0, "bcvnst", 0.0, 0);
  pllsty(3);   
  pljoin(Set[1].wlg-1,0,Set[fit.nDatasets].wlg+1,0);
  plcol0(2);
  plmtex("b", 3.2, 0.5, 0.5, "Wavelength, nm");
  plmtex("l", 5.0, 0.5, 0.5, "Amplitude, a.u. ");
  plmtex("t", 3.0, 0.5, 0.5, "Decay-associated Spectra");
  plschr(0,0.5);
  plmtex("t", 1.5, 0.5, 0.5, filename);  
  plschr(0,1.0);
  // Draw symbols
  for(i=1;i<=fit.n;i++)
    if(ctrl.save_DAS[i])
      { 	
	plcol0(i+1);
	plsym(fit.nDatasets+1,plp_wlg,plp_DAS[i],base+i);
	pllsty(i);
	if(fit.nDatasets>2)
	  plline(npts,plp_wlg_spl,plp_DAS_spl[i]);
	if(fit.nDatasets==2)
	  pljoin(Set[1].wlg,model.DAS[i][1],Set[2].wlg, model.DAS[i][2]);
	
	x= Set[fit.nDatasets].wlg;
	y=0.7-0.08*i;
	plschr(0,0.8);
	sprintf(leg_str_begin,"#(%i) #gt#d%i#u = %6.3f ns", base+i,i, 1/model.TransferMatrix[i][i]);
	plmtex("rv", 2.0, 1-0.06*i, 0, leg_str_begin );
      }
	plcol0(2);
 	plmtex("rv", 5.0, 1-0.065*i, 0, Chi );
	plschr(0,0.5);
 	plmtex("b", 7.0, 1.3, 0.5, date);
	plschr(0,1.0);


  free_vector(wlg,1,fit.nDatasets);
  free_dvector(yield,1,fit.nDatasets);
  free_vector(yield_n,1,fit.nDatasets);
  free_matrix(spl,1,fit.n,1,fit.nDatasets);
  free_matrix(das_spl,1,fit.n,1,npts);
  free_matrix(comp_yield,1,fit.n,1,fit.nDatasets);
  printf("\x1B[31m Right click on the graphics window to continue \x1B[0m\n");
  return;
}

void Print(char *object)
{
  PLINT cur_strm, new_strm;
  if(!strncmp(object,"data",4))
    Display_data();
  if(!strncmp(object,"res",3))
    Display_residuals();
  if(!strncmp(object,"das",3))
    Display_DAS();
  plgstrm(&cur_strm);    /* get current stream */
  plmkstrm(&new_strm);   /* create a new one */ 
  
  plsfnam("plout.ps");    
  plsdev("ps");        
  
  plcpstrm(cur_strm, 0);  /* copy old stream parameters to new stream */
  plreplot();	           /* do the save by replaying the plot buffer */
  plend1();               /* finish the device */
  
  plsstrm(cur_strm);      /* return to previous stream */
  plbop();
  system("lpr plout.ps");
}
