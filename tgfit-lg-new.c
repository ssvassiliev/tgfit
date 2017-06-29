#include "tgfit.h"
#include "plcdemos.h"

int main( int argc, char *argv[])
{
  int i,j,k,h,NEW_FIT=1;
  FILE *fp;
  char  chs[80];


  // INITIALIZE PLPLOT 
  plsdev("xwin");
  plspage (600, 480, 600, 480, 100, 100);
  plinit();
  plfont(1);
  
  // PRINT WELCOME MESSAGE
  fprintf(stdout,"\n"
"\x1B[31m                     :->   TARGET GLOBAL DATA FIT   <-:\x1B[0m\n"
"\x1B[31m                         :->   VERSION 2.2.06   <-:\x1B[0m\n"
"           (c) 2000-2006 S.Vasil'ev, Brock University, Biology\n\n");
  if (argc<2){
    printf(
	   "SYNOPSIS\n" 
	   "        tgfit [INPUT.ans FILE] [OUTPUT.ans FILE - optional]\n\n" 
	   "FILES\n"
	   "       ~/.tgfitrc           - configuration file\n"
	   "       ./[your file].ans    - model setup file\n"
	   "       ./spec.dat           - fluorescence spectrum \n\n");
    return(0);
  }
    
 new_fit:
  // LOAD CONTROL PARAMETERS
  Load_Ans(argv[1]);
  Load_Tbl();

  // LOAD ANSWER FILE
  Load_Ans(argv[1]);


  printf("conv=%i\n", ctrl.convlv);

  // LOAD DATA
  if(NEW_FIT)Load_Data();

  // PROCESS LINKS
  Parse_Links(fit.GFitVect,1);

  // SET NAME OF THE OUTPUT FILE
  strcpy(fit.ans_in,argv[1]);
  if(argc==3)strcpy(fit.ans_out,argv[2]);
  if(argc==2)strcpy(fit.ans_out,argv[1]);

  // ALLOCATE MEMORY
  if(NEW_FIT)Alloc_Memory();
  
  // CLEAN BUFFERS
  for(i=ctrl.NData+1,h=1;i<=2*ctrl.NData;i++,h++){
    model.XTime[h]=(h-1)*Set[1].TCal;
    model.Fluor[i]=0.0;
    for(j=1;j<=fit.nDatasets;j++){
      Set[j].irf[i]=0.0; 
      Set[j].irf_spl[i]=0.0;
    }
  }
  
  // LOAD CORRECTION DATA
  if(ctrl.correct){
    fp=fopen("corr.dat","rt");
    if(fp==NULL)nrerror("correction file could not be opened for writing");  
    for(i=1;i<=ctrl.NData;i++)
      fscanf(fp, "%f",  &fit.corr[i]);
    fclose(fp);
  }
  
  // FIND DATA MINIMUM AND MAXIMUM 
  for(j=1;j<=fit.nDatasets;j++){
    for(i=1;i<=ctrl.NData;i++){
      if(ctrl.correct)Set[j].data[i]=Set[j].data[i]*fit.corr[i];
      Set[j].data_max=FMAX(Set[j].data_max,Set[j].data[i]);
      if(Set[j].data[i]!=.0)
	Set[j].data_min=FMIN(Set[j].data_min,Set[j].data[i]);
    }
    spline(model.XTime,Set[j].irf,ctrl.NData,Set[j].spline);
  }
  printf("Data peak counts:\n" );
  for(j=1,k=0;j<=fit.nDatasets;j++)
    {
      printf("%7i",(int)Set[j].data_max);
      if(k++==4)
	{
	  printf("\n");
	  k=0;
	}
    }
  printf("\n\x1B[31m** I will scale decay amplitudes using peak count values \x1B[0m\n"
	 "** Use \x1B[31m SAS_scale_factor \x1B[0m for functional links: SAS[i] = k * SAS[j]  \x1B[0m\n\n");
  pladv(0);

  // DO FITTING OR SIMULATION
  if(ctrl.simulate)
    {
      Simulate();
      Display_data();
      //     if(ctrl.trace)Save_trace(ctrl.trace_number);
      printf("simulation completed press < Return > to exit\n");
      getchar(); 
      return(1);
    }  
  else 
    {
      fit.lamda=-1;
    more:
      for(h=1;h<=ctrl.max_iterrations;h++)
	{
	  fprintf(stdout,"\r%i ",h);fflush(stdout);
	  mrqmin();    
	  if((fit.delta_chisq<ctrl.tolerance)&&(fit.delta_chisq!=0)){
	    printf("\nTerminated - minimal allowed delta chisquare reached\n");
	    break;
	  }
	  if(fit.lamda>ctrl.max_lamda){
	    printf("\nTerminated - max lamda reached\n");
	    break;
	  }
	  if(fit.delta_vect==0.0){
	    printf("\nTerminated - minimal allowed delta fit vector reached\n");
	    break;
	  }
	  if(h==ctrl.max_iterrations)
	    printf("\nTerminated - max iterrations reached\n");
	}
      
      // SAVE RESULTS

      if(ctrl.resid)Save_residuals();
      if(ctrl.autocorr)Save_autocorr();
      if(ctrl.fit)Save_fit();
      if(ctrl.DAS)Save_DAS();
      if(ctrl.trace)Save_trace(ctrl.trace_number);
      if(ctrl.derv)Save_derivs(ctrl.derv_number);
      printf("Final chisquare:\n");
      Save_Ans(fit.ans_out);
      printf("\n");
      
      NEW_FIT=0;
      chs[0]='2';
      
      // READ COMMANDS TO CHANGE DISPLAY MODE OR DO MORE FITTING
      while(1) {
	if(!strncmp(chs,"?",1)){
	  printf(
		 " 1 - display data\n"
		 " 2 - display residuals\n"
		 " 3 - display DAS\n"
		 " it - more ITerrations\n" 
		 " fi - FIt from scratch\n" 
		 " print data \n"
		 " print das \n"
		 " print res \n"
		 " !  - lines starting with a ! are dispatched to the OS\n"
		 " ex - EXit\n"
		 " ? - show options\n"
		 " NOTE: to calculate fl. yield use option 3\n");goto done;}
	if(!strncmp(chs,"it",2))
	  goto more;
	if(!strncmp(chs,"fi",2))
	  goto new_fit;
	if(!strncmp(chs,"!",1))
	  {system(&chs[1]);goto done;}
	if(!strncmp(chs,"print data",10))
	  {
	    Print("data");
	    goto done;
	  }
	if(!strncmp(chs,"print das",9))
	  {
	    Print("das");
	    goto done;
	  }
	if(!strncmp(chs,"print res",9))
	  {
	    Print("res");
	    goto done;
	  }

	if(!strncmp(chs,"1",1))
	  {Display_data();pladv(0);goto done;} 
	if(!strncmp(chs,"2",1))
	  {Display_residuals();pladv(0); goto done;}
	if(!strncmp(chs,"3",1))
	  {Display_DAS(); pladv(0);goto done;}
	if(!strncmp(chs,"ex",2))
	  break;
	printf(" Unknown command\n Type < ? > for list of options\n");
      done:
	printf("TgFit>");
	fgets(chs,MAX_STRING,stdin);
      }
    }
  plend();
  return(1);
}       


void Alloc_Memory(void)
{
int i;
  model.Fluor=vector(1,2*ctrl.NData);
  model.TmpMatrix=matrix(1,fit.n,1,fit.n);
  model.EigenvecMatrix=matrix(1,fit.n,1,fit.n);
  model.EigenvecMatrixInverse=matrix(1,fit.n,1,fit.n);
  model.gindx=ivector(1,fit.nGFitParams);

  model.EigenvalRe=vector(1,fit.n);
  model.trace=matrix(1,fit.n,1,2*ctrl.NData);
  model.DAS=matrix(1,fit.n,1,fit.nDatasets);
  model.XTime=vector(1,ctrl.NData);
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
    Set[i].ymod=vector(1,ctrl.NData);
    Set[i].Derivs=matrix(1,fit.nGFitParams,1,ctrl.NData);
    Set[i].auto_corr=vector(1,(Set[i].End-Set[i].Start)/2);
    Set[i].resid=vector(1,ctrl.NData);
    Set[i].irf_spl=vector(1,2*ctrl.NData);
    Set[i].spline=vector(1,ctrl.NData);
  } 
}


void Exc_Flow_m(float* buf)
{
  // Euler
  int i, j, p;
  double v[fit.n];
  char UPLO = 'L';
  char TRANS = 'N';
  int LDA = fit.n;
  int INCX = 1, INCY = 1;
  double BETA = 1.0;
  double dt = 0.001;
  double k;
  int steps = 10;
  float M[fit.n][fit.n];


  for(i=1; i <= fit.n; i++){ /*  RECIPES -> LAPACK  */
    for(j = 1; j <= fit.n; j++){
      M[i-1][j-1] = -1*model.TransferMatrix[j][i]; /*TRANSPOSE MATRIX */
      /* LAPACK ROUTINES REQUIRE NOTATION [column][row] */
    }
  }
    
  /* for(i=0; i < fit.n; i++){ */
  /*   for(j = 0; j < fit.n; j++){ */
  /*     printf("%f ", M[i][j]);  */
  /*   } */
  /*   printf("\n"); */
  /* } */
    



  for(i = 1; i <= ctrl.NData; i++) {  /* RESET TRACE MATRIX */
    buf[i]=0.0;
    for(j=1; j<=fit.n; j++)
      model.trace[j][i]=0.0;
  }


  for (i = 0; i < fit.n; i++) {
    v[i] = fit.Bound[i];
    model.trace[i+1][1] = v[i];
    if(fit.SAS[i+1])
      buf[i+1] += fit.SAS[i+1]*model.trace[i+1][1];
  }
    
  dt = Set[1].TCal / (double)steps;

  /* for (i = 0; i < fit.n; i++) */
  /*     printf("%f ",v[i]);	   */
  /*   printf("\n"); */

  for (i = 2; i <= ctrl.NData; i++) 
    {
      for (p = 0; p < steps; p++)            
	dgemv_ (&TRANS, &fit.n, &fit.n, &dt, M, &LDA, &v, &INCX, &BETA, &v, &INCY);

      /* for (i = 0; i < fit.n; i++) */
      /*   printf("%f ",v[i]);	   */
      /* printf("\n"); */

      for (j = 1; j <= fit.n; j++) 
	{ //Copy v into trace
	  model.trace[j][i] = v[j-1];
	  //	printf("%f ",model.trace[j][i]);
	  if(fit.SAS[j])
	    buf[i] += fit.SAS[j]*model.trace[j][i];
	}
    }
}

 

void Exc_Flow(float* buf) /*** COMPUTE MODEL DECAY CURVES ***/
     /* Solves system of rate equations, computes fluorescence decay and
	stores it in vector model.Fluor pointed.
	User supply model parameters in struct MODEL */
{
  extern void sgeev_();
  extern void sgetrf_();
  extern void sgetri_();
  float TmpMatrix[fit.n][fit.n],diag=1.;
  char jobvl = 'N', jobvr = 'V';
  int info, lwork = 4*fit.n;
  float vl[fit.n][fit.n], vr[fit.n][fit.n],
    wi[fit.n], wr[fit.n],
    work[lwork], ipiv[fit.n];
  float S[fit.n+1][fit.n+1],A[fit.n+1],Rho_0[fit.n+1],u,du;  

  int  i,j,k,l,m;

  for(i=1;i<=ctrl.NData;i++){  /* RESET TRACE MATRIX */
    buf[i]=0.0;
    for(j=1;j<=fit.n;j++)
      model.trace[j][i]=0.0;
  } 
 

  for(i=1;i<=fit.n;i++){ /*  RECIPES -> LAPACK  */
    for(j=1;j<=fit.n;j++){
      if((i != j) && (model.TransferMatrix[i][j] != .0))diag=0;
      TmpMatrix[i-1][j-1]=model.TransferMatrix[j][i]; /*TRANSPOSE MATRIX */
      S[i][j]=0.0;
      /* LAPACK ROUTINES REQUIRE NOTATION [column][row] */
    } 
  } 

  /* IF MATRIX DIAGONAL  */
  if(diag){  
    for(i=1;i<=ctrl.NData;i++)
      for(j=1;j<=fit.n;j++){
	model.trace[j][i]+=exp(-model.TransferMatrix[j][j]*model.XTime[i]);
	if(fit.SAS[j])buf[i]+=fit.SAS[j]*fit.Bound[j]*model.trace[j][i];
      }
return;
  }

  /* MATRIX FACTORIZATION  */
  sgeev_(&jobvl, &jobvr, &fit.n, &TmpMatrix, &fit.n, &wr, &wi, &vl, &fit.n, &vr, &fit.n, &work, &lwork, &info);
  if(info)nrerror("Eigensystem solution failed");   /* FIND EIGENVALUES AND RIGHT EIGENVECTORS */
  for(i=0;i<fit.n;i++)  /* MAKE A COPY OF RIGHT EIGENVECTORS */
    for(j=0;j<fit.n;j++)vl[j][i]=vr[j][i]; 
  sgetrf_(&fit.n, &fit.n, &vl, &fit.n, &ipiv, &info);  /* COMPUTE LU FACTORIZATION OF THE MATRIX VR */
  if(info)nrerror("Matrix inversion failed (1)");  
  sgetri_(&fit.n, &vl, &fit.n, &ipiv, &work, &lwork, &info); /* COMPUTE INVERSE MATRIX: RESULT IN VL */  
  if(info)nrerror("Matrix inversion failed (2)");
 
 /* LAPACK -> RECIPES */
/*    for(i=0;i<fit.n;i++){  */
/*      model.EigenvalRe[i+1]=wr[i]; */
/*      for(j=0;j<fit.n;j++){ */
/*        model.EigenvecMatrix[i+1][j+1]=vr[j][i];  */
/*        model.EigenvecMatrixInverse[i+1][j+1]=vl[j][i];  */
/*      } */
/*    } */

/* CALCULATE EXCITATION PROBABILITY FLOW */
/*    for(l=1;l<=ctrl.NData;l++)  */
/*      { */
/*        for(i=1;i<=fit.n;i++) */
/*  	for (k=1;k<=fit.n;k++) */
/*  	  for(j=1;j<=fit.n;j++) */
/*  	    { */
/*  	      if(fit.Bound[j]) */
/*  		model.trace[i][l]+=model.EigenvecMatrix[k][i]* */
/*  		  exp(-model.EigenvalRe[k]*model.XTime[l])* */
/*  		  model.EigenvecMatrixInverse[j][k]*fit.Bound[j]; */
/*  	    } */
/*        for(m=1;m<=fit.n;m++) */
/*  	if(fit.SAS[m])buf[l]+=fit.SAS[m]*model.trace[m][l]; */
/*      } */

  for(i=0;i<fit.n;i++){ 
    model.EigenvalRe[i+1]=wr[i];
    for(j=0;j<fit.n;j++){
      model.EigenvecMatrix[i+1][j+1]=vr[i][j]; 
      model.EigenvecMatrixInverse[i+1][j+1]=vl[i][j]; 
    }
  }

  for(k=0;k<=fit.n;k++){Rho_0[k]= A[k]=0.0;}

    for(k=1;k<=fit.n;k++){
      u = 0;
      for(j=1;j<=fit.n;j++)
	u += model.EigenvecMatrixInverse[j][k]*fit.Bound[j];
      S[k][i]=0;
      for(i=1;i<=fit.n;i++)       
	S[k][i]+= model.EigenvecMatrix[k][i]*u;
    }
    
    for(i=1;i<=fit.n;i++){
      du = exp(-model.EigenvalRe[i]*Set[1].TCal);
      for(k=1;k<=fit.n;k++){
	u = S[i][k];
	model.trace[k][1] += u;
	for(l=2;l<=ctrl.NData;l++){
	  u *= du;
	  model.trace[k][l] += u;
	}
      }
    }  
    for(l=1;l<=ctrl.NData;l++)
      for(m=1;m<=fit.n;m++)
	if(fit.SAS[m])buf[l]+=fit.SAS[m]*model.trace[m][l];

/*    for(k=1;k<=fit.n;k++){ */
/*      for(j=1;j<=fit.n;j++) */
/*        A[k]+=model.EigenvecMatrixInverse[j][k]*fit.Bound[j]; */
/*      for(j=1;j<=fit.n;j++){ */
/*        S[k][j] = fit.SAS[j] * model.EigenvecMatrix[k][j] * A[k]; */
/*        Rho_0[k] += S[k][j]; */
/*      } */
/*      du=exp(-model.EigenvalRe[k]*Set[1].TCal); */
/*      u = Rho_0[k]; */
/*      buf[1] += u; */
/*      for(l=2;l<=ctrl.NData;l++){ */
/*        u=u*du; */
/*        buf[l] += u; */
/*      } */
/*    } */
    
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
  convlv(model.Fluor,2*ctrl.NData,Set[ds].irf_spl,1,fit.Fluor);
   for(i=1;i<=ctrl.NData;i++)buf[i]=fit.Fluor[i]/Set[ds].irf_int+fit.bkgnd;
}

void Expand_Params2Fit(float *params, float *sas_scaling, float d_max)
{
  int i,j,k;
  /* EXPAND LIST OF PARAMETERS TO STRUCT FIT */
  for(i=1;i<=fit.n;i++)fit.Bound[i]=params[i];
  for(j=1;j<=fit.n;j++,i++)fit.SAS[j]=sas_scaling[j]*params[i]*d_max/2.4e4;
  for(k=1;k<=fit.n;k++) /* Fill off-diagonal */
    for(j=1;j<=fit.n;j++,i++)model.TransferMatrix[k][j]=params[i];
  /*	for(k=1;k<=fit.n;k++){
	for(j=1;j<=fit.n;j++)printf("%f  ", model.TransferMatrix[k][j]);
	printf("\n");
	}  */
  for(k=1,i=1;k<=fit.n;k++,i++){  /* Fill diagonal */
    model.TransferMatrix[k][k]*=2;
    for(j=1;j<=fit.n;j++)model.TransferMatrix[k][k]-=params[(j+1)*fit.n+i];
  }

  fit.bkgnd=params[fit.nParams-1];
  fit.shift=params[fit.nParams];

/*    printf("\n); */
/*    for(i=1;i<=fit.n;i++){ */
/*      for(j=1;j<=fit.n;j++)printf("%f  ", model.TransferMatrix[i][j]); */
/*      printf("\n"); */
/*      } */
/*    getchar(); */
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
}


void Expand_Lnk2Params(float *a)
{
  /*** EXPAND GLOBAL FIT PARAMETER VECTOR TO ARRAY DATASET[].PARAMETERS[] ***/
  int i,k;
  for(k=1;k<=fit.nDatasets;k++)
    for(i=1;i<=fit.nParams;i++)
      if(Set[k].FitList[i])Set[k].Params[i]=a[Set[k].LnkList[i]];
}


void Chisquare(float* a) /*** COMPUTE CHISQUARE ***/
  /* Given global parameter vector *a this function evaluates global
     chisquare and displays chisquares of each dataset */
{
  int i,ds;
  float dy,tmp_chisq;
  
  Expand_Lnk2Params(a);
  fprintf(stdout,"* "); 
  fit.chisq=0.0;
  tmp_chisq=0.0;

  for(ds=1;ds<=fit.nDatasets;ds++)
    {
      Set[ds].res_max=-1e8;
      Set[ds].res_min=1e8;
    Expand_Params2Fit(Set[ds].Params,Set[ds].SAS_scale_factor,Set[ds].data_max);
    Exc_Flow(model.Fluor);
    
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
    fprintf(stdout,"%.3f ",Set[ds].chisq);
    tmp_chisq=0.0;
  }
 fprintf(stdout,"*  %.3f [%.3f] L=%.0e",
	 fit.chisq/(fit.Glb_NData-fit.nGFitVarParams),
	 fit.old_chisq/(fit.Glb_NData-fit.nGFitVarParams), fit.lamda);
 fflush(stdout);
}


void Check_Constraints(void)
{
  int i,j,k,ds;
  
  for(ds=1;ds<=fit.nDatasets;ds++){
    for(i=1;i<=fit.n;i++)  /* check bond */
      if(Set[ds].Params[i]<0){
	Set[ds].Params[i]=-Set[ds].Params[i]/5;
	printf("     ** C ** %i-Bond\n",ds);
      }
    
    if(ctrl.SAS_positive)
      {
	for(;i<=2*fit.n;i++)  /* check SAS */
	  if(Set[ds].Params[i]<0){
	    Set[ds].Params[i]=-Set[ds].Params[i]/5;
	    printf("     ** C ** %i-SAS \n",ds);
	  }
      }
	  for(;i<=2*fit.n;i++);    
    
    for(k=1;k<=fit.n;k++) /* check off-diagonal matrix elements*/
      for(j=1;j<=fit.n;j++,i++)
	if(k!=j&&Set[ds].Params[i]>0){
	  Set[ds].Params[i]=-Set[ds].Params[i]/5;
	  printf("     ** C ** %i-MEL\n",ds);
	}
    
    for(i=1;i<=fit.n;i++) /* check diagonal matrix elements*/
      if(Set[ds].Params[(i+1)*fit.n+i]<0){
	Set[ds].Params[(i+1)*fit.n+i]=-Set[ds].Params[(i+1)*fit.n+i];
	printf("     ** C ** %i-DMEL\n",ds);
      }
  }
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
    Exc_Flow(model.Fluor);
    
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


































































