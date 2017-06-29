#include "tgfit.h"
//#include <plcdemos.h>
#include <plConfig.h>
#include <plplot.h>




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
	   "       ./spec.dat           - fluorescence spectrum \n\n"
	   "** Order of matrix elements: **\n"
	   "| K11 K12 |\n"
	   "| K21 K22 |\n");
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
  ctrl.verb=1; // allow to print timing 
  
  // CLEAN BUFFERS
  for(i=ctrl.NData+1,h=1;i<=5*ctrl.NData;i++,h++){
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
      if(ctrl.trace)Save_trace(ctrl.trace_number);
      
      //    Display_data();
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
	  fprintf(stdout,"\n%i ",h);fflush(stdout);

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
		 " 4 - display model kinetics\n"
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
	if(!strncmp(chs,"4",1))
	  {Display_trace(); pladv(0);goto done;}
	if(!strncmp(chs,"q",1))
	  {Save_Model();break;}
	printf(" Unknown command\n Type < ? > for list of options\n");
      done:
	printf("TgFit>");
	fgets(chs,MAX_STRING,stdin);
      }
    }

  plend();
  return(1);
}       

