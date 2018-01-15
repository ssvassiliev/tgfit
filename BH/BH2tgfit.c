#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <ctype.h>
#include <string.h>

#define OPTIONS       "f:o:n:h"
#define default_spcm_datafile "test.sdt"
#define default_tgfitfile_prefix "block_"
#define INFO_LENGTH 5000
#define LABEL_LENGTH 100
#define MAX_NCURVES 100

int main(int argc, char *argv[])
{
  FILE *tgfitfile;
  int i, j, opt, spcfile, no_of_irf_block, n_label_tgfile;
  float tcal;
  unsigned short databuf[MAX_NCURVES][4096];
  float databuff[MAX_NCURVES][4096];
  char info[INFO_LENGTH];
  float DATA_MAX, PEAK_COUNT;
  int normalize_flag;
 

  char label_tgfitfile[LABEL_LENGTH];
  char label_spcfile[LABEL_LENGTH];
  char peak[20];


  extern char *optarg;           /* option argument from getopt */ 
  static char usage[]="\nThe program to convert Becker & Hickl SPCM .sdt files into format required"
"by TGFIT target global fit data analysis software\n\n"
"Options:\n"
"  -f [FILE]         SPCM data file, default -  test.sdt\n"
"  -o [FILE]         prefix for output TGFIT files, default - block\n"
"  -n [PEAK COUNT]   normalize data,  all datasets will have the same counts in the peak channel\n";

  /* SPC file HEADER */
  short revision;
  long info_offset;
  short info_length;
  long setup_offs;
  short setup_length;
  long data_block_offset;
  short no_of_data_blocks;
  long data_block_length;
  long meas_desc_block_offset;
  short no_of_meas_desc_blocks;
  short meas_desc_block_length;
  unsigned short header_valid;
  unsigned long reserved1;
  unsigned short reserved2;
  unsigned short chksum;

  /* Data block header */
  short block_no[MAX_NCURVES];
  long data_offs;
  long next_block_offs;
  unsigned short block_type;
  short meas_desc_block_no;
  unsigned long lblock_no;
  unsigned long block_length;

  /* Measurement description block headers */
  short adc_re;
  float tac_r;
  short tac_g;
  float coll_time[MAX_NCURVES];
  short stopt[MAX_NCURVES];

  /* defaults */
  strcpy(label_spcfile, default_spcm_datafile);
  strcpy(label_tgfitfile, default_tgfitfile_prefix);
  normalize_flag=0;

 while ((opt = getopt(argc, argv, OPTIONS)) != EOF) {
   switch (opt) {
   case 'f':
     strcpy(label_spcfile,optarg);
       break;
   case 'o':
     strcpy(label_tgfitfile,optarg);
     break;
   case 'n':
     normalize_flag=1;
     sscanf(optarg,"%f",&PEAK_COUNT);
     break;
   case 'h':
     fprintf(stderr,"%s\n", usage);
     exit(1);
   }  /* end switch */
 } /* end while */


 n_label_tgfile=strlen(label_tgfitfile);

if((spcfile = open(label_spcfile, O_RDONLY )) == -1)  {
  fprintf(stderr,"Error: Could not find SPCM data file %s\n",  label_spcfile);
  _exit(3);
}
  /* Read file HEADER */
  read(spcfile,&revision,2);
  read(spcfile,&info_offset,4);
  read(spcfile,&info_length,2);
  read(spcfile,&setup_offs,4);
  read(spcfile,&setup_length,2);
  read(spcfile,&data_block_offset,4);
  read(spcfile,&no_of_data_blocks,2);
  read(spcfile,&data_block_length,4);
  read(spcfile,&meas_desc_block_offset,4);
  read(spcfile,&no_of_meas_desc_blocks,2);
  read(spcfile,&meas_desc_block_length,2);

  /* Read file INFO */
  lseek(spcfile,info_offset,SEEK_SET);
  read(spcfile,&info,info_length);
  for(i=0;i<info_length;i++)
    if((info[i]==13)|(info[i]==4))
      info[i]=32;

  /* Read Measurement description block */
  lseek(spcfile,meas_desc_block_offset+64,SEEK_SET);
  read(spcfile,&tac_r,4);
  read(spcfile,&tac_g,2);
  lseek(spcfile,meas_desc_block_offset+82,SEEK_SET);
  read(spcfile,&adc_re,2);
  tcal=1e9*tac_r/(tac_g*adc_re);


/*   for(j=0,i=meas_desc_block_offset;j<no_of_meas_desc_blocks;j++,i+=meas_desc_block_length) */
/*     {  */
/*       lseek(spcfile,i+92,SEEK_SET); */
/*       read(spcfile, &coll_time[j],4); */
/*       printf("Collection time: %f\n",coll_time[j]); */
/*     } */


  printf("  File info:\n", info);
  printf("%s", info);
  printf("  Number of datablocks = %i\n", no_of_data_blocks);
  printf("  Number of measurement description blocks = %i \n", no_of_meas_desc_blocks);  
  printf("  Data block length  = %li\n", data_block_length);
  printf("  TAC resolution = %g ns\n", tac_r*1e9);
  printf("  TAC gain = %i\n", tac_g);
  printf("  ADC resolution = %i\n", adc_re);
  printf("  Tcal = %g ns\n",tcal);
     
  if(normalize_flag)
    printf("\nNormalization is ON, Peak Count will be %.0f\n", PEAK_COUNT);  
  else
    printf("\nNormalization is OFF");  

  printf("\n  Number of the measurement page containing IRF: ");
  scanf("%i",&no_of_irf_block);


  /* Read data */
  lseek(spcfile,data_block_offset,SEEK_SET);
  int  irfnum;
  for(j=0;j<no_of_data_blocks;j++)
    {
      /* Data block header */
      read(spcfile,&block_no[j],2);
      if(block_no[j]==no_of_irf_block)
        irfnum=j;   
      read(spcfile,&data_offs,4);
      read(spcfile,&next_block_offs,4);
      read(spcfile,&block_type,2);
      read(spcfile,&meas_desc_block_no,2);
      read(spcfile,&lblock_no,4);
      read(spcfile,&block_length,4);

      /* Data */
      for(i=0;i<data_block_length/2;i++)
      read(spcfile,&databuf[j][i],2);  
    }


  // NORMALIZE DATA
      printf("Multiplying datasets by:\n");
      for(j=0;j<no_of_data_blocks;j++)
	{
	  // FIND DATA MAXIMUM
	  DATA_MAX=0;
	  for(i=0;i<data_block_length/2;i++)
	    if(databuf[j][i]>DATA_MAX)
	      DATA_MAX=databuf[j][i];
	  DATA_MAX=PEAK_COUNT/DATA_MAX;
	  if(!normalize_flag)
	    DATA_MAX=1;
	  // NORMALIZE DATA
	  printf(" %f\n",DATA_MAX);
	  for(i=0;i<data_block_length/2;i++)
	    databuff[j][i]= DATA_MAX*(float)databuf[j][i];
	  // NORMALIZE DATA
	} 
    
  

    printf("Block order:");
    for(j=0;j<no_of_data_blocks;j++)printf("%i ",block_no[j]);
    printf("\n");

  for(j=0;j<no_of_data_blocks;j++)
    if(j!=irfnum){
      sprintf(&label_tgfitfile[n_label_tgfile],"%i.asc",block_no[j]);
      tgfitfile=fopen(label_tgfitfile,"w");
      fprintf(tgfitfile,"%f     1    1000\n", tcal);
      for(i=0;i<data_block_length/2;i++)
	fprintf(tgfitfile,"%12.2f %12.2f\n",(float)databuff[irfnum][i], (float)databuff[j][i]);
      fclose(tgfitfile); 
    }
  close(spcfile);
  return(0);
}
