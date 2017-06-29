#include <math.h>
#include "tgfit.h" 
#include "mkl_vsl.h"

int convlv(float* data, int n, float* respns, float *ans)
{
 VSLConvTaskPtr task;
 int status;

 //  Create task descriptor 
 status = vslsConvNewTask1D(&task,VSL_CONV_MODE_FFT,n,n,n);
 if( status != VSL_STATUS_OK ){
   printf("ERROR: creation of job failed, exit with %d\n", status);
   return 1;
 }
 // Calculate 1 dimensional convolution of arrays h and x   
 status = vslsConvExec1D(task, &respns[1],1, &data[1],1, &ans[1],1);
 if( status != VSL_STATUS_OK ){
   printf("ERROR: job status bad, exit with %d\n", status);
   return 1;
 }
 //  Delete task descriptor 
 status = vslConvDeleteTask(&task);
 if( status != VSL_STATUS_OK ){
   printf("ERROR: failed to delete task object, exit with %d\n", status);
   return 1;
 }
 return 0;
}




