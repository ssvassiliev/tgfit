#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <fcntl.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
//#include "mkl_vsl.h"

#define max_datasets  50
#define DTC 500 /* COLOR OF AXES */
#define MAX_STRING 80
#define FMAX(a,b) ( (a)>(b) ? (a) : (b))
#define FMIN(a,b) ((a)<(b) ? (a) : (b))
#define SIGN(a,b) ((b)>=0.0 ? fabs(a) : -fabs(a))

/********************* GLOBAL VARIABLES BEGIN  *********************/
enum gplot_colours {	BLACK, RED, BLUE, MAGENTA, GREEN, YELLOW, CYAN, WHITE, CORAL,
			RED_MAGENTA, GREEN_CYAN, BLUE_CYAN	};
struct CONTROLS {
  int save_DAS[20],
    max_iterrations,
    NData,
    NTAB,
    denom_param,
    SAS_positive,
    simulate,
    noise,
    convlv,
    weight,
    resid,
    autocorr,
    fit,
    correct,
    DAS,
    trace,
    trace_number,
    derv,
    derv_number,
    diag,
    verb; 
// 0 - no print
// 1 - print only chisq
// 2 - print chisq, constraints and timing  
  float max_lamda,
    tolerance,
    noise_factor;
} ctrl;

struct GLOBAL {
  float **Derivs,
    *Params,
    *SAS_scale_factor,
    *irf,
    *irf_spl,
    *spline,
    *ymod,
    *data,
    data_max,
    data_min,
    res_max,
    res_min,
    *auto_corr,
    *resid,
    irf_int,
    DW,
    chisq,
    wlg,
    TCal;
  int   *LnkList,
    *FitList,
    Start,
    End,
    nvarParams;
char label[100];
}Set[max_datasets];

struct FIT {
//  VSLConvTaskPtr task;
  float	**alpha,
    **covar,
    *oneda,
    *atry,
    *delta_params,
    *beta,
    *Bound,
    *GFitVect,
    *SAS,
    *SAS_scale_factor,
    *Fluor,
    *sig,
    *corr,
    chisq,
    old_chisq,
    delta_chisq,
    bkgnd,
    shift,
    delta_vect,
    lamda,
    auto_corr_min,
    auto_corr_max,
    resid_min,
    resid_max;
  int  nDatasets,
    n,
    nGFitParams,
    nGFitVarParams,
    nParams,
    *GFitList,
    Glb_NData;
  char ans_in[60],
    ans_out[60];
} fit;

struct MODEL {
  float **TransferMatrix,
    **TmpMatrix,
    **EigenvecMatrix,
    **EigenvecMatrixInverse,
    **trace,
    **DAS,
    *EigenvalRe,
    *Fluor,
    *XTime,
    wmax,
    wmin,
    d;
  int  *gindx;
} model;

struct NDIFF {
  float	**ymod_p,
    **ymod_n,
    **a,
    *u;
} ndiff;
/********************** GLOBAL VARIABLES END  **********************/

/***********************  PROTOTYPES BEGIN ************************/
/************************* tgfit ****************************/
void Load_Tbl(void);
void Load_Ans(char*);
void Load_Data(void);
void Save_Ans(char*);
void Grow_Ans (char *, char *, float); 
void Save_derivs(int);
void Save_results(void);
void Save_residuals(void);
void Save_autocorr(void);
void Save_fit(void);
void Save_Model(void);
void Save_trace(int);
void Save_DAS(void);
void Simulate(void);
void Parse_Links(float*,int);
void Expand_Lnk2Params(float*);
void Expand_Params2Fit(float*,float*,float);
void Exc_Flow(float*, int);
void Derivs(void);
void Derivs_diag(void);
void Chisquare(float*);
void Check_Constraints(void);
void Alloc_Memory(void);
void Display_data(void);
void Display_residuals(void);
void Display_DAS(void);
void draw_text(void);
void gplot_setnam(char *name, float f);
void gplot_setlab(char *name, char *text);
void Print(char *);
void Random_init_params(void);

/*********************** recipes ****************************/
void hqr(float**, int, float*, float*);
void svdcmp(float**, int, int, float*,float**);
void svbksb(float**,float*,float**,int,int,float*,float*);
void lubksb(float**,int,int*,float*);
void ludcmp(float**,int,int*,float*);
void elmhes(float**,int);
void four1(float*,int,int);
void realft(float*,int,int);
void twofft(float*,float*,float*,float*,int);
int  convlv(float*,int, float*, float*);
void convolve(int, float*);
void splint(float*,float*,float*,int,float,float*);
void spline(float*,float*,int,float*);
float *vector(int,int);
double *dvector(int,int);
float **matrix(int,int,int,int);
int *ivector(int,int);
void free_vector(float *,int,int);
void free_dector(double *,int,int);
void free_ivector(int*,int,int);
void free_matrix(float **, int,int,int,int);
void nrerror(char *);
/********************* modified recipes *********************/
void mrqmin(void);
void mrqcof(float*,float**,float*);
/************************* X-windows ************************/
void XDisplay_Residuals(Display*,Window,GC,unsigned long);
void XDisplay_Autocorr(Display*,Window,GC,unsigned long); 
void XDisplay_Data(Display*,Window,GC,int,unsigned long); 
/********************* PROTOTYPES END *********************/










