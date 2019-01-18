#include <stdio.h>

//GLOBAL CONSTANTS
#define NVAR 7   //number of variables
#define NPAR 28   //number of parameters

static const int neq = NVAR;  //number of equations (same as NVAR)

//GLOBAL VARIABLES
static double par[NPAR];        //vector of model parameters
static double iniValue[NVAR];   //vector of initial variable values
static double iniTime;          //initial time

//FUNCTION PROTOTYPES (these three lines are very important)
int fex(int*, double*, double*, double*);

typedef int (*funptr)(int *neq, double *t, double *y, double *ydot);

extern void dlsoda_( funptr f, const int *NEQ,double *y,double *T,double *TOUT,int *ITOL,double *RTOL,double *ATOL,int *ITASK,int *ISTATE,int *IOPT,double *RWORK,int *LRW,int *IWORK,int *LIW,int *JAC,int *MF); //from libodepack.a (fortran codes)


//==============
//MAIN FUNCTION
//==============
int ODEmodel(double *par0, double *iniValue0, double iniTime0, int ntimepoints, double *timescale, double *sol)    
{
  double y[NVAR]; //variable values at certain time (used for both input and output)
  double atol[NVAR]; //absolute tolerance (scalar or array)
  double rtol = 1.E-2; //relative tolerance
  double tout; //time point of output
  int lrw = 134; //size of rwork (no less than 22 + NVAR * max(16, NVAR + 9))
  double rwork[134]; //real work array 
  int liw = 27; //size of iwork (no less than 20 + NVAR)
  int iwork[27]; //integer work array
  int itol = 2; //indicator for the type of error control: 1 if atol is scalar or 2 if atol is array 
  int itask = 1; //indicator specifying the task to be performed: 1 if normal
  int istate = 1; //indicator specifying the state of the calculation (used for both input and output): 1 if normal
  int iopt = 0; //indicator specifying optinal input: 0 if no optional inputs or 1 if one or more optional inputs
  int jt = 2; //indicator for Jacobian type (2 if a dummy argument is used)
  int jdum; //dummy argument
  int nsave = 0; //indicator for number of save in sol
  int i, iout; //cycle indicators

  //assign par[.], iniValue[.], and iniTime
  for (i = 0; i < NPAR; i++) par[i] = *(par0 + i);
  for (i = 0; i < NVAR; i++) iniValue[i] = *(iniValue0 + i); 
  iniTime = iniTime0;

  //assign initial variable values to y[.] used as input to dlsoda_()
  for (i = 0; i < NVAR; i++) y[i] = iniValue[i];

  //assign atol[.]
  for (i = 0; i < NVAR; i++) atol[i] = 1.E-2;
  
  //solve ODE and deposit solutions to sol[.]
  for (iout = 0; iout < ntimepoints; iout++) 
    {
      tout = *(timescale + iout);
      dlsoda_(fex,&neq,y,&iniTime,&tout,&itol,&rtol,atol,&itask,&istate,&iopt,rwork,&lrw,iwork,&liw,&jdum,&jt);
      if (istate <= 0) {printf("error istate = %d\n",istate); return 0;}
      for (i = 0; i < NVAR; i++)
        {
	  if (isnan(y[i])) {printf("solution has NaN values!\n"); return 0;}
	  *(sol + nsave) = y[i]; //y[.] used as output from dlsoda_()
	  nsave = nsave + 1;
        }         
    }  
  
  return 1;
}



//===================
//Define ODEs system
//===================
int fex(int *neq, double *t, double *y, double *ydot)
{
  //model parameters (size of NPAR)
  double s; 
  
  //initial variables (size of NVAR)
  double x1, x2, x3, x4, x5, x6, x7;
  
  //updated variables
  double d_x1, d_x2, d_x3, d_x4, d_x5, d_x6, d_x7;
  
  //assign parameter values
   s = par[0];  //SDF1;  par[1-27], parameters 
   
  //assign initial variable values
  x1 = y[0];
  x2 = y[1];  
  x3 = y[2];
  x4 = y[3];
  x5 = y[4];
  x6 = y[5];
  x7 = y[6];
  
  //update of varible values by ODEs system
  d_x1 = par[1]*s/(par[2]+s) - par[3]*x1;
  d_x2 = par[4]*s/(par[5]+s) - par[6]*x2;
  d_x3 = par[7]*x1/(par[8]+x1) - par[9]*x3;
  d_x4 = par[10]*x1/(par[11]+x1) + par[12]*x2/(par[13]/x2) - par[14]*x4;
  d_x5 = par[15]*x2/(par[16]+x2) - par[17]*x5;
  d_x6 = par[18]*x4/(par[19]+x4) + par[20]*x5/(par[21]+x5) - par[22]*x6;
  d_x7 = par[23]*x3/(par[24]+x3) + par[25]*x6/(par[26]+x6) - par[27]*x7;
  
  //output
  ydot[0] = d_x1;
  ydot[1] = d_x2;
  ydot[2] = d_x3;
  ydot[3] = d_x4;
  ydot[4] = d_x5;
  ydot[5] = d_x6;
  ydot[6] = d_x7;
  
  return 1;  
}

