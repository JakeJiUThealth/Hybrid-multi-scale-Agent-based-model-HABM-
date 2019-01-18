#include <stdio.h>

//GLOBAL CONSTANTS
#define N_VAR 9   //number of variables
#define N_PAR 35   //number of parameters

static const int neq = N_VAR;  //number of equations (same as NVAR)

//GLOBAL VARIABLES
static double par2[N_PAR];        //vector of model parameters
static double iniValue2[N_VAR];   //vector of initial variable values
static double iniTime2;          //initial time

//FUNCTION PROTOTYPES (these three lines are very important)
int fex2(int*, double*, double*, double*);

typedef int (*funptr)(int *neq, double *t, double *y, double *ydot);

extern void dlsoda_( funptr f,const int *NEQ,double *y,double *T,double *TOUT,int *ITOL,double *RTOL,double *ATOL,int *ITASK,int *ISTATE,int *IOPT,double *RWORK,int *LRW,int *IWORK,int *LIW,int *JAC,int *MF); //from libodepack.a (fortran codes)


//==============
//MAIN FUNCTION
//==============
int ODEmodel2(double *par0, double *iniValue0, double iniTime0, int ntimepoints, double *timescale, double *sol)    
{
  double y[N_VAR]; //variable values at certain time (used for both input and output)
  double atol[N_VAR]; //absolute tolerance (scalar or array)
  double rtol = 1.E-2; //relative tolerance
  double tout; //time point of output
  int lrw = 184; //size of rwork (no less than 22 + NVAR * max(16, NVAR + 9))
  double rwork[184]; //real work array 
  int liw = 29; //size of iwork (no less than 20 + NVAR)
  int iwork[29]; //integer work array
  int itol = 2; //indicator for the type of error control: 1 if atol is scalar or 2 if atol is array 
  int itask = 1; //indicator specifying the task to be performed: 1 if normal
  int istate = 1; //indicator specifying the state of the calculation (used for both input and output): 1 if normal
  int iopt = 0; //indicator specifying optinal input: 0 if no optional inputs or 1 if one or more optional inputs
  int jt = 2; //indicator for Jacobian type (2 if a dummy argument is used)
  int jdum; //dummy argument
  int nsave = 0; //indicator for number of save in sol
  int i, iout; //cycle indicators

  //assign par[.], iniValue[.], and iniTime
  for (i = 0; i < N_PAR; i++) par2[i] = *(par0 + i);
  for (i = 0; i < N_VAR; i++) iniValue2[i] = *(iniValue0 + i); 
  iniTime2 = iniTime0;

  //assign initial variable values to y[.] used as input to dlsoda_()
  for (i = 0; i < N_VAR; i++) y[i] = iniValue2[i];

  //assign atol[.]
  for (i = 0; i < N_VAR; i++) atol[i] = 1.E-2;
  
  //solve ODE and deposit solutions to sol[.]
  for (iout = 0; iout < ntimepoints; iout++) 
    {
      tout = *(timescale + iout);
      dlsoda_(fex2,&neq,y,&iniTime2,&tout,&itol,&rtol,atol,&itask,&istate,&iopt,rwork,&lrw,iwork,&liw,&jdum,&jt);
      if (istate <= 0) {printf("error istate = %d\n",istate); return 0;}
      for (i = 0; i < N_VAR; i++)
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
int fex2(int *neq, double *t, double *y, double *ydot)
{
  //model parameters (size of NPAR)
  double s,d; 
  
  //initial variables (size of NVAR)
  double x1, x2, x3, x4, x5, x6, x7, x8 ,x9;
  
  //updated variables
  double d_x1, d_x2, d_x3, d_x4, d_x5, d_x6, d_x7, d_x8, d_x9;
  
  //assign parameter values
   s = par2[0];  //Stiffness;  par[1-27], parameters 
   d = par2[34];   //dose of BTZ
  //assign initial variable values
  x1 = y[0];
  x2 = y[1];  
  x3 = y[2];
  x4 = y[3];
  x5 = y[4];
  x6 = y[5];
  x7 = y[6];
  x8 = y[7];
  x9 = y[8];
  
  //update of varible values by ODEs system
  d_x1 = par2[1]*s/(par2[2]+s) - par2[3]*x1;
  d_x2 = par2[4]*x1/(par2[5]+x1) + par2[6]*x3/(par2[7]+x3) - par2[8]*x2;
  d_x3 = par2[9]*x1/(par2[10]+x1) - par2[11]*x3;
  d_x4 = par2[12]*x3/(par2[13]+x3) - par2[14]*x4;
  d_x5 = par2[15]*x2/(par2[16]+x2) - par2[17]*x5;
  d_x6 = par2[18]*x5/(par2[19]+x5) - par2[20]*x6;
  d_x7 = par2[21]*x2/(par2[22]+x2) + par2[23]*x4/(par2[24]+x4) - par2[25]*x7;
  d_x8 = par2[26]*x6/(par2[27]+x6) - par2[28]*x8; 
  d_x9 = par2[29]*x7/(par2[30]+x7) - par2[31]*x9 - par2[32]*d/(par2[33]+d);
  
  //output
  ydot[0] = d_x1;
  ydot[1] = d_x2;
  ydot[2] = d_x3;
  ydot[3] = d_x4;
  ydot[4] = d_x5;
  ydot[5] = d_x6;
  ydot[6] = d_x7;
  ydot[7] = d_x8;
  ydot[8] = d_x9;
  
  return 1;  
}

