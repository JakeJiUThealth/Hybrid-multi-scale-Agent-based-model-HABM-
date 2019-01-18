//#include "StdAfx.h"
//#include ".\microenv.h"
#include "MicroEnv.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <fstream>
#include <string.h>
#include <iostream>
#include <math.h> 
#include "par.hpp" 

//--should be comment selection in Linux----------------------------------------------
//#include "mclmcr.h"
//#include "matrix.h"
//#include "mclcppclass.h"
//#include "libODEs.h"

//#pragma comment(lib,"mclmcrrt.lib")
//#pragma comment(lib,"libmx.lib")
//#pragma comment(lib,"libmat.lib")
//#pragma comment(lib,"mclmcr.lib")
//#pragma comment(lib,"libODEs.lib")
//-----------------------------------------------------

#define NVAR 7 //number of variables in the model
#define NPAR 28 //number of parameters in the model

#define N_VAR 9 //number of variables in the model
#define N_PAR 35 //number of parameters in the model
//------------------------------------------------
// PARAMETERS FOR PROLIFERATION AND MIGRATION
#define K_s 3.2e-5;
#define r_r0 0.58;  //>=0.525 for MIC (may important)
#define r_a 0.0;
#define r_d0 0.42;  //should consistent with r_r0
#define r_pathway 0.15;  
#define FE0 0.5;
#define FEMAX 0.5;

#define DMIC 2.0;
#define DMM 2.0;
#define DCTL 2.0;
#define DTREG 2.0;
#define DTMM 2.0;
#define DPC 2.0;

#define E0 100;
#define EMAX 400;
#define IV 1.0;
#define mux 1.0;
#define muy 1.0;
#define muz 1.0;
#define theta 1.0;


#define PI 3.14159265358979323846
//-------------------------------------------------------------
//using namespace par;

MicroEnv::MicroEnv()
{  
    
     time_t t = time(0);
     struct tm * now =localtime(&t);
     stringstream ss;
     ss<<par::iT<<now->tm_mday<<now->tm_hour<<now->tm_min;
     string str=ss.str();

     int  seed = atoi(str.c_str());     
     srand(seed);
//-------------------------------------------------------------
     TGFb = new double **[GridLength+1];
     IL6 = new double **[GridLength+1];
     tmp= new double **[GridLength+1];
     SDF1 = new double **[GridLength+1];
     STIFF = new double **[GridLength+1];
     ECMHigh = new double **[GridLength+1];
     ECMCan = new double **[GridLength+1];
     ECM3 = new double **[GridLength+1];
     Protain =new double **[GridLength+1];
    
     string sFilename1="";
     sFilename1 = par::sOutDir + "/SDF/SDF1_" + string(itostr(par::iT)) + ".txt";
     foutsdf1.open(sFilename1,ios::app);

     string sFilename2="";
     sFilename2 = par::sOutDir + "/STIFF/STIFFNESS_" + string(itostr(par::iT)) + ".txt";
     foutstiff.open(sFilename2,ios::app);

     string sFilename3="";
     sFilename3 = par::sOutDir + "/TGF/TGFB_" + string(itostr(par::iT)) + ".txt";
     fouttgfb.open(sFilename3,ios::app);

     double drug=0;  //BTZ: chemotherapy drug   0,0.1,0.2,...1   (maximal: 5nM)
     double codrug=0; //LEN: immunotherapy drug 0,0.1,0.2,...1   (maximal: ?)
	 double trddrug=0; //new drug to block the concentration of SDF1

	 sw_mic=par::sw_mic;
	 sw_mm=par::sw_mm;
	 sw_cd8=par::sw_cd8;  //swith to control the immune system
	 sw_treg=par::sw_treg;

     MatrixGrid = new double **[GridLength+1];

     WNT2= new double **[GridLength+1];
     DKK2= new double **[GridLength+1];

  //   double podrug=1.0e-1;
  //   double powdrug=pow(10,0.2);
  //   double codrugfac=2.0e-3;
  //// cout << "i" << "\t" << "codrug " << "\t" << "drug" << "\n";
  //   for ( int i=0;i<11;i++)
  //   {
  //      drug[i]=podrug*pow(powdrug,i);// drug: from 0.1 to 10, multiply e+0.2 every step
  //      codrug[i]=codrugfac*drug[i]; // codrug: from 2.0e-4 to 2.0e-2, multiply e+0.2 every step
  //    // cout << i << "\t" << codrug[i] << "\t" << drug[i] << "\n";
  //   }

	 for (int i=0;i<=GridLength;i++)
	 {
	     TGFb[i] = new double *[GridWidth+1];
		//WNT2[i] = new double *[GridWidth+1];
		//DKK2[i] = new double *[GridWidth+1];
		 SDF1[i] = new double *[GridWidth+1];
		 IL6[i] = new double *[GridWidth+1];
		 
		 /*drug[i] = new double *[GridWidth];
		 A[i] = new double *[GridWidth];*/
		 
		 STIFF[i] = new double *[GridWidth+1]; 
		 ECMHigh[i] = new double *[GridWidth+1];
		 MatrixGrid[i] = new double *[GridWidth+1];
		 ECMCan[i] = new double *[GridWidth+1];
		 ECM3[i] = new double *[GridWidth+1];
		 tmp[i] = new double *[GridWidth+1];
	     for (int j=0;j<=GridWidth;j++)
	     {
		    	TGFb[i][j] = new double [GridHeight+1];
			   //WNT2[i][j] = new double [GridHeight+1];
			   //DKK2[i][j] = new double [GridHeight+1];
			    SDF1[i][j] = new double [GridHeight+1];
			    IL6[i][j] = new double [GridHeight+1];
			    STIFF[i][j] = new double [GridHeight+1]; 
			    ECMHigh[i][j] = new double [GridHeight+1];
			    ECMCan[i][j] = new double [GridHeight+1];
				ECM3[i][j] = new double [GridHeight+1];
			    MatrixGrid[i][j] = new double [GridHeight+1];
			    //A[i][j] = new double [GridHeight+1];
			    tmp[i][j] = new double [GridHeight+1];
	     }
	}	   
}


void MicroEnv::MicroEnvInit()  //initialize the BM
{	
	int i,j,k,l;
	i1=0;
	for(i=0;i<=GridWidth;i++)
	{
	  for(j=0;j<=GridHeight;j++)
	  {
		  for (k=0;k<=GridLength;k++)
		  {
			  if (BoundaryCheck(i,j,k)==0)
				  {
					  // ========>  The initial values should be determined/given   <=========
					  TGFb[i][j][k]=3.0*pow(10.0,-7);
					  SDF1[i][j][k]=par::iniSDF1;   //default SDF1 concentration
					  tmp[i][j][k]=0;
					  IL6[i][j][k]=3.0*pow(10.0,-7);
					  STIFF[i][j][k]=par::inistiffness;  
					  ECMHigh[i][j][k]=NOCELL;
					  ECMCan[i][j][k]=NOCELL;
					  ECM3[i][j][k]=NOCELL;
					  MatrixGrid[i][j][k]=0;
					  
				   }
			  else
			  {
				  TGFb[i][j][k]=0;
				  SDF1[i][j][k]=par::basicSDF1; 
				  IL6[i][j][k]=0;
				  tmp[i][j][k]=0;
				  STIFF[i][j][k]=par::basicstiffness; // on boundary, the values are low
				  ECMHigh[i][j][k]=NOCELL;
				  ECMCan[i][j][k]=NOCELL;	
				  ECM3[i][j][k]=NOCELL;
				  MatrixGrid[i][j][k]=0;
				  //A[i][j][k]=0;//2.0e-3: co-drug; 0: no co-drug
			  }
		  }
	  }
	}

   
   for(i=0;i<=GridWidth;i+=4)
	{
	  for(j=0;j<=GridHeight;j+=4)
	  {
		  for (k=0;k<=GridLength;k+=4)
		  {
			  for (l=i;l<=i+4;l++)
			  {
			  if (BoundaryCheck(l,j+4,k)==0)
				  if(BoundaryCheck(l,j,k)==0)
					  if(BoundaryCheck(l,j,k+4)==0)
						  if(BoundaryCheck(l,j+4,k+4)==0)
						   {
							   if(ECMHigh[l][j+4][k]==NOCELL)
							   {
								   BMSC b1(i1++,l,j+4,k);	
								   BMSCList.insert(BMSCList.end(),b1);
							       ECMHigh[l][j+4][k]=BMSCCELL;
								   ECMCan[l][j+4][k]=BMSCCELL;
								   STIFF[l][j+4][k]=par::inistiffness;
							   }
							   if(ECMHigh[l][j][k]==NOCELL)
							   {
								   BMSC b2(i1++,l,j,k);	
								   BMSCList.insert(BMSCList.end(),b2);
								   ECMHigh[l][j][k]=BMSCCELL;
								   ECMCan[l][j][k]=BMSCCELL;
								   STIFF[l][j][k]=par::inistiffness;
							   }
							   if(ECMHigh[l][j][k+4]==NOCELL)
							   {
								   BMSC b3(i1++,l,j,k+4);	
								   BMSCList.insert(BMSCList.end(),b3);
								   ECMHigh[l][j][k+4]=BMSCCELL;
								   ECMCan[l][j][k+4]=BMSCCELL;
								   STIFF[l][j][k+4]=par::inistiffness;
							   } 
							   if(ECMHigh[l][j+4][k+4]==NOCELL)	
							   {
								   BMSC b4(i1++,l,j+4,k+4);	
								   BMSCList.insert(BMSCList.end(),b4);
								   ECMHigh[l][j+4][k+4]=BMSCCELL;
								   ECMCan[l][j+4][k+4]=BMSCCELL;
								   STIFF[l][j+4][k+4]=par::inistiffness;
							   }
						   }
			     }
			  }
		  }
   }	


   for(i=0;i<=GridWidth;i+=4)
	{
	  for(j=0;j<=GridHeight;j+=4)
	  {
		  for (k=0;k<=GridLength;k+=4)
		  {
			  for (l=j;l<=j+4;l++)
			  {
			  if (BoundaryCheck(i,l,k)==0)
				  if(BoundaryCheck(i+4,l,k)==0)
					  if(BoundaryCheck(i+4,l,k+4)==0)
						  if(BoundaryCheck(i,l,k+4)==0)
						   {
							   if(ECMHigh[i][l][k+4]==NOCELL)
							   {
								   BMSC b1(i1++,i,l,k+4);	
								   BMSCList.insert(BMSCList.end(),b1);
							 	   ECMHigh[i][l][k+4]=BMSCCELL;
								   ECMCan[i][l][k+4]=BMSCCELL;
								   STIFF[i][l][k+4]=par::inistiffness;
							   }
							   if(ECMHigh[i+4][l][k+4]==NOCELL)
							   {
								   BMSC b2(i1++,i+4,l,k+4);	
								   BMSCList.insert(BMSCList.end(),b2);
								   ECMHigh[i+4][l][k+4]=BMSCCELL;
								   ECMCan[i+4][l][k+4]=BMSCCELL;
								   STIFF[i+4][l][k+4]=par::inistiffness;
							   }
							   if(ECMHigh[i+4][l][k]==NOCELL)
							   {
								   BMSC b3(i1++,i+4,l,k);	
								   BMSCList.insert(BMSCList.end(),b3);
								   ECMHigh[i+4][l][k]=BMSCCELL;
								   ECMCan[i+4][l][k]=BMSCCELL;
								   STIFF[i+4][l][k]=par::inistiffness;
							   }
							   if(ECMHigh[i][l][k]==NOCELL)	
							   {
								   BMSC b4(i1++,i,l,k);	
								   BMSCList.insert(BMSCList.end(),b4);
								   ECMHigh[i][l][k]=BMSCCELL;
								   ECMCan[i][l][k]=BMSCCELL;
								   STIFF[i][l][k]=par::inistiffness;
							   }

						   }
                 }
			  }
		  }
   }



    for(i=0;i<=GridWidth;i+=4)
	{
	  for(j=0;j<=GridHeight;j+=4)
	  {
		  for (k=0;k<=GridLength;k+=4)
		  {
			  for (l=k;l<=k+4;l++)
			  {
			  if (BoundaryCheck(i,j+4,l)==0)
				  if(BoundaryCheck(i,j,l)==0)
					  if(BoundaryCheck(i+4,j,l)==0)
						  if(BoundaryCheck(i+4,j+4,l)==0)
						   {
							   if(ECMHigh[i][j+4][l]==NOCELL)
							   {
								   BMSC b1(i1++,i,j+4,l);	
								   BMSCList.insert(BMSCList.end(),b1);
								   ECMHigh[i][j+4][l]=BMSCCELL;
								   ECMCan[i][j+4][l]=BMSCCELL;
								   STIFF[i][j+4][l]=par::inistiffness;							   
								   
							   }
							   if(ECMHigh[i][j][l]==NOCELL)
							   {
								   BMSC b2(i1++,i,j,l);	
								   BMSCList.insert(BMSCList.end(),b2);
								   ECMHigh[i][j][l]=BMSCCELL;
								   ECMCan[i][j][l]=BMSCCELL;
								   STIFF[i][j][l]=par::inistiffness;
							   }
							   if(ECMHigh[i+4][j][l]==NOCELL)
							   {
								  BMSC b3(i1++,i+4,j,l);	
								  BMSCList.insert(BMSCList.end(),b3);
								  ECMHigh[i+4][j][l]=BMSCCELL;
								  ECMCan[i+4][j][l]=BMSCCELL;
								  STIFF[i+4][j][l]=par::inistiffness;
							   }
							   if(ECMHigh[i+4][j+4][l]==NOCELL)	
							   {
								  BMSC b4(i1++,i+4,j+4,l);	
								  BMSCList.insert(BMSCList.end(),b4);
								  ECMHigh[i+4][j+4][l]=BMSCCELL;
								  ECMCan[i+4][j+4][l]=BMSCCELL;
								  STIFF[i+4][j+4][l]=par::inistiffness;		 
							   }
						   }
                }
			  }
		  }
   }


	// INITIAL SEEDINGS ======>
	
	for(i=GridWidth/2-Radius;i<=GridWidth/2+Radius;i++)
	 {
	  for(j=GridHeight/2-Radius;j<=GridHeight/2+Radius;j++)
	   {
		   for(k=GridLength/2-Radius;k<=GridLength/2+Radius;k++)
		   { 
		     if(Distance(i,j,k,GridWidth/2,GridHeight/2,GridLength/2)<=Radius)  
				 if (ECMHigh[i][j][k]==NOCELL)			 
				 {
						double ran2=getRandom();
																				
						if ((ran2>0 && ran2<=0.5) && (sw_mic==1))
						   {	
							   if (MICLivingList.size()<100)
						       {
						    	MIC b(i1++,i,j,k);														
							    MICLivingList.insert(MICLivingList.end(),b);
							    ECMHigh[i][j][k]=CTYPEMIC;	
								ECM3[i][j][k]=CTYPEMIC;
							    SDF1[i][j][k]+=MICSDF1;										
						       }
						   }
						if ((ran2>0.5 && ran2<=1) && (sw_mm==1))
							{
						      if (MMLivingList.size()<100)
						      {							
							
								MM b(i1++,i,j,k);						
							    MMLivingList.insert(MMLivingList.end(),b);
							    ECMHigh[i][j][k]=CTYPEMM;
								ECM3[i][j][k]=CTYPEMM;
							    TGFb[i][j][k]+=MMTGFb;			
							  }
						}

				 }

		   }
	  }
	}

	for(i=GridWidth/2-Radius;i<=GridWidth/2+Radius;i++)
	 {
	  for(j=GridHeight/2-Radius;j<=GridHeight/2+Radius;j++)
	   {
		   for(k=GridLength/2-Radius;k<=GridLength/2+Radius;k++)
		   {
		     if(Distance(i,j,k,GridWidth/2,GridHeight/2,GridLength/2)<=Radius)  
				 if (ECMHigh[i][j][k]==NOCELL)			   
					{
						double ran3=getRandom();	
						if(ran3>0 && ran3<=0.5 && CTLLivingList.size() <20 && (sw_cd8==1))   //For CD8+
						{							 
							CTL b(i1++,i,j,k);						
							CTLLivingList.insert(CTLLivingList.end(),b);
							ECMHigh[i][j][k]=CTYPECTL;								
						}		
						
						double ran4=getRandom();
						if (ran3>0.5 && ran3<=1.0 && TregLivingList.size()<5 && (sw_treg==1))   //For CD8+ Treg
						{														
							Treg b(i1++,i,j,k);
							TregLivingList.insert(TregLivingList.end(),b);
							ECMHigh[i][j][k]=CTYPETREG;		
						}
					}
			}
	   }
	 }

	 //fout4.open("CellInfo.txt",ios::out);
	 //fout41.open("DeadCell.txt",ios::out);
	 //fout4<<"T MIC PC   MM   TMM---Living"<<endl;
	 //fout41<<"T MIC PC   MM   TMM---Dead"<<endl;
	 //foutsdf.open("SdfInfo.txt",ios::out);
	 //foutsdf<<"T	 SDF1total"<<endl;
	
}

void MicroEnv::Microdrug(int drugindex, int coindex, string isM)
{
    string sFilename = "";
    sFilename = par::sOutDir + "/" + isM + "_A_" + string(itostr(par::iA)) + "_B_" + string(itostr(par::iB)) + ".txt";
	//cout<<sFilename;
	//sFilename="c:\\file.txt";
	//fout4.open(sFilename,ios::out);
	//fout4<<"T\tMIC\tPC\tMM\tTMM"<<endl;
	
	//getchar();
	//FILE *fp_output2=fopen(sFilename,"w");
	//fprintf(fp_output2, "%d " ,3);
	//fclose(sFilename);
	//getchar();
}

void MicroEnv::MicroAddDrug(int Tstep, double drugindex, double coindex,double thirdidx)
{	
	//double maxsdf=0;
	//int a,b,c,coun;
	//double sum1=0;
	//coun=0;
	//for(int i=0;i<=GridWidth;i++)
	//{
	//  for(int j=0;j<=GridHeight;j++)
	//  {
	//	  for (int k=0;k<=GridLength;k++)
	//	  {
	//			 				 
	//		   //---------------------------------
	//			  if (BoundaryCheck(i,j,k)==0)
	//			  {
	//			   sum1=sum1+SDF1[i][j][k];
	//			   coun=coun+1;
	//			   if (SDF1[i][j][k]>maxsdf)
	//			     {
	//					 maxsdf=SDF1[i][j][k];	
	//					 a=i;
	//					 b=j;
	//					 c=k;
	//			     }	
	//			  }
	//	  }
	//  }
	//}
	//printf("drug SDF1:%f  [pos]:(%d,%d,%d)  %f  %f  %f  \n",maxsdf,a,b,c, ECMHigh[a][b][c],ECMCan[a][b][c],sum1/coun);
	//getchar();

            if ((Tstep>=100)&&(Tstep <= 172))  //the time to add drug // 6 days
		     { 
				   drug=drugindex;
				   codrug=coindex;
				   trddrug=thirdidx;				   
			  }
	        else
		     {
		        drug=0;
				codrug=0;
				trddrug=0;
	         }
   
}

double MicroEnv::diffusionGauss(int x, int y, int z)// IV, double mu_x, double mu_y, double mu_z, double theta are defined in micro
{
	double p,pdom,pnum1,pnum2,x1,x2,x3;
	x1=x-mux;
	x1=x1*x1;
	x2=y-muy;
	x2=x2*x2;
	x3=z-muz;
	x3=x3*x3;
	pnum1=-x1-x2-x3;

	pnum2=2*theta;
	pnum2=pnum2*theta;
	//pnum1=-(x-mux)*(x-mux)-(y-muy)*(y-muy)-(z-muz)*(z-muz);
	pdom=sqrt(2*PI);
	pdom=pdom*theta;
	p=IV;
	p=p/pdom;
	p=p*exp (  pnum1/pnum2  ) ; 
	//p= IV/(sqrt(2*PI)* theta) * exp (  ( -((double x)-mu_x)*(x-mu_x)-(y-mu_y)*(y-mu_y)-(z-mu_z)*(z-mu_z))/(2*theta*theta)  ) ; 
	return p;	
}

// THE RANDOM GENERATOR
double MicroEnv::getRandom()// get a die number between 0 and 1
{
	double b=rand()/(RAND_MAX+1.0);	
	return b;	
}

double MicroEnv::getMigD(int x, int y, int z) // get the RADIUS of MSearch function
{
	  double migD,migD1,kprol,migkd,migmax,d0;
	  migkd=migK_D;
	  migD1=(TGFb[x][y][z])/migkd;
	  migD1=pow(migD1,2);
	  kprol=K_prol;
	  kprol=pow(kprol,2);
	  migmax=migD_max;
	  d0=migD0;
	  migD=migD1/(1+migD1+kprol)*migmax+d0;
	  migD=migD*2/10; // convert to unit grid/timestep
	  return migD;
}

void MicroEnv::DiffTGFb(int Tstep, double lamda, double ***imtmp, double ***U)
{
	int i,j,k,l;
	double s=0,x=0,maxtgf=0;
	for(i=0;i<=GridWidth;i++)
	{
	  for(j=0;j<=GridHeight;j++)
	  {
		  for (k=0;k<=GridLength;k++)
		  {
			  if ((BoundaryCheck(i,j,k)==0))
			  {
     			  imtmp[i][j][k]=U[i][j][k];

    			  if (imtmp[i][j][k]>1)
	     			  imtmp[i][j][k]=1;
		    	  s=s+U[i][j][k];
			      x=x+1;
				   if (imtmp[i][j][k]>maxtgf)
				     maxtgf=imtmp[i][j][k];	
			  }
		  }
	  }
	}
	//-------------------------------------------------
	
	/*printf("\nDiffTGFB SDF1:%f  position:%d %d %d  %f  %f  %f %d\n",maxsdf,a,b,c, ECMHigh[a][b][c],ECMCan[a][b][c],sum1/coun,coun);
	getchar();*/
        fouttgfb<<Tstep<<" "<<maxtgf<<" "<<s/x<<endl;	
	fouttgfb<<"\n";
	//--------------------------------------------------

	for (l=1;l<=20;l++)
	{
	 for(i=0;i<=GridWidth;i++)
	  {
	    for(j=0;j<=GridHeight;j++)
	    {
		   for (k=0;k<=GridLength;k++)
		   {			 
			   if ((BoundaryCheck(i,j,k)==0))//
				 { 
					double localsum=0;
					int localnum=0;
					if (BoundaryCheck(i+1,j,k)==0)
					{
					 localsum=localsum+imtmp[i+1][j][k];
					 localnum=localnum+1;
					}
					if (BoundaryCheck(i-1,j,k)==0)
					{
					 localsum=localsum+imtmp[i-1][j][k];
					 localnum=localnum+1;
					}
					if (BoundaryCheck(i,j+1,k)==0)
					{
					 localsum=localsum+imtmp[i][j+1][k];
					 localnum=localnum+1;
					}
					if (BoundaryCheck(i,j-1,k)==0)
					{
					 localsum=localsum+imtmp[i][j-1][k];
					 localnum=localnum+1;
					}
					if (BoundaryCheck(i,j,k+1)==0)
					{
					 localsum=localsum+imtmp[i][j][k+1];
					 localnum=localnum+1;
					}
					if (BoundaryCheck(i,j,k-1)==0)
					{
					 localsum=localsum+imtmp[i][j][k-1];
					 localnum=localnum+1;
					}
					U[i][j][k]=((1-lamda)*imtmp[i][j][k]+lamda*(localsum/localnum))*(1-DEGT);
					//U[i][j][k]=((1-lamda)*imtmp[i][j][k]+lamda*(imtmp[i+1][j][k]+imtmp[i-1][j][k]+imtmp[i][j+1][k]+imtmp[i][j-1][k]+imtmp[i][j][k+1]+imtmp[i][j][k-1])/6)*(1-DEG);
		            if (U[i][j][k]<0.000001)
						U[i][j][k]=0.000001;
					if (U[i][j][k]>1)
						U[i][j][k]=1;
		         }
            }
	     }
	  }
	}
	printf(", AvgTGFb:%f\n",s/x);
}

void MicroEnv::Diff(int Tstep, double lamda, double ***imtmp, double ***U)  //for SDF1 diffusion
{	
	int i,j,k,l;

	double s=0;
	int x=0;
	double largesdf=0;
	double maxsdf=0;
	int a,b,c;
	for(i=0;i<=GridWidth;i++)
	{
	  for(j=0;j<=GridHeight;j++)
	  {
		  for (k=0;k<GridLength;k++)
		   {
				  imtmp[i][j][k]=U[i][j][k];

				  if (imtmp[i][j][k]>0.99)
					  imtmp[i][j][k]=0.99;

			   //---------------------------------
				  if (BoundaryCheck(i,j,k)==0)
				   {
				     s=s+U[i][j][k];
				     x=x+1;
				     if (imtmp[i][j][k]>0.1)
						largesdf=largesdf+1;	
				     if (imtmp[i][j][k]>maxsdf)
				     {
						 maxsdf=imtmp[i][j][k];	
						 a=i;
						 b=j;
						 c=k;
				     }	
				   }
		   }
	  }
	}

	//cout<<x;
    //getchar();
	printf("Dif[M_SDF]:%f [p]:%d,%d,%d %f, %f. %f",maxsdf,a,b,c, ECMHigh[a][b][c],ECMCan[a][b][c],s/x);	
	
	foutsdf1<<Tstep<<" "<<maxsdf<<" "<<s/x<<endl;	
	foutsdf1<<"\n";

	for (l=1;l<=20;l++)
	{
	 for(i=0;i<=GridWidth;i++)
	  {
	    for(j=0;j<=GridHeight;j++)
	    {
		   for (k=0;k<=GridLength;k++)
		   {			 
			   if ((BoundaryCheck(i,j,k)==0))//
				 { 
					
					double localsum=0;
					int localnum=0;
					if (BoundaryCheck(i+1,j,k)==0)
					{
					 localsum=localsum+imtmp[i+1][j][k];
					 localnum=localnum+1;					
					}
					if (BoundaryCheck(i-1,j,k)==0)
					{
					 localsum=localsum+imtmp[i-1][j][k];
					 localnum=localnum+1;				
					}
					if (BoundaryCheck(i,j+1,k)==0)
					{
					 localsum=localsum+imtmp[i][j+1][k];
					 localnum=localnum+1;				
					}
					if (BoundaryCheck(i,j-1,k)==0)
					{
					 localsum=localsum+imtmp[i][j-1][k];
					 localnum=localnum+1;					
					}
					if (BoundaryCheck(i,j,k+1)==0)
					{
					 localsum=localsum+imtmp[i][j][k+1];
					 localnum=localnum+1;					
					}
					if (BoundaryCheck(i,j,k-1)==0)
					{
					 localsum=localsum+imtmp[i][j][k-1];
					 localnum=localnum+1;					
					}
					
					U[i][j][k]=((1-lamda)*imtmp[i][j][k]+lamda*(localsum/localnum))*(1-DEG);
					//U[i][j][k]=((1-lamda)*imtmp[i][j][k]+lamda*(imtmp[i+1][j][k]+imtmp[i-1][j][k]+imtmp[i][j+1][k]+imtmp[i][j-1][k]+imtmp[i][j][k+1]+imtmp[i][j][k-1])/6)*(1-DEG);
		            if (U[i][j][k]<0.001)
						U[i][j][k]=0.001;
					if (U[i][j][k]>1)
						U[i][j][k]=1;



		         }
            }
	     }
	  }
	}

	//printf("\n");	
}	

void MicroEnv::DiffTrd(int Tstep, double lamda, double ***imtmp, double ***U)  //new drug for SDF1 inhibition
{	
	int i,j,k,l;
	double s=0;
	int x=0;
	double largesdf=0;
	double maxsdf=0;
	int a,b,c;
	for(i=0;i<=GridWidth;i++)
	{
	  for(j=0;j<=GridHeight;j++)
	  {
		  for (k=0;k<GridLength;k++)
		   {
				  imtmp[i][j][k]=U[i][j][k];  //get the value

				  if (imtmp[i][j][k]>0.99)
					  imtmp[i][j][k]=0.99;

			   //---------------------------------
				  if (BoundaryCheck(i,j,k)==0)
				   {
				     s=s+U[i][j][k];
				     x=x+1;
				     if (imtmp[i][j][k]>0.1)
						largesdf=largesdf+1;	
				     if (imtmp[i][j][k]>maxsdf)
				     {
						 maxsdf=imtmp[i][j][k];	
						 a=i;
						 b=j;
						 c=k;
				     }	
				   }
		   }
	  }
	}

	printf("degr[M_SDF1]:%f [p]%d,%d,%d %f, %f. %f ",maxsdf,a,b,c, ECMHigh[a][b][c],ECMCan[a][b][c],s/x);	
	//getchar();
	
	for (l=1;l<=20;l++)
	{
	 for(i=0;i<=GridWidth;i++)
	  {
	    for(j=0;j<=GridHeight;j++)
	    {
		   for (k=0;k<=GridLength;k++)
		   {			 
			   if ((BoundaryCheck(i,j,k)==0))//
				 { 			
					 double thmp,r,t;
                                         r=0.1212/(trddrug+0.1211);
                                         if (r>1)
                                             r=1;
                                         if (r<0)
                                             r=0;
                                         t=Thal_deg;//degration constant
                                         thmp=1-t*(1-r)/24;  //the effect will be occurred after 48hours					  
					
					 U[i][j][k]=imtmp[i][j][k]*thmp;  //drug-induced SDF1 degration
					
		             if (U[i][j][k]<0.001)
						  U[i][j][k]=0.001;
					 if (U[i][j][k]>1)
						  U[i][j][k]=1;
		         }
            }
	     }
	  }
	}		
}	

int MicroEnv::JudgeStatus(int CycleTimer)
{
    if(CYCLEG01BEGIN<=CycleTimer && CycleTimer<=CYCLEG01END)
		return StaG01;
	else
	{
		if (CYCLESBEGIN<=CycleTimer && CycleTimer<=CYCLESEND)
	     	  return StaS;
		else
		{
			if (CYCLEG2BEGIN<=CycleTimer && CycleTimer<=CYCLEG2END)
				return StaG2;
			else
			{
				if (CYCLEMBEGIN<=CycleTimer && CycleTimer<=CYCLEMEND)
					return StaM;
				else
					return StaN;
			}
		}
	}
}

double MicroEnv::Calculation_Stiff(int Tstep, double s)  //calculation of stiffness with developed lib
{
//       double e;	  
//	   //--initialize the lib of ODEs----------------------
//	   if (Tstep<=1)
//	   {	
//		if (!(libODEsInitialize()))
//			{
//			  cout<<"could not init lib!";
//			}
//	   }
//	 //  //-------------------------------------------------	   
// 	   
//	   if (s>=1)
//		  s=1;
//	   if (s<0.0012)
//		  s=0.0012;
//
//	   mwArray mwtype(1,1, mxDOUBLE_CLASS);
//	   mwArray mwsdf1(1,1, mxDOUBLE_CLASS);
//	   mwArray mwstiffness(1,1, mxDOUBLE_CLASS);
//						
//	   if (par::sM=="M")
//		    mwtype(1,1)=1;
//	   else
//		    mwtype(1,1)=0;						
//
//	   mwsdf1(1,1)=s;
//						
//	   cal_stiffness(1,mwstiffness,mwtype,mwsdf1);
//	   e=mwstiffness(1,1);		
//	   
//	   return e;
	return 1;
}

double MicroEnv::Calculation_Stiff2(int Tstep,double SDF1)
{
	  double par_m[NPAR]={0,0.5609,0.2244,0.2246,0.5004,0.5020,0.3456,0.9800,0.1000,0.5441,0.8905,0.1920,0.0910,0.3744,0.2199,0.7677,0.9960,0.3101,0.7731,0.9231,0.6381,0.6379,0.5342,0.1729,0.1232,0.9854,0.6788,0.6791};
	  double par_n[NPAR]={0,0.8594,0.1494,0.3855,0.4673,0.5308,0.3801,0.9072,0.9329,0.5899,0.1492,0.5239,0.5608,0.6446,0.4861,0.0924,0.0183,0.5503,0.2880,0.9817,0.7184,0.1092,0.6066,0.1791,0.8115,0.6934,0.3548,0.2909};
	  double *par;
	  int types;
	  double iniValue[NVAR]={1,1,1,1,1,1,1};
	  double iniTime = 0;
	  int ntimepoints = 2;
	  double timescale[2]={0,60};
	  double *sol;
	  int nsol = NVAR * ntimepoints;
	  int i,j;
	  double rst;
	 
	  if (SDF1>=1)
	      SDF1=1;
	  if (SDF1<0.0012)
		  SDF1=0.0012;
	  //----------------------------------------------------
	  types=1;  //parameters
	  if (types==0)
	  {
		  par=par_n;
	  }
	  else
	  {  
		  par=par_m; 	 
	  }	  
	  par[0]=SDF1; //given SDF1 concentration
	  //-------------------------------------------------------
	  
	  sol = (double*)malloc(nsol * sizeof(double));
	  if (!ODEmodel(par, iniValue, iniTime, ntimepoints, timescale, sol))
	    {
	      printf("Something wrong with ODEs solution!\n");
	      return 0;
	    }
	  ////-calculate stiffness---------------------------------
	  if (types==0)
		  rst=200*sol[2*NVAR-1];
	  else
		  rst=400*sol[2*NVAR-1];	  
	  
	  
	  if (sol) free(sol);
	  
	  return rst;
	   
}

void MicroEnv::BMSCCellCheck(int Tstep)
{
	double sumsti=0,ks=250,ss=0,sumsdf=0;
	int csti=0;
	double e=0;
	int a,b,c;
	list <BMSC>::iterator LivingListpo1;
	//---------------------------------------------------------------------------------------		   
	for (LivingListpo1=BMSCList.begin();LivingListpo1!=BMSCList.end();)
	{
		int TempTimer=LivingListpo1->BMSCCurrentTimer+1;
		if(TempTimer==Tstep)
		{
			LivingListpo1->BMSCCurrentTimer=LivingListpo1->BMSCCurrentTimer+1;
			//-------------------------------------------------------------------------------			
			
			//e=Calculation_Stiff(Tstep, SDF1[LivingListpo1->BMSCLx][LivingListpo1->BMSCLy][LivingListpo1->BMSCLz]);	
			e=Calculation_Stiff2(Tstep, SDF1[LivingListpo1->BMSCLx][LivingListpo1->BMSCLy][LivingListpo1->BMSCLz]);
                        //printf("%f  %f\n",SDF1[LivingListpo1->BMSCLx][LivingListpo1->BMSCLy][LivingListpo1->BMSCLz],e);
                       // getchar();
			STIFF[LivingListpo1->BMSCLx][LivingListpo1->BMSCLy][LivingListpo1->BMSCLz]=e;			
			//-------------------------------------------------------------------------------
			if (BoundaryCheck(LivingListpo1->BMSCLx,LivingListpo1->BMSCLy,LivingListpo1->BMSCLz)==0)
			{
				sumsti=sumsti+e;
			    csti=csti+1;
				sumsdf=sumsdf+SDF1[LivingListpo1->BMSCLx][LivingListpo1->BMSCLy][LivingListpo1->BMSCLz];
			    if (e>ks)
				    ks=e;
				if(SDF1[LivingListpo1->BMSCLx][LivingListpo1->BMSCLy][LivingListpo1->BMSCLz]>ss)
					{
						ss=SDF1[LivingListpo1->BMSCLx][LivingListpo1->BMSCLy][LivingListpo1->BMSCLz];
						a=LivingListpo1->BMSCLx;
						b=LivingListpo1->BMSCLy;
						c=LivingListpo1->BMSCLz;
				    } 
			}
			//-------------------------------------------------------------------------------			
			//BMSC also can secret TGFB
			double sti1,sti2,ksti,stif;
			sti1=e-par::basicstiffness;
			sti2=par::maxstiffness-par::basicstiffness;  //530-250pa
			ksti=0.5*sti2;
			stif=pow(sti1/ksti,2)/(1+pow(sti1/ksti,2));					
			TGFb[LivingListpo1->BMSCLx][LivingListpo1->BMSCLy][LivingListpo1->BMSCLz]=TGFb[LivingListpo1->BMSCLx][LivingListpo1->BMSCLy][LivingListpo1->BMSCLz]+stif*BMSCTGFb;
		}
		LivingListpo1++;		
	}	
	//cout<<csti;
	//getchar();
	printf("\nSDF1_BMSC %f  [POS]:(%d,%d,%d) %f %f\n",ss,a,b,c,ECMHigh[a][b][c],ECMCan[a][b][c]);	
	printf("AVG_SDF1:%f  AVG STIFF:%f, Max STIFF:%f\n",sumsdf/csti,sumsti/csti,ks);	
	

	foutstiff<<Tstep<<" "<<ks<<" "<<sumsti/csti<<endl;	
	foutstiff<<"\n";
}
  
double MicroEnv::Sur_MIC(int Tstep,double sti,double btz)
{
	 /* if (Tstep<=1)
	  {	
		if (!(libODEsInitialize()))
			{
			  cout<<"could not init lib!";
			}
	  }

	  double surdrug;
	  mwArray mwsti(1,1,mxDOUBLE_CLASS);
	  mwArray mwbtz(1,1,mxDOUBLE_CLASS);
	  mwArray mwsuv(1,1,mxDOUBLE_CLASS);
	  mwArray mwad(1,1,mxDOUBLE_CLASS);
	  	  	 
	  if(sti<300)
			sti=0.001+0.009*(sti-200)/100;
	  else if(sti<400)
			sti=0.01+0.09*(sti-300)/100;
	  else if (sti>400)
			sti=0.1+0.9*(sti-400)/130;
	  else
		    sti=1;
	    
	  mwbtz(1,1)=drug;  //dose of BTZ
	  mwsti(1,1)=sti;

	  Suv_Adh_MIC(2,mwsuv,mwad,mwsti,mwbtz);
	  surdrug=mwsuv(1,1);
	  return surdrug;*/
	  return 1;
}

double MicroEnv::Sur_MIC2(int Tstep,double sti,double btz)
{
	double par[N_PAR]={0,0.5932, 0.2164, 0.1284, 0.3480, 0.8460, 0.0310, 0.9893, 0.5240, 0.8107, 0.8363, 0.6183, 0.3256, 0.0447, 0.1347, 0.4187, 0.4953, 0.2110, 0.6662, 0.1271, 0.6820, 0.7570, 0.3366, 0.5238, 0.1182, 0.5872, 0.3867, 0.5586, 0.3591, 0.8945, 0.3076, 0.8035,0.7326,0.1023,0};
	
	if (sti<300)
	    sti=0.001+0.009*(sti-200)/100;
	else if(sti<400)
		sti=0.01+0.09*(sti-300)/100;
	else if(sti>400)
		sti=0.1+0.9*(sti-400)/130;
	
	par[0]=sti;
	par[34]=btz;
	  
    double iniValue[N_VAR]={1,1,1,1,1,1,1,0.1749,0.6897};  //ok
    double iniTime = 0;
    int ntimepoints = 3;
    double timescale[3]={0,48,192};
	double *sol;
	int nsol = N_VAR * ntimepoints;
	int i,j;
	double rst;
	//--------------------------------------------------------
	sol = (double*)malloc(nsol * sizeof(double));
	if (!ODEmodel2(par, iniValue, iniTime, ntimepoints, timescale, sol))
	{
	    printf("Something wrong with ODEs solution!\n");
	    return 0;
	 }
	   
	rst=sol[3*N_VAR-1];  //survival rate
	
	if ((btz==1) && (rst<0))
	       rst=0.022;  
	
	return rst;
	
}

void MicroEnv::MICCellCheck(int Tstep, int drugindex)
{		
  list <MIC>::iterator LivingListpo1;
  list <MIC>::iterator ApopListpo1;
  
  for (LivingListpo1=MICLivingList.begin();LivingListpo1!=MICLivingList.end();)
  {	  
	  int TempTimer=LivingListpo1->MICCurrentTimer+1;// every loop, the CurrentTimer is one less than the TimeStep 
	  //printf("NIC_ID,%d %d %d %f %f\n",LivingListpo1->MICLx,LivingListpo1->MICLy,LivingListpo1->MICLz,ECMHigh[LivingListpo1->MICLx][LivingListpo1->MICLy][LivingListpo1->MICLz],ECMCan[LivingListpo1->MICLx][LivingListpo1->MICLy][LivingListpo1->MICLz]);
     //---------------------------------------------------------------------------
	  if (TempTimer==Tstep) // to guarantee that every cell is handled only once in each timestep.
	  {
		  LivingListpo1->MICCurrentTimer=LivingListpo1->MICCurrentTimer+1;// Here just plus 1 inside the loop for considering about the new generated cells that added in each timestep
		  //printf("MIC___ID:%d\n",LivingListpo1->ID);
		  if((LivingListpo1->MICApopStatus==APOPY) || (ECM3[LivingListpo1->MICLx][LivingListpo1->MICLy][LivingListpo1->MICLz]==CTYPECTL))// if the cell is beginning to die
		  {   
			   if ((ECM3[LivingListpo1->MICLx][LivingListpo1->MICLy][LivingListpo1->MICLz]==CTYPECTL))
				{  
					ECM3[LivingListpo1->MICLx][LivingListpo1->MICLy][LivingListpo1->MICLz]=NOCELL;
			    }
			 
			  //printf("Kill___MIC_ID:%d\n",LivingListpo1->ID);
			  LivingListpo1->MICApopTimer=LivingListpo1->MICApopTimer+1;// after 5 timesteps, put it to MICAbsorbList			  
			  MICApopList.insert(MICApopList.end(),*LivingListpo1); // put it to being apoptosis cell list
			  LivingListpo1=MICLivingList.erase(LivingListpo1); // then delete it from the living cell list.
			 
		  }
		  else // the cell is not apoptosis
		  {
			  // To judge if apoptosis rate < threshold
			  double thresh=getRandom(); // generate threshold	
			  //----------------------------------------------------------------------------	 
	          double sti, ratiobtz,ratiocd8,surdrug,surcd8;	
	  
	          sti=STIFF[LivingListpo1->MICLx][LivingListpo1->MICLy][LivingListpo1->MICLz];
	          //ratiobtz=Sur_MIC(Tstep,sti,drug);  //calulate the survival rate based on our develop lib     
			  ratiobtz=Sur_MIC2(Tstep,sti,drug);

	         // printf("%f  %f  %f  %f", sti,ratiobtz, surdrug, ECMHigh[LivingListpo1->MICLx][LivingListpo1->MICLy][LivingListpo1->MICLz]);
                 //printf("Sur:%f",Sur_MIC2(100,530,0));
                 //getchar();	  	  
	  //---------------------------------------------------------------------------	  
	          double sizecd8,sizemic, b1, b2, rctlmic, apop0, apop1,apop;
	          sizecd8=CTLLivingList.size();
	          sizemic=MICLivingList.size()+MMLivingList.size();
	          rctlmic=sizecd8/sizemic;
	          ratiocd8=0.2998/(rctlmic+0.2354);
			  if (ratiocd8>1)
				  ratiocd8=1;
			  //-get the defined parameters--------
			  b1=myeloma_btz1;
			  if (par::sw_cd8==0)
                          {
                           b2=myeloma_off_cd8;
			  }
                          else
                          {
                           b2=myeloma_on_cd8;
                          }
                          apop0=MICapop0;
			  apop1=MICapopdrug;
			  //-----------------------------------
			  surdrug=pow(ratiobtz/b1,3)/(1+pow(ratiobtz/b1,3));
			  surcd8=pow(ratiocd8/b2,4)/(1+pow(ratiocd8/b2,4));   //CD8
	          			 
			  apop=apop0+apop1*(1-surdrug*surcd8);	 
	        // printf("MIC: ratio_btz:%f, ratio_cd8: %f, btz: %f, cd8: %f, apop:%f  \n",ratiobtz, ratiocd8, surdrug, surcd8, apop);	//very important          
                // getchar();	          
			  //-----------------------------------------------------------------------------
			  if (apop>thresh) // if apoptosis rate < Threshold  ===> apoptosis
			  {
				  LivingListpo1->MICApopTimer=LivingListpo1->MICApopTimer+1;// after 5 timesteps, put it to MICAbsorbList
				  //LivingListpo1->MICApopStatus=APOPY;
				  MICApopList.insert(MICApopList.end(),*LivingListpo1); // put it to being apoptosis cell list
				  LivingListpo1=MICLivingList.erase(LivingListpo1); // then delete it from the living cell list.
			  }
			  else // go on according to the flowchart
			  {				  
			   if (LivingListpo1->MICCellStatus!=StaN)// enter the cell cycle
			   {
				  LivingListpo1->MICCellStatus=JudgeStatus(LivingListpo1->MICCycleTimer);
				  switch (LivingListpo1->MICCellStatus)// Judge which stage the cell is in cell cycle
				  {
					  case StaG01: // in the G0/G1 stage
						  {
							  // after five time steps-> one time step is 2 hours
								 //* migration1 *//
							  Location *a;  
							  // THE SEARCH FUNCTION NEEDS TO BE CHANGED
							  //double getradius=getMigD(LivingListpo1->MICLx,LivingListpo1->MICLy,LivingListpo1->MICLz);
							  double micd=DMIC;	
							  a=MSearch(LivingListpo1->MICLx,LivingListpo1->MICLy,LivingListpo1->MICLz,micd);   //use the adhesion rate in this function
							  if (a!=null)// There are free spaces to migrate
							  {
			                                          if  (ECMCan[LivingListpo1->MICLx][LivingListpo1->MICLy][LivingListpo1->MICLz]==BMSCCELL)				
                                                               	       ECMHigh[LivingListpo1->MICLx][LivingListpo1->MICLy][LivingListpo1->MICLz]=BMSCCELL;		
							          else
                                                                       ECMHigh[LivingListpo1->MICLx][LivingListpo1->MICLy][LivingListpo1->MICLz]=NOCELL;  

                                                          	  ECM3[LivingListpo1->MICLx][LivingListpo1->MICLy][LivingListpo1->MICLz]=NOCELL;	
								  LivingListpo1->MICLx=a->PosX;
								  LivingListpo1->MICLy=a->PosY;
								  LivingListpo1->MICLz=a->PosZ;
								  ECMHigh[a->PosX][a->PosY][a->PosZ]=CTYPEMIC;	
								  ECM3[a->PosX][a->PosY][a->PosZ]=CTYPEMIC;
								  // MIC only secrete SDF1
								  SDF1[a->PosX][a->PosY][a->PosZ]=SDF1[a->PosX][a->PosY][a->PosZ]+MICSDF1;
								  LivingListpo1->MICCycleTimer=LivingListpo1->MICCycleTimer+1;// after one time step, should plus 1
								  LivingListpo1++;
							  } 
							  else // There are not free spaces for migration
							  {
								  LivingListpo1->MICCycleTimer=LivingListpo1->MICCycleTimer+1;
								  SDF1[LivingListpo1->MICLx][LivingListpo1->MICLy][LivingListpo1->MICLz]=SDF1[LivingListpo1->MICLx][LivingListpo1->MICLy][LivingListpo1->MICLz]+MICSDF1;
								  LivingListpo1++;
							  }

							 break;
 
						  }
					  case StaS: // in the S stage
						  {
							  //* migration1 *//
							  /* ===To set NEW LOCATION, migD is given*/
							  Location *a;
							  double micd=DMIC;
							  a=MSearch(LivingListpo1->MICLx,LivingListpo1->MICLy,LivingListpo1->MICLz,micd);  //use the adhesion rate in this function
							  if (a!=null)
							  {
								  if  (ECMCan[LivingListpo1->MICLx][LivingListpo1->MICLy][LivingListpo1->MICLz]==BMSCCELL)
                                                                       ECMHigh[LivingListpo1->MICLx][LivingListpo1->MICLy][LivingListpo1->MICLz]=BMSCCELL;	
                                                                  else
                                                                       ECMHigh[LivingListpo1->MICLx][LivingListpo1->MICLy][LivingListpo1->MICLz]=NOCELL;								 

                                                                  ECM3[LivingListpo1->MICLx][LivingListpo1->MICLy][LivingListpo1->MICLz]=NOCELL;
								  LivingListpo1->MICLx=a->PosX;
								  LivingListpo1->MICLy=a->PosY;
								  LivingListpo1->MICLz=a->PosZ;
								  ECMHigh[a->PosX][a->PosY][a->PosZ]=CTYPEMIC;	
								  ECM3[a->PosX][a->PosY][a->PosZ]=CTYPEMIC;
								  SDF1[a->PosX][a->PosY][a->PosZ]=SDF1[a->PosX][a->PosY][a->PosZ]+MICSDF1;
								  LivingListpo1->MICCycleTimer=LivingListpo1->MICCycleTimer+1;
								  LivingListpo1++;
							  }
							  else
							  {
								  // IF THERE ARE NO FREE SPACES, REMAIN IN THE SAME PLACE WITH SECRETING SDF1
								  LivingListpo1->MICCycleTimer=LivingListpo1->MICCycleTimer+1;
								  SDF1[LivingListpo1->MICLx][LivingListpo1->MICLy][LivingListpo1->MICLz]=SDF1[LivingListpo1->MICLx][LivingListpo1->MICLy][LivingListpo1->MICLz]+MICSDF1;
								  LivingListpo1++;
							  }
							  break;


						  }
					  case StaG2:  // in the G2 stage
						  {
							  //* migration1 *//
							  Location *a;
							  double micd=DMIC;
							  a=MSearch(LivingListpo1->MICLx,LivingListpo1->MICLy,LivingListpo1->MICLz,micd);  //use the adhesion rate in this function
							  if (a!=null)
							  {
								  if  (ECMCan[LivingListpo1->MICLx][LivingListpo1->MICLy][LivingListpo1->MICLz]==BMSCCELL)
                                                                       ECMHigh[LivingListpo1->MICLx][LivingListpo1->MICLy][LivingListpo1->MICLz]=BMSCCELL;	
								  else
                                                                       ECMHigh[LivingListpo1->MICLx][LivingListpo1->MICLy][LivingListpo1->MICLz]=NOCELL;                                                                      

                                                                  ECM3[LivingListpo1->MICLx][LivingListpo1->MICLy][LivingListpo1->MICLz]=NOCELL;
								  LivingListpo1->MICLx=a->PosX;
								  LivingListpo1->MICLy=a->PosY;
								  LivingListpo1->MICLz=a->PosZ;
								  ECMHigh[a->PosX][a->PosY][a->PosZ]=CTYPEMIC;
								  ECM3[a->PosX][a->PosY][a->PosZ]=CTYPEMIC;
								  SDF1[a->PosX][a->PosY][a->PosZ]=SDF1[a->PosX][a->PosY][a->PosZ]+MICSDF1;
								  LivingListpo1->MICCycleTimer=LivingListpo1->MICCycleTimer+1;
								  LivingListpo1++;
							  }
							  else
							  {
								  LivingListpo1->MICCycleTimer=LivingListpo1->MICCycleTimer+1;
								  SDF1[LivingListpo1->MICLx][LivingListpo1->MICLy][LivingListpo1->MICLz]=SDF1[LivingListpo1->MICLx][LivingListpo1->MICLy][LivingListpo1->MICLz]+MICSDF1;
								  LivingListpo1++;
			  				  }
							  break;

						  }
					  case StaM:
						  {
							  LivingListpo1->MICCycleTimer=CYCLEG01BEGIN; // first reset the Cell Cycle Timer to CYCLEG01BEGIN.
							 /* ===Search whether there are free spaces */
							  Location *a;
							  int r=radx;
							  
							  a=ProSearch(LivingListpo1->MICLx,LivingListpo1->MICLy,LivingListpo1->MICLz,r);  //find a neighbour position to proliferation
							  
							  if (a!=null)
								  LivingListpo1->MICFreeSpace=FREESPAY;
							  else
								  LivingListpo1->MICFreeSpace=FREESPAN;

							  if(LivingListpo1->MICFreeSpace==FREESPAY)// if there is free space, to proliferate (copy themselves)
							  {
								double dr = getRandom(); 
								double r,e0,emax,kre,r1,r2,r3,deltar;///p1,p2,p3;
								int count1=0;
								//------------------------------------------------------------
								r1=r_pathway;
								e0=par::basicstiffness;
								emax=par::maxstiffness;
								kre=(emax-e0)/2;
								r=0;
								for(int x=-1;x<=1;x++)
								{
								 for(int y=-1;y<=1;y++)
								  {
									for(int z=-1;z<=1;z++)
									{
										if( (BoundaryCheck(a->PosX+x,a->PosY+y,a->PosZ+z)==0)&& (ECMHigh[a->PosX+x][a->PosY+y][a->PosZ+z]==BMSCCELL) )	
									    {									 
										   r=r+STIFF[a->PosX+x][a->PosY+y][a->PosZ+z];
										   count1++;									
									    }
									}
								 }
							   }
							    if (count1!=0)
								    r=r/count1;															    
								r=(r-e0)/kre;
								r=pow(r,2);
								r=r/(1+r);								
								deltar=r1*r;
								r1=r_r0;								
								r1=r1+deltar;
								r2=r_a;
								r2=r1+r2;
								r3=r_d0;
								r3=r3-deltar;
								r3=r2+r3;				
								
								//---------------------------------------------------------------							
								if ( dr>=0 && dr<r1 )
									LivingListpo1->ProWhat=MICPROMIC;   //self renewal. 1MIC--2MIC
								if (dr>=r1 && dr<r2 )//r1==r2, no chance to go into
									LivingListpo1->ProWhat=MICPROPC;	// 1MIC---1MIC+1MM
								if (dr>=r2 && dr<r3 )
									LivingListpo1->ProWhat=MICCHANGEPC;  //1MIC---2MM
								
								//printf("%f, %d\n", dr, LivingListpo1->ProWhat);
								//getchar();
							 	/*if ((ECMCan[LivingListpo1->MICLx][LivingListpo1->MICLy][LivingListpo1->MICLz]==BMSCCELL)&&(a!=null))
							    {
								   printf("%d   %d   %d  %f  %d   %d   %d  %f\n",LivingListpo1->MICLx,LivingListpo1->MICLy,LivingListpo1->MICLz, ECMHigh[LivingListpo1->MICLx][LivingListpo1->MICLy][LivingListpo1->MICLz], a->PosX,a->PosY,a->PosZ,ECMHigh[a->PosX][a->PosY][a->PosZ]);
							       getchar();
							    }*/
								switch (LivingListpo1->ProWhat)
								{
								case MICPROMIC:  //1MIC--2MIC
									{
										//cout<<"1MIC--2MIC";
										//getchar();
										//printf("%d, %d, %d, %d\n",MICLivingList.size(),MICApopList.size(), MMLivingList.size(), MMApopList.size());

										MIC c(i1++,a->PosX,a->PosY,a->PosZ);
										c.MICCurrentTimer=Tstep;// to ensure that the new cell c will not be checked in this time step
										MICLivingList.insert(MICLivingList.end(),c);
										ECMHigh[a->PosX][a->PosY][a->PosZ]=CTYPEMIC;	
										ECM3[a->PosX][a->PosY][a->PosZ]=CTYPEMIC;
										SDF1[a->PosX][a->PosY][a->PosZ]+=MICSDF1;// new cell secrete SDF1										
										SDF1[LivingListpo1->MICLx][LivingListpo1->MICLy][LivingListpo1->MICLz]+=MICSDF1; // the SDF1 that original cell secrete
										LivingListpo1++;
										//getchar();
									    //printf("%d, %d, %d, %d\n",MICLivingList.size(),MICApopList.size(), MMLivingList.size(), MMApopList.size());
																				
										break;
									}
								case MICPROPC: //generate a new MM
									{
										//cout<<"generate a new MM";
										//getchar();
										//printf("%d, %d, %d, %d\n",MICLivingList.size(),MICApopList.size(), MMLivingList.size(), MMApopList.size());

										MM c(i1++,a->PosX,a->PosY,a->PosZ);
										c.MMCurrentTimer=Tstep;// to ensure that the new cell c will not be checked in this time step
										MMLivingList.insert(MMLivingList.end(),c);
										ECMHigh[a->PosX][a->PosY][a->PosZ]=CTYPEMM;	
										ECM3[a->PosX][a->PosY][a->PosZ]=CTYPEMM;
										IL6[a->PosX][a->PosY][a->PosZ]+=PCIL6;
										TGFb[a->PosX][a->PosY][a->PosZ]+=PCTGFb;																				
										SDF1[LivingListpo1->MICLx][LivingListpo1->MICLy][LivingListpo1->MICLz]+=MICSDF1; 
										LivingListpo1++;
										
										//getchar();
										//printf("%d, %d, %d, %d\n",MICLivingList.size(),MICApopList.size(), MMLivingList.size(), MMApopList.size());
										break;
									}
								case MICCHANGEPC:  //1MIC-->2MM
									{										
										//---------------------------------------------------------------------
										//cout<<"1MIC-->2MM";
										//getchar();
										//printf("%d, %d, %d, %d\n",MICLivingList.size(),MICApopList.size(), MMLivingList.size(), MMApopList.size());
										//---------------------------------------------------------------------
										MM c1(i1++,a->PosX,a->PosY,a->PosZ);
										c1.MMCurrentTimer=Tstep;// to ensure that the new cell c will not be checked in this time step
									    MMLivingList.insert(MMLivingList.end(),c1);													
									    ECMHigh[a->PosX][a->PosY][a->PosZ]=CTYPEMM;	
										ECM3[a->PosX][a->PosY][a->PosZ]=CTYPEMM;
										TGFb[a->PosX][a->PosY][a->PosZ]+=MMTGFb;
									    IL6[a->PosX][a->PosY][a->PosZ]+=MMIL6;									    									
									    
										//--------------------------------------------------------------------								
									    MM c2(LivingListpo1->ID,LivingListpo1->MICLx,LivingListpo1->MICLy,LivingListpo1->MICLz);
									    c2.MMCurrentTimer=Tstep;// to ensure that the cell changed from mother cell will not be checked in this time step
									    MMLivingList.insert(MMLivingList.end(),c2);										
									    ECMHigh[LivingListpo1->MICLx][LivingListpo1->MICLy][LivingListpo1->MICLz]=CTYPEMM;
										ECM3[LivingListpo1->MICLx][LivingListpo1->MICLy][LivingListpo1->MICLz]=CTYPEMM;
										TGFb[LivingListpo1->MICLx][LivingListpo1->MICLy][LivingListpo1->MICLz]+=MMTGFb;
									    IL6[LivingListpo1->MICLx][LivingListpo1->MICLy][LivingListpo1->MICLz]+=MMIL6;									   									
									    									
									    LivingListpo1=MICLivingList.erase(LivingListpo1);//delete MIC
										//---------------------------------------------------------------------								
										
										//getchar();
										//printf("%d, %d, %d, %d\n",MICLivingList.size(),MICApopList.size(), MMLivingList.size(), MMApopList.size());
										break;
									}
								}
								

							  }
							  else // there are no free spaces for proliferation ==> then stand by
							  {	
                                                                  printf("MM no room to generate\n");							  
								  LivingListpo1->MICCellStatus=StaM;
								  LivingListpo1->MICCycleTimer=CYCLEMBEGIN; // to ensure the cell is still in M Status for the next time step
								  SDF1[LivingListpo1->MICLx][LivingListpo1->MICLy][LivingListpo1->MICLz]+=MICSDF1;
								  LivingListpo1++;
							  }
						  }
				  }
			    }
			  else // Not cell cycle on
			  {
				  //Diff(D_TGFb,Kd_TGFb,theta_d,TGFb); // to get new values of TGFb at each grid point
				  double pathdo=par::basicstiffness;
				  double e0,emax,Ke,pathway,p0,pprol;
				  int count1=0;
				  e0=par::basicstiffness;				 
				  emax=par::maxstiffness;
				  Ke=(emax-e0)/2;
				  for(int x=-1;x<=1;x++)
					{
					 for(int y=-1;y<=1;y++)
					  {
						for(int z=-1;z<=1;z++)
						 {
							if( (BoundaryCheck(LivingListpo1->MICLx+x,LivingListpo1->MICLy+y,LivingListpo1->MICLz+z)==0)&& (ECMHigh[LivingListpo1->MICLx+x][LivingListpo1->MICLy+y][LivingListpo1->MICLz+z]==BMSCCELL) )	
						     {
							 pathdo+=STIFF[LivingListpo1->MICLx+x][LivingListpo1->MICLy+y][LivingListpo1->MICLz+z];
							 count1++;							 					
						     }
						 }
					  }
				    }
				  
				  if (count1!=0)
					  pathdo=pathdo/count1;	
				  
				  pathdo=(pathdo-e0)/Ke;
				  pathdo=pow(pathdo,2);
				  pathdo=pathdo/(1+pathdo);
				  p0=MICp0;
				  pathway=MICP_pathway;
				  pprol=p0+pathway*pathdo;  //proliferationrate for go into cell cycle
				 				 
				  double dr = getRandom();
				  
				  if (dr>=0 && dr<pprol) // enter the cell cycle
				  {
					  LivingListpo1->MICCycleTimer=CYCLEG01BEGIN;
					  LivingListpo1->MICCellStatus=StaG01;
					  LivingListpo1++;
				  }
				  else
				  {
					  Location *a;
					  double micd=DMIC;
					  a=MSearch(LivingListpo1->MICLx,LivingListpo1->MICLy,LivingListpo1->MICLz,micd);
					  if (a!=null)
					  {
						  if (ECMCan[LivingListpo1->MICLx][LivingListpo1->MICLy][LivingListpo1->MICLz]==BMSCCELL)
                                                      ECMHigh[LivingListpo1->MICLx][LivingListpo1->MICLy][LivingListpo1->MICLz]=BMSCCELL;		
						  else
                                                      ECMHigh[LivingListpo1->MICLx][LivingListpo1->MICLy][LivingListpo1->MICLz]=NOCELL;

                                                  ECM3[LivingListpo1->MICLx][LivingListpo1->MICLy][LivingListpo1->MICLz]=NOCELL;	
						  LivingListpo1->MICLx=a->PosX;
						  LivingListpo1->MICLy=a->PosY;
						  LivingListpo1->MICLz=a->PosZ;
						  ECMHigh[a->PosX][a->PosY][a->PosZ]=CTYPEMIC;		
						  ECM3[a->PosX][a->PosY][a->PosZ]=CTYPEMIC;
						  LivingListpo1->MICCycleTimer=LivingListpo1->MICCycleTimer+1;
						  SDF1[a->PosX][a->PosY][a->PosZ]+=MICSDF1;
						  LivingListpo1++;
						  //LivingListpo1->HHCurrentTimer=LivingListpo1->HHCurrentTimer+1;
					  }
					  else
					  {
						  LivingListpo1->MICCycleTimer=LivingListpo1->MICCycleTimer+1;
						  SDF1[LivingListpo1->MICLx][LivingListpo1->MICLy][LivingListpo1->MICLz]+=MICSDF1;
						  LivingListpo1++;
					  }

				  }
			  
			  }
		  }
		  }

	  }
	  else
		  LivingListpo1++;
  }
 
  for (ApopListpo1=MICApopList.begin();ApopListpo1!=MICApopList.end();)// to search Apoplist
  {
	  int Timer1=ApopListpo1->MICCurrentTimer+1;//every loop, the CurrentTimer is the same to the TimeStep
		  
	  if (Timer1==Tstep)// if the cell is inserted from the above process, then its MICCurrentTimer+1!=Tstep
		                                     // therefore it will not be considered in this process
	  {
		   ApopListpo1->MICApopTimer = ApopListpo1->MICApopTimer+1;
		   ApopListpo1->MICCurrentTimer=ApopListpo1->MICCurrentTimer+1;
		  if (ApopListpo1->MICApopTimer >=5)
		  {
			  ApopListpo1->MICApopStatus=APOPA;
			  ECMHigh[ApopListpo1->MICLx][ApopListpo1->MICLy][ApopListpo1->MICLz]=NOCELL;
			  ECM3[ApopListpo1->MICLx][ApopListpo1->MICLy][ApopListpo1->MICLz]=NOCELL;
			  MICAbsorbList.insert(MICAbsorbList.end(),*ApopListpo1);
			  ApopListpo1=MICApopList.erase(ApopListpo1);
		  }
		  else
			  ApopListpo1++;
	  }
	  else
		  ApopListpo1++;
  }

}

void MicroEnv::MMCellCheck(int Tstep, int drugindex)
{
		
  list <MM>::iterator LivingListpo1;
  list <MM>::iterator ApopListpo1;

  for (LivingListpo1=MMLivingList.begin();LivingListpo1!=MMLivingList.end();)
  {
	  int TempTimer=LivingListpo1->MMCurrentTimer+1;//every loop, the CurrentTimer is one less than the TimeStep	  
	  if (TempTimer==Tstep) // to guarantee that every cell is handled only once in every timestep.
	  {
		  //printf("MM___ID:%d\n",LivingListpo1->ID);
		  LivingListpo1->MMCurrentTimer=LivingListpo1->MMCurrentTimer+1;// Here just plus 1 inside the loop for considering about the new MM cells that added in this timestep
		  if((LivingListpo1->MMApopStatus==APOPY)||(ECM3[LivingListpo1->MMLx][LivingListpo1->MMLy][LivingListpo1->MMLz]==CTYPECTL))// if the cell is beginning to die
		  {
			  //printf("Kill___MM___ID:%d\n",LivingListpo1->ID);
			  if (ECM3[LivingListpo1->MMLx][LivingListpo1->MMLy][LivingListpo1->MMLz]==CTYPECTL)
			  {
				 ECM3[LivingListpo1->MMLx][LivingListpo1->MMLy][LivingListpo1->MMLz]=NOCELL;					
			  }
			
			  LivingListpo1->MMApopTimer=LivingListpo1->MMApopTimer+1;// after 5 timesteps, put it to MMAbsorbList
			  MMApopList.insert(MMApopList.end(),*LivingListpo1); // put it to being apoptosis cell list
			  LivingListpo1=MMLivingList.erase(LivingListpo1); // then delete it from the living cell list.
			  			  
		  }
		  else // the cell is not apoptosis
		  {
			  double thresh=getRandom(); // generate threshold
			  //--------------------------------------------------------------
			  double apop0,apop, apopdrug;
			  double rctlmm, pcsize,mmsize;
			  double surb, surc, drugbtz, ratio_mm_btz, ratio_mm_cd8,b3,b4;
			  drugbtz=drug*5;  //the dose of btz;
			  pcsize=CTLLivingList.size();
			  mmsize=MMLivingList.size()+MICLivingList.size();
			  rctlmm=pcsize/mmsize;			  
			  
			  //ratio_mm_btz=(423.4*drugbtz+1734)/(drugbtz*drugbtz*drugbtz-106.4*drugbtz*drugbtz+1313*drugbtz+1703);
			  ratio_mm_btz=(0.4632*drugbtz*drugbtz-0.3839*drugbtz+0.6247)/(drugbtz*drugbtz-0.4306*drugbtz+0.6246);
			  if (ratio_mm_btz>1)
				  ratio_mm_btz=1;
			  ratio_mm_cd8=0.2998/(rctlmm+0.2354);
			  if (ratio_mm_cd8 >1)
				  ratio_mm_cd8 = 1;
			  b3=myeloma_btz2;
			  if (par::sw_cd8==0)
                          {
                            b4=myeloma_off_cd8;
			  }
                          else
                          {
                            b4=myeloma_on_cd8;
                          }
                          apop0=MMapop0;
			  apopdrug=MMapopdrug;
			  surb= pow(ratio_mm_btz/b3,4)/(1+pow(ratio_mm_btz/b3,4));  //BTZ-induced apoptosis
			  surc= pow(ratio_mm_cd8/b4,4)/(1+pow(ratio_mm_cd8/b4,4));  //CD8-induced lysis
			  apop=apop0+apopdrug*(1-surb*surc);	

			  //printf("MM: ratio_btz:%f, ratio_cd8: %f, btz: %f, cd8: %f, apop: %f\n",ratio_mm_btz,ratio_mm_cd8,surb,surc,apop); //very important
			  
			  //---------------------------------------------------------------------------
			  if (apop>thresh) // if apoptosis rate < Threshold  ===> apoptosis
			  {
				  LivingListpo1->MMApopTimer=LivingListpo1->MMApopTimer+1;// after 5 timesteps, put it to MICAbsorbList
				  //LivingListpo1->MICApopStatus=APOPY;
				  MMApopList.insert(MMApopList.end(),*LivingListpo1); // put it to being apoptosis cell list
				  LivingListpo1=MMLivingList.erase(LivingListpo1); // then delete it from the living cell list.
			  }
			  else
			  {
			  
			  if (LivingListpo1->MMCellStatus!=StaN)// enter the cell cycle
			  {
				  LivingListpo1->MMCellStatus=JudgeStatus(LivingListpo1->MMCycleTimer);
				  switch (LivingListpo1->MMCellStatus)// Judge which stage the cell is in cell cycle
				  {
					  case StaG01: // in the G0/G1 stage
						  {
							  // after five time steps-> one time step is 2 hours
								 //* migration1 *//
							  Location *a;  
							  double micd=DMM;
							  a=MSearch_mm(LivingListpo1->MMLx,LivingListpo1->MMLy,LivingListpo1->MMLz,micd);
							 
							  if (a!=null)// There are free space to migrate
							  {
					
                                                                  if (ECMCan[LivingListpo1->MMLx][LivingListpo1->MMLy][LivingListpo1->MMLz]==BMSCCELL)
                                         			      ECMHigh[LivingListpo1->MMLx][LivingListpo1->MMLy][LivingListpo1->MMLz]=BMSCCELL;
					                          else
                                                                      ECMHigh[LivingListpo1->MMLx][LivingListpo1->MMLy][LivingListpo1->MMLz]=NOCELL;
                                                                      
                                         			  ECM3[LivingListpo1->MMLx][LivingListpo1->MMLy][LivingListpo1->MMLz]=NOCELL;
								  LivingListpo1->MMLx=a->PosX;
								  LivingListpo1->MMLy=a->PosY;
								  LivingListpo1->MMLz=a->PosZ;
								  ECMHigh[a->PosX][a->PosY][a->PosZ]=CTYPEMM;	
								  ECM3[a->PosX][a->PosY][a->PosZ]=CTYPEMM;
								  TGFb[a->PosX][a->PosY][a->PosZ]=TGFb[a->PosX][a->PosY][a->PosZ]+MMTGFb;
								  LivingListpo1->MMCycleTimer=LivingListpo1->MMCycleTimer+1;// after one time step, should plus 1
								  LivingListpo1++;
							  }
							  else // There are no free spaces for migration, then stay in the same place
							  {
								  LivingListpo1->MMCycleTimer=LivingListpo1->MMCycleTimer+1;
								  TGFb[LivingListpo1->MMLx][LivingListpo1->MMLy][LivingListpo1->MMLz]+=MMTGFb;		
								  LivingListpo1++;
							  }

							 break;
 
						  }
					  case StaS: // in the S stage
						  {
							  //* migration1 *//
							  /* ===To set NEW LOCATION, migD is given*/
							  Location *a;
							  double micd=DMM;
							  a=MSearch_mm(LivingListpo1->MMLx,LivingListpo1->MMLy,LivingListpo1->MMLz,micd);
							  
							  if (a!=null)
							  {
								  if  (ECMCan[LivingListpo1->MMLx][LivingListpo1->MMLy][LivingListpo1->MMLz]==BMSCCELL)
                                                                       ECMHigh[LivingListpo1->MMLx][LivingListpo1->MMLy][LivingListpo1->MMLz]=BMSCCELL;	
                                                                  else								 
                                                                       ECMHigh[LivingListpo1->MMLx][LivingListpo1->MMLy][LivingListpo1->MMLz]=NOCELL;
                                                                    
                                                                  ECM3[LivingListpo1->MMLx][LivingListpo1->MMLy][LivingListpo1->MMLz]=NOCELL;
								  LivingListpo1->MMLx=a->PosX;
								  LivingListpo1->MMLy=a->PosY;
								  LivingListpo1->MMLz=a->PosZ;
								  ECMHigh[a->PosX][a->PosY][a->PosZ]=CTYPEMM;	
								  ECM3[a->PosX][a->PosY][a->PosZ]=CTYPEMM;
								  TGFb[a->PosX][a->PosY][a->PosZ]=TGFb[a->PosX][a->PosY][a->PosZ]+MMTGFb;									  
								  LivingListpo1->MMCycleTimer=LivingListpo1->MMCycleTimer+1;
								  LivingListpo1++;
							  }
							  else
							  {
								  LivingListpo1->MMCycleTimer=LivingListpo1->MMCycleTimer+1;
								  TGFb[LivingListpo1->MMLx][LivingListpo1->MMLy][LivingListpo1->MMLz]=TGFb[LivingListpo1->MMLx][LivingListpo1->MMLy][LivingListpo1->MMLz]+MMTGFb;
								  LivingListpo1++;
							  }
							  break;
						  }
					  case StaG2:  // in the G2 stage
						  {
							  //* migration1 *//
							  Location *a;
							  double micd=DMM;
							  a=MSearch_mm(LivingListpo1->MMLx,LivingListpo1->MMLy,LivingListpo1->MMLz,micd);
							 
							  if (a!=null)
							  {
								  if (ECMCan[LivingListpo1->MMLx][LivingListpo1->MMLy][LivingListpo1->MMLz]==BMSCCELL)
                                                                      ECMHigh[LivingListpo1->MMLx][LivingListpo1->MMLy][LivingListpo1->MMLz]=BMSCCELL;
								  else
                                                                      ECMHigh[LivingListpo1->MMLx][LivingListpo1->MMLy][LivingListpo1->MMLz]=NOCELL;                            
   
                                                                  ECM3[LivingListpo1->MMLx][LivingListpo1->MMLy][LivingListpo1->MMLz]=NOCELL;
								  LivingListpo1->MMLx=a->PosX;
								  LivingListpo1->MMLy=a->PosY;
								  LivingListpo1->MMLz=a->PosZ;
								  ECMHigh[a->PosX][a->PosY][a->PosZ]=CTYPEMM;
								  ECM3[a->PosX][a->PosY][a->PosZ]=CTYPEMM;								  
								  TGFb[a->PosX][a->PosY][a->PosZ]=TGFb[a->PosX][a->PosY][a->PosZ]+MMTGFb;								
								  LivingListpo1->MMCycleTimer=LivingListpo1->MMCycleTimer+1;
								  LivingListpo1++;
							  }
							  else
							  {
								  LivingListpo1->MMCycleTimer=LivingListpo1->MMCycleTimer+1;
								  TGFb[LivingListpo1->MMLx][LivingListpo1->MMLy][LivingListpo1->MMLz]=TGFb[LivingListpo1->MMLx][LivingListpo1->MMLy][LivingListpo1->MMLz]+MMTGFb;
								  LivingListpo1++;
							  }
							  break;

						  }
					  case StaM:
						  {
							  LivingListpo1->MMCycleTimer=CYCLEG01BEGIN; // first reset the Cell Cycle Timer to CYCLEG01BEGIN.
							  //HHMStatus=StaG01; 
								/* ===Search whether there are free spaces */
							  Location *a;
							  int r=radx;
							  a=ProSearch(LivingListpo1->MMLx,LivingListpo1->MMLy,LivingListpo1->MMLz,r);
							  
							  if (a!=null)
								  LivingListpo1->MMFreeSpace=FREESPAY;
							  else
								  LivingListpo1->MMFreeSpace=FREESPAN;

							  if(LivingListpo1->MMFreeSpace==FREESPAY)
							  {
								  if (LivingListpo1->MMGenNo>MAXMM)
								  {							
								        MMApopList.insert(MMApopList.end(),*LivingListpo1);
                                                                 	LivingListpo1=MMLivingList.erase(LivingListpo1);  //delete MM
								  }
								  else
								  {
								    MM c(i1++,a->PosX,a->PosY,a->PosZ);
									c.MMGenNo=LivingListpo1->MMGenNo+1;
									c.MMCurrentTimer=Tstep;// to ensure that the new cell c will not be checked in this time step
									ECMHigh[a->PosX][a->PosY][a->PosZ]=CTYPEMM;	
									ECM3[a->PosX][a->PosY][a->PosZ]=CTYPEMM;
									TGFb[a->PosX][a->PosY][a->PosZ]+=MMTGFb;																	
									MMLivingList.insert(MMLivingList.end(),c);

									TGFb[LivingListpo1->MMLx][LivingListpo1->MMLy][LivingListpo1->MMLz]+=MMTGFb;								
									LivingListpo1->MMGenNo=LivingListpo1->MMGenNo+1;
																		
									LivingListpo1++;

									  
								  }
							  }
							  else // there are no free spaces for pro, stand by
							  {
								 LivingListpo1->MMCellStatus=StaM;
								 LivingListpo1->MMCycleTimer=CYCLEMBEGIN; // to ensure the cell is still in M Status for the next time step
								 LivingListpo1++;
							  }

						  }
				  }
			  }
			  else // Not cell cycle on
			  {				  			  	

				  double pathdo=0;
				  double e0,emax,Ke,pathway,p0,pprol;
				  int count1=0;
				  e0=E0;
				  emax=EMAX;
				  Ke=(emax-e0)/2;

				  for(int x=-1;x<=1;x++)
				   {
					 for(int y=-1;y<=1;y++)
					  {
						for(int z=-1;z<=1;z++)
						{
						  if( (BoundaryCheck(LivingListpo1->MMLx+x,LivingListpo1->MMLy+y,LivingListpo1->MMLz+z)==0)&& (ECMHigh[LivingListpo1->MMLx+x][LivingListpo1->MMLy+y][LivingListpo1->MMLz+z]==BMSCCELL) )	
						  {						 
							 pathdo+=STIFF[LivingListpo1->MMLx+x][LivingListpo1->MMLy+y][LivingListpo1->MMLz+z];
							 count1++;						 
						  }
						}
					 }
				   }
	
				  if (count1!=0)
					  pathdo=pathdo/count1;
				  pathdo=(pathdo-e0)/Ke;
				  pathdo=pow(pathdo,2);
				  pathdo=pathdo/(1+pathdo);
				  p0=MMp0;
				  pathway=MMP_pathway;
				  pprol=p0+pathway*pathdo;				
				  
				  double dr =getRandom();
				  
				  if (dr>=0 && dr<pprol)
				  {
					  LivingListpo1->MMCycleTimer=CYCLEG01BEGIN;
					  LivingListpo1->MMCellStatus=StaG01;
					  LivingListpo1++;
				  }
				  else
				  {
					  Location *a; // SEARCH FUNCTION SHOULD ALSO BE CHANGED
					  double micd=DMM;
					  a=MSearch_mm(LivingListpo1->MMLx,LivingListpo1->MMLy,LivingListpo1->MMLz,micd);  //should be checked
					  if (a!=null)
					  {
						  if (ECMCan[LivingListpo1->MMLx][LivingListpo1->MMLy][LivingListpo1->MMLz]==BMSCCELL)
                                                       ECMHigh[LivingListpo1->MMLx][LivingListpo1->MMLy][LivingListpo1->MMLz]=BMSCCELL;	
						  else
                                                       ECMHigh[LivingListpo1->MMLx][LivingListpo1->MMLy][LivingListpo1->MMLz]=NOCELL;                    

                                                  ECM3[LivingListpo1->MMLx][LivingListpo1->MMLy][LivingListpo1->MMLz]=NOCELL;
						  LivingListpo1->MMLx=a->PosX;
						  LivingListpo1->MMLy=a->PosY;
						  LivingListpo1->MMLz=a->PosZ;
						  ECMHigh[a->PosX][a->PosY][a->PosZ]=CTYPEMM;		
						  ECM3[a->PosX][a->PosY][a->PosZ]=CTYPEMM;
						  LivingListpo1->MMCycleTimer=LivingListpo1->MMCycleTimer+1;
						  TGFb[a->PosX][a->PosY][a->PosZ]+=MMTGFb;
						  //IL6[a->PosX][a->PosY][a->PosZ]+=MMIL6;
						  
						  LivingListpo1++;
						  //LivingListpo1->HHCurrentTimer=LivingListpo1->HHCurrentTimer+1;
					  }
					  else
					  {
						  LivingListpo1->MMCycleTimer=LivingListpo1->MMCycleTimer+1;
						  TGFb[LivingListpo1->MMLx][LivingListpo1->MMLy][LivingListpo1->MMLz]+=MMTGFb;
						  //IL6[LivingListpo1->MMLx][LivingListpo1->MMLy][LivingListpo1->MMLz]+=MMIL6;						 
						  LivingListpo1++;
					  }

				  }

			  }
			  }
		  }
	  }
	  else
		  LivingListpo1++;
  }
  //getchar();
  for (ApopListpo1=MMApopList.begin();ApopListpo1!=MMApopList.end();)// to search Apoplist
  {
	  int Timer1=ApopListpo1->MMCurrentTimer+1;//every loop, the CurrentTimer is the same to the TimeStep
	  if (Timer1==Tstep)// if the cell is inserted from the above process, then its MMCurrentTimer+1!=Tstep
		                                     // therefore it will not be considered in this process
	  {
		   ApopListpo1->MMApopTimer = ApopListpo1->MMApopTimer+1;
		   ApopListpo1->MMCurrentTimer=ApopListpo1->MMCurrentTimer+1;
		  if (ApopListpo1->MMApopTimer >=5)
		  {
			  ApopListpo1->MMApopStatus=APOPA;
			  ECMHigh[ApopListpo1->MMLx][ApopListpo1->MMLy][ApopListpo1->MMLz]=NOCELL;
			  ECM3[ApopListpo1->MMLx][ApopListpo1->MMLy][ApopListpo1->MMLz]=NOCELL;
			  MMAbsorbList.insert(MMAbsorbList.end(),*ApopListpo1);
			  ApopListpo1=MMApopList.erase(ApopListpo1);
		  }
		  else
			  ApopListpo1++;
	  }
	  else
		  ApopListpo1++;
  }

}

void MicroEnv::CTLCellCheck(int Tstep, int drugindex)
{
  list <CTL>::iterator LivingListpo1;
  list <CTL>::iterator ApopListpo1;
  
  for (LivingListpo1=CTLLivingList.begin();LivingListpo1!=CTLLivingList.end();)
  {
	  int TempTimer=LivingListpo1->CTLCurrentTimer+1;//every loop, the CurrentTimer is one less than the TimeStep	  
	  
	  if (TempTimer==Tstep) // to guarantee that every cell is dealt with only once in every timestep.
	  {
		  LivingListpo1->CTLCurrentTimer=LivingListpo1->CTLCurrentTimer+1;// Here just plus 1 inside the loop for considering about the new PC cells that added in this timestep
		  //printf("ID:%d, Cycletimer:%d, Apopst:%d\n",LivingListpo1->ID,LivingListpo1->CTLCycleTimer,LivingListpo1->CTLApopStatus);
		  if ((ECM3[LivingListpo1->CTLLx][LivingListpo1->CTLLy][LivingListpo1->CTLLz]==CTYPETREG) && (LivingListpo1->CTLCycleTimer>=1))  //cell cycle arrest or apoptosis
			  {			  
				  double th=getRandom();
				  if (th<0.25)//randomly choose: cell cycle arrest or apoptosis induced by Treg
			  	  {
					  LivingListpo1->CTLCycleTimer=CYCLEG01BEGIN;		
				  }
				  else
				  {				   
					  LivingListpo1->CTLApopStatus=APOPY;				    
				  }
				  ECM3[LivingListpo1->CTLLx][LivingListpo1->CTLLy][LivingListpo1->CTLLz]=NOCELL;
				  //printf("[%d, %d, %d]Inhibited by TREG: ID:%d, Cycletimer:%d, Apopst:%d",LivingListpo1->CTLLx,LivingListpo1->CTLLy,LivingListpo1->CTLLz,LivingListpo1->ID,LivingListpo1->CTLCycleTimer,LivingListpo1->CTLApopStatus);	
				  //getchar();
			  }
		  //-------------------------------------------------------------------------------------------------------------		 
		  if(LivingListpo1->CTLApopStatus==APOPY)// if the cell is beginning to die
		  {
			  LivingListpo1->CTLApopTimer=LivingListpo1->CTLApopTimer+1;// after 5 timesteps, put it to HHAbsorbList
			  CTLApopList.insert(CTLApopList.end(),*LivingListpo1); // put it to being apoptosis cell list
			  LivingListpo1=CTLLivingList.erase(LivingListpo1); // then delete it from the living cell list.
		  }
		  else // the cell is not apoptosis
		  {			  
			  // To judge if apoptosis rate < threshold
			  double thresh=getRandom(); // generate threshold
			  double cur_tgfbval,tmp,apop;
			  cur_tgfbval=TGFb[LivingListpo1->CTLLx][LivingListpo1->CTLLy][LivingListpo1->CTLLz];
			  tmp=pow(cur_tgfbval/CTLTGFB0,1);			 
			  apop=CTLapop0+CTLTGFBapop*(1+tmp/(1+tmp));	//

			  			  
			  if (apop>thresh) // if apoptosis rate < Threshold  ===> apoptosis
			  {
				  LivingListpo1->CTLApopTimer=LivingListpo1->CTLApopTimer+1;// after 5 timesteps, put it to MICAbsorbList
				  CTLApopList.insert(CTLApopList.end(),*LivingListpo1); // put it to being apoptosis cell list
				  LivingListpo1=CTLLivingList.erase(LivingListpo1); // then delete it from the living cell list.
			  }
			  else // apoptosis rate > Threshold, go on following the flowchart
			  {
			  	 
			  if (LivingListpo1->CTLCellStatus!=StaN)// enter the cell cycle
			  {
				  LivingListpo1->CTLCellStatus=JudgeStatus(LivingListpo1->CTLCycleTimer);
				  //printf("CTL, ID: %d  apopotosis rate:%f   %d\n",LivingListpo1->ID, apop,LivingListpo1->CTLCellStatus);
				  switch (LivingListpo1->CTLCellStatus)// Judge which stage the cell is in cell cycle
				  {
					  case StaG01: // in the G0/G1 stage
						  {
							  Location *a;  
							  double micd=DPC;
							  a=MSearch_ctl(LivingListpo1->CTLLx,LivingListpo1->CTLLy,LivingListpo1->CTLLz,micd);  //will change
							  if (a!=null)// There are free space to migrate
							  {
								  if ((ECMHigh[a->PosX][a->PosY][a->PosZ]==CTYPEMM) || (ECMHigh[a->PosX][a->PosY][a->PosZ]==CTYPEMIC))
								  {									  
									  ECM3[a->PosX][a->PosY][a->PosZ]=CTYPECTL;  //want to kill Myeloma cell
									  //printf("Stage G01: %d  %d  %d  %f  %f\n",a->PosX,a->PosY,a->PosZ,ECMHigh[a->PosX][a->PosY][a->PosZ],ECM3[a->PosX][a->PosY][a->PosZ]);
								  }
								  
								  ECMHigh[LivingListpo1->CTLLx][LivingListpo1->CTLLy][LivingListpo1->CTLLz]=NOCELL;								  
								  LivingListpo1->CTLLx=a->PosX;
								  LivingListpo1->CTLLy=a->PosY;
								  LivingListpo1->CTLLz=a->PosZ;
								  ECMHigh[a->PosX][a->PosY][a->PosZ]=CTYPECTL;								  
								  LivingListpo1->CTLCycleTimer=LivingListpo1->CTLCycleTimer+1;// after one time step, should plus 1
								  LivingListpo1++;
							  }
							  else // There are no free spaces for migration,stay there
							  {
								  LivingListpo1->CTLCycleTimer=LivingListpo1->CTLCycleTimer+1;
								  LivingListpo1++;
							  }
							  break;
 
						  }
					  case StaS: // in the S stage
						  {							  
							  Location *a;
							  double micd=DPC;
							  a=MSearch_ctl(LivingListpo1->CTLLx,LivingListpo1->CTLLy,LivingListpo1->CTLLz,micd);  //will change
							  if (a!=null)
							  {
								  if ((ECMHigh[a->PosX][a->PosY][a->PosZ]==CTYPEMM) || (ECMHigh[a->PosX][a->PosY][a->PosZ]==CTYPEMIC))
								  {									 
									  ECM3[a->PosX][a->PosY][a->PosZ]=CTYPECTL;  //want to kill Myeloma cell	
									  //printf("stage s:%d  %d  %d  %f  %f\n",a->PosX,a->PosY,a->PosZ,ECMHigh[a->PosX][a->PosY][a->PosZ],ECM3[a->PosX][a->PosY][a->PosZ]);
								  }

								  ECMHigh[LivingListpo1->CTLLx][LivingListpo1->CTLLy][LivingListpo1->CTLLz]=NOCELL;
								  LivingListpo1->CTLLx=a->PosX;
								  LivingListpo1->CTLLy=a->PosY;
								  LivingListpo1->CTLLz=a->PosZ;
								  ECMHigh[a->PosX][a->PosY][a->PosZ]=CTYPECTL;								  
								  LivingListpo1->CTLCycleTimer=LivingListpo1->CTLCycleTimer+1;
								  LivingListpo1++;
							  }
							  else
							  {
								  LivingListpo1->CTLCycleTimer=LivingListpo1->CTLCycleTimer+1;
								  LivingListpo1++;
							  }
							  break;


						  }
					  case StaG2:  // in the G2 stage
						  {							 
							  Location *a;
							  double micd=DPC;
							  a=MSearch_ctl(LivingListpo1->CTLLx,LivingListpo1->CTLLy,LivingListpo1->CTLLz,micd);
							  if (a!=null)
							  {
								  if ((ECMHigh[a->PosX][a->PosY][a->PosZ]==CTYPEMM) || (ECMHigh[a->PosX][a->PosY][a->PosZ]==CTYPEMIC))
								  {									
									  ECM3[a->PosX][a->PosY][a->PosZ]=CTYPECTL;  //want to kill Myeloma cell
									  //printf("stage G2: %d  %d  %d  %f  %f\n",a->PosX,a->PosY,a->PosZ,ECMHigh[a->PosX][a->PosY][a->PosZ],ECM3[a->PosX][a->PosY][a->PosZ]);
								  }

								  ECMHigh[LivingListpo1->CTLLx][LivingListpo1->CTLLy][LivingListpo1->CTLLz]=NOCELL;
								  LivingListpo1->CTLLx=a->PosX;
								  LivingListpo1->CTLLy=a->PosY;
								  LivingListpo1->CTLLz=a->PosZ;
								  ECMHigh[a->PosX][a->PosY][a->PosZ]=CTYPECTL;
								  LivingListpo1->CTLCycleTimer=LivingListpo1->CTLCycleTimer+1;
								  LivingListpo1++;
							  }
							  else
							  {
								  LivingListpo1->CTLCycleTimer=LivingListpo1->CTLCycleTimer+1;
								  LivingListpo1++;
							  }
							  break;

						  }
					  case StaM:
						  {
							  LivingListpo1->CTLCycleTimer=CYCLEG01BEGIN; // first reset the Cell Cycle Timer to CYCLEG01BEGIN.
							/* ===Search whether there are free spaces */
							  Location *a;
							  int r=radx;							 

							  a=ProSearch_ctl(LivingListpo1->CTLLx,LivingListpo1->CTLLy,LivingListpo1->CTLLz,r); 
							 
							  if (a!=null)
								  LivingListpo1->CTLFreeSpace=FREESPAY;
							  else
								  LivingListpo1->CTLFreeSpace=FREESPAN;
							  
							  if(LivingListpo1->CTLFreeSpace==FREESPAY)
							  {
								  if (LivingListpo1->CTLGenNo > MAXCTL)
								  {								
			                                                CTLApopList.insert(CTLApopList.end(),*LivingListpo1);
                           						LivingListpo1=CTLLivingList.erase(LivingListpo1);//delete CTL
								  }
								  else
								  {
									CTL c(i1++,a->PosX,a->PosY,a->PosZ);
							                c.CTLCurrentTimer=Tstep;// to ensure that the new cell c will not be checked in this time step
									ECMHigh[a->PosX][a->PosY][a->PosZ]=CTYPECTL;									
									CTLLivingList.insert(CTLLivingList.end(),c);
									
									LivingListpo1->CTLGenNo=LivingListpo1->CTLGenNo+1;
									LivingListpo1++;
								  }

							  }
							  else // there are no free spaces for proliferation, stay there
							  {								  
								  LivingListpo1->CTLCellStatus=StaM;
								  LivingListpo1->CTLCycleTimer=CYCLEMBEGIN;
								  LivingListpo1++;
							  }							  
						  }
				  }
			  }
			  else // Not cell cycle on //
			  {	
				  //--prolifeartion rate--------------------------------------------------------
				  double len0,pro_len,proratio,pprol;
				  double ratio_C_T,size1,size2;
				  double prodyn,tmp;
				  double tgfbvalue;				
				  double pro_treg;
				  size1=TregLivingList.size();
				  size2=CTLLivingList.size(); 
				  tgfbvalue=TGFb[LivingListpo1->CTLLx][LivingListpo1->CTLLy][LivingListpo1->CTLLz];
				  ratio_C_T=size1/size2;
				  len0=0.7152*codrug/(codrug+0.5217);   //LEN
				  proratio=0.82-3.03*((ratio_C_T-0.0659)/3.5609)/(1+((ratio_C_T-0.0659)/3.5609));   //the proliferation rate induced by the ratio of Treg/CD8
				  if (proratio<0.1)
					  proratio=0.1;
				  
                                  if (par::sw_treg==1)
				       pro_treg=pow(proratio/Td_Treg,3)/(1+pow(proratio/Td_Treg,3));  //the effect of Treg
				  else
                                       pro_treg=pow(proratio/Td_off_Treg,3)/(1+pow(proratio/Td_off_Treg,3));

				  pro_len=pow(len0/Td_LEN,4)/(1+pow(len0/Td_LEN,4));  //LEN
				  tmp=pow(tgfbvalue/CTLTGFB0,2); //the contribution of TGFb
				  prodyn=(1+pro_len)*pro_treg*(1/(1+tmp));  //Three elements: LEN, Treg/CD8, TGFb
				  pprol=CTLproa+CTLprob*prodyn;
				 
				 // printf("ratio:%f,Pro_treg:%f, Len:%f, Pro_len:%f, prodyn:%f, P:%f, A:%f\n",ratio_C_T, pro_treg, len0, pro_len,prodyn, pprol,apop);
				  //--------------------------------------------------------------------------------
				  double dr =getRandom();
				  //printf(" %f  %f\n",dr,pprol);
				  //getchar();
				  if (dr>=0 && dr<pprol)  //if can enter to cell cycle 
				  {
					  LivingListpo1->CTLCycleTimer=CYCLEG01BEGIN;
					  LivingListpo1->CTLCellStatus=StaG01;
					  LivingListpo1++;
				  }
				  else //continute to migrate
				  {
					  Location *a;
					  double micd=DPC;
					  a=MSearch_ctl(LivingListpo1->CTLLx,LivingListpo1->CTLLy,LivingListpo1->CTLLz,micd);  //will change
					  if (a!=null)
					  {						  
						  if ((ECMHigh[a->PosX][a->PosY][a->PosZ]==CTYPEMM) || (ECMHigh[a->PosX][a->PosY][a->PosZ]==CTYPEMIC))  //will directly lysis myeloma cell
						  {
						       ECM3[a->PosX][a->PosY][a->PosZ]=CTYPECTL;  //want to kill Myeloma cell	
							   //printf("No cell cycle: %d  %d  %d  %f  %f\n",a->PosX,a->PosY,a->PosZ,ECMHigh[a->PosX][a->PosY][a->PosZ],ECM3[a->PosX][a->PosY][a->PosZ]);
						  }
						 
						  ECMHigh[LivingListpo1->CTLLx][LivingListpo1->CTLLy][LivingListpo1->CTLLz]=NOCELL;						 
						  LivingListpo1->CTLLx=a->PosX;
						  LivingListpo1->CTLLy=a->PosY;
						  LivingListpo1->CTLLz=a->PosZ;
						  
						  ECMHigh[a->PosX][a->PosY][a->PosZ]=CTYPECTL;						  
						  LivingListpo1->CTLCycleTimer=LivingListpo1->CTLCycleTimer+1;
						  	
						  LivingListpo1++;						  
						  					  
					  }
					  else
					  {
						  LivingListpo1->CTLCycleTimer=LivingListpo1->CTLCycleTimer+1;						  						 
						  LivingListpo1++;
					  }
				  }		  
			  }
			  }
		  }
	  
	  }
	  else
		  LivingListpo1++;
  }

  for (ApopListpo1=CTLApopList.begin();ApopListpo1!=CTLApopList.end();)// to search Apoplist
  {
	  int Timer1=ApopListpo1->CTLCurrentTimer+1;//every loop, the CurrentTimer is the same to the TimeStep
	  if (Timer1==Tstep)// if the cell is inserted from the above process, then its PCCurrentTimer+1!=Tstep
	                                     // therefore it will not be considered in this process
	  {
		   ApopListpo1->CTLApopTimer = ApopListpo1->CTLApopTimer+1;
		   ApopListpo1->CTLCurrentTimer=ApopListpo1->CTLCurrentTimer+1;
		  if (ApopListpo1->CTLApopTimer >=5)
		  {
			  ApopListpo1->CTLApopStatus=APOPA;
			  ECMHigh[ApopListpo1->CTLLx][ApopListpo1->CTLLy][ApopListpo1->CTLLz]=NOCELL;			  
			  CTLAbsorbList.insert(CTLAbsorbList.end(),*ApopListpo1);
			  ApopListpo1=CTLApopList.erase(ApopListpo1);
		  }
		  else
			  ApopListpo1++;
	  }
	  else
		  ApopListpo1++;
  }
}

void MicroEnv::TregCellCheck(int Tstep, int drugindex)
{
  list <Treg>::iterator LivingListpo1;
  list <Treg>::iterator ApopListpo1;
  
  for (LivingListpo1=TregLivingList.begin();LivingListpo1!=TregLivingList.end();)
  {
	  int TempTimer=LivingListpo1->TregCurrentTimer+1;//every loop, the CurrentTimer is one less than the TimeStep	  
	  if (TempTimer==Tstep) // to guarantee that every cell is dealt with only once in every timestep.
	  {
		  LivingListpo1->TregCurrentTimer=LivingListpo1->TregCurrentTimer+1;// Here just plus 1 inside the loop for considering about the new PC cells that added in this timestep
		  if(LivingListpo1->TregApopStatus==APOPY)// if the cell is beginning to die
		  {
			  LivingListpo1->TregApopTimer=LivingListpo1->TregApopTimer+1;// after 5 timesteps, put it to HHAbsorbList
			  TregApopList.insert(TregApopList.end(),*LivingListpo1); // put it to being apoptosis cell list
			  LivingListpo1=TregLivingList.erase(LivingListpo1); // then delete it from the living cell list.
		  }
		  else // the cell is not apoptosis
		  {
			  // To judge if apoptosis rate < threshold
			  double thresh=getRandom(); // generate threshold
			  double cur_tgfbval,apoptgfb,len,apop,tmp;
			  cur_tgfbval=TGFb[LivingListpo1->TregLx][LivingListpo1->TregLy][LivingListpo1->TregLz];
			  apoptgfb=1/(1+pow(cur_tgfbval/CTLTGFB0,2));

			  tmp=Treg_LEN;
			  tmp=pow((codrug/tmp), 3);
			  len=tmp/(1+tmp);
			  apop=Tregapop0+TregTGFBapop*(1+len)*apoptgfb;	//TGFb,len may reduce apoptosis of Treg
			  			  
			  //printf("AP_T:%f,AP_L:%f,AP:%f\n", apoptgfb,len,apop);

			  if (apop>thresh) // if apoptosis rate < Threshold  ===> apoptosis
			  {
				  LivingListpo1->TregApopTimer=LivingListpo1->TregApopTimer+1;// after 5 timesteps, put it to MICAbsorbList
				  TregApopList.insert(TregApopList.end(),*LivingListpo1); // put it to being apoptosis cell list
				  LivingListpo1=TregLivingList.erase(LivingListpo1); // then delete it from the living cell list.
			  }
			  else // apoptosis rate > Threshold, go on following the flowchart
			  {
			  	 
			  if (LivingListpo1->TregCellStatus!=StaN)// enter the cell cycle
			  {
				  LivingListpo1->TregCellStatus=JudgeStatus(LivingListpo1->TregCycleTimer);
				  //printf("Treg, ID: %d  apopotosis rate:%f   %d\n",LivingListpo1->ID, apop,LivingListpo1->TregCellStatus);
				  switch (LivingListpo1->TregCellStatus)// Judge which stage the cell is in cell cycle
				  {
					  case StaG01: // in the G0/G1 stage
						  {
							  Location *a;  
							  double micd=DPC;
							  a=MSearch_treg(LivingListpo1->TregLx,LivingListpo1->TregLy,LivingListpo1->TregLz,micd);  //will change
							  if (a!=null)// There are free space to migrate
							  {						
								  if (ECMHigh[a->PosX][a->PosY][a->PosZ]==NOCELL)
								  {
									  ECMHigh[LivingListpo1->TregLx][LivingListpo1->TregLy][LivingListpo1->TregLz]=NOCELL;								  
								      LivingListpo1->TregLx=a->PosX;
								      LivingListpo1->TregLy=a->PosY;
								      LivingListpo1->TregLz=a->PosZ;
								      ECMHigh[a->PosX][a->PosY][a->PosZ]=CTYPETREG;							  
								   }
								  else
								  {
									  ECM3[a->PosX][a->PosY][a->PosZ]=CTYPETREG;
								  }
								  LivingListpo1->TregCycleTimer=LivingListpo1->TregCycleTimer+1;// after one time step, should plus 1
								  LivingListpo1++;
								  
							  }
							  else // There are no free spaces for migration,stay there
							  {
								  LivingListpo1->TregCycleTimer=LivingListpo1->TregCycleTimer+1;
								  LivingListpo1++;
							  }
							  break;
 
						  }
					  case StaS: // in the S stage
						  {							  
							  Location *a;
							  double micd=DPC;
							  a=MSearch_treg(LivingListpo1->TregLx,LivingListpo1->TregLy,LivingListpo1->TregLz,micd);  //will change
							  if (a!=null)
							  {					
								  if (ECMHigh[a->PosX][a->PosY][a->PosZ]==NOCELL)
								  {
									  ECMHigh[LivingListpo1->TregLx][LivingListpo1->TregLy][LivingListpo1->TregLz]=NOCELL;
								      LivingListpo1->TregLx=a->PosX;
								      LivingListpo1->TregLy=a->PosY;
								      LivingListpo1->TregLz=a->PosZ;
								      ECMHigh[a->PosX][a->PosY][a->PosZ]=CTYPETREG;									      
								  }
								  else
								  {
									  ECM3[a->PosX][a->PosY][a->PosZ]=CTYPETREG;
								  }
								  LivingListpo1->TregCycleTimer=LivingListpo1->TregCycleTimer+1;
								  LivingListpo1++;
							  }
							  else
							  {
								  LivingListpo1->TregCycleTimer=LivingListpo1->TregCycleTimer+1;
								  LivingListpo1++;
							  }
							  break;


						  }
					  case StaG2:  // in the G2 stage
						  {							 
							  Location *a;
							  double micd=DPC;
							  a=MSearch_treg(LivingListpo1->TregLx,LivingListpo1->TregLy,LivingListpo1->TregLz,micd);
							  if (a!=null)
							  {								  
								  if (ECMHigh[a->PosX][a->PosY][a->PosZ]==NOCELL)
								  {
									  ECMHigh[LivingListpo1->TregLx][LivingListpo1->TregLy][LivingListpo1->TregLz]=NOCELL;
								      LivingListpo1->TregLx=a->PosX;
								      LivingListpo1->TregLy=a->PosY;
								      LivingListpo1->TregLz=a->PosZ;
								      ECMHigh[a->PosX][a->PosY][a->PosZ]=CTYPETREG;									     
								  }
								  else
								  {
								      ECM3[a->PosX][a->PosY][a->PosZ]=CTYPETREG;
								  }
								  LivingListpo1->TregCycleTimer=LivingListpo1->TregCycleTimer+1;
								  LivingListpo1++;
							  }
							  else
							  {
								  LivingListpo1->TregCycleTimer=LivingListpo1->TregCycleTimer+1;
								  LivingListpo1++;
							  }
							  break;

						  }
					  case StaM:
						  {
							  LivingListpo1->TregCycleTimer=CYCLEG01BEGIN; // first reset the Cell Cycle Timer to CYCLEG01BEGIN.
							/* ===Search whether there are free spaces */
							  Location *a;
							  int r=radx;							 
							  
							  a=ProSearch_treg(LivingListpo1->TregLx,LivingListpo1->TregLy,LivingListpo1->TregLz,r); 
							  
							  if (a!=null)
								  LivingListpo1->TregFreeSpace=FREESPAY;
							  else
								  LivingListpo1->TregFreeSpace=FREESPAN;
							  
							  if(LivingListpo1->TregFreeSpace==FREESPAY)
							  {
								  if (LivingListpo1->TregGenNo > MAXTreg)
								  {								
									LivingListpo1=TregLivingList.erase(LivingListpo1);//delete Treg
								  }
								  else
								  {
									Treg c(i1++,a->PosX,a->PosY,a->PosZ);
							        c.TregCurrentTimer=Tstep;// to ensure that the new cell c will not be checked in this time step
									ECMHigh[a->PosX][a->PosY][a->PosZ]=CTYPETREG;									
									TregLivingList.insert(TregLivingList.end(),c);
									
									LivingListpo1->TregGenNo=LivingListpo1->TregGenNo+1;
									LivingListpo1++;
								  }

							  }
							  else // there are no free spaces for proliferation, stay there
							  {								  
								  LivingListpo1->TregCellStatus=StaM;
								  LivingListpo1->TregCycleTimer=CYCLEMBEGIN;
								  LivingListpo1++;
							  }							  
						  }
				  }
			  }
			  else // Not cell cycle on //
			  {	
				  //--prolifeartion rate--------------------------------------------------------
				  double prolen,proratio,pprol;
				  double tgfbvalue,protgfb,tmp,tmp2;		
				  
				  tgfbvalue=TGFb[LivingListpo1->TregLx][LivingListpo1->TregLy][LivingListpo1->TregLz];				  
				  tmp=pow(tgfbvalue/CTLTGFB0,2);
				  protgfb=1+tmp/(1+tmp);
				  				  
				  tmp2=Treg_LEN;
				  prolen=1/(1+pow(codrug/tmp2,3)); 
				  pprol=Tregpro0+Tregprotgfb*protgfb*prolen;
				  
				  //printf("Pr_T:%f,Pr_L:%f,pro:%f Ap_L:%f, Ap:%f \n",protgfb, prolen,pprol,len, apop);
				  //--------------------------------------------------------------------------------
				  double dr =getRandom();				  
				  
				  if (dr>=0 && dr<pprol)  //if can enter to cell cycle 
				  {
					  LivingListpo1->TregCycleTimer=CYCLEG01BEGIN;
					  LivingListpo1->TregCellStatus=StaG01;
					  LivingListpo1++;
				  }
				  else //continute to migrate
				  {
					  Location *a;
					  double micd=DPC;
					  a=MSearch_treg(LivingListpo1->TregLx,LivingListpo1->TregLy,LivingListpo1->TregLz,micd);  //will change
					  if (a!=null)
					  {						  						 					 
						 if ((ECMHigh[a->PosX][a->PosY][a->PosZ]==NOCELL))
						  {
							  ECMHigh[LivingListpo1->TregLx][LivingListpo1->TregLy][LivingListpo1->TregLz]=NOCELL;
							  LivingListpo1->TregLx=a->PosX;
							  LivingListpo1->TregLy=a->PosY;
							  LivingListpo1->TregLz=a->PosZ;
							  ECMHigh[a->PosX][a->PosY][a->PosZ]=CTYPETREG;								 
						   }
						 else
						  {
						      ECM3[a->PosX][a->PosY][a->PosZ]=CTYPETREG;
						  }
						  LivingListpo1->TregCycleTimer=LivingListpo1->TregCycleTimer+1;
						  LivingListpo1++;						  					  
					  }
					  else
					  {
						  LivingListpo1->TregCycleTimer=LivingListpo1->TregCycleTimer+1;						  						 
						  LivingListpo1++;
					  }
				  }		  
			  }
			  }
		  }
	  
	  }
	  else
		  LivingListpo1++;
  }

  for (ApopListpo1=TregApopList.begin();ApopListpo1!=TregApopList.end();)// to search Apoplist
  {
	  int Timer1=ApopListpo1->TregCurrentTimer+1;//every loop, the CurrentTimer is the same to the TimeStep
	  if (Timer1==Tstep)// if the cell is inserted from the above process, then its PCCurrentTimer+1!=Tstep
	                                     // therefore it will not be considered in this process
	  {
		   ApopListpo1->TregApopTimer = ApopListpo1->TregApopTimer+1;
		   ApopListpo1->TregCurrentTimer=ApopListpo1->TregCurrentTimer+1;
		  if (ApopListpo1->TregApopTimer >=5)
		  {
			  ApopListpo1->TregApopStatus=APOPA;
			  ECMHigh[ApopListpo1->TregLx][ApopListpo1->TregLy][ApopListpo1->TregLz]=NOCELL;			  
			  TregAbsorbList.insert(TregAbsorbList.end(),*ApopListpo1);
			  ApopListpo1=TregApopList.erase(ApopListpo1);
		  }
		  else
			  ApopListpo1++;
	  }
	  else
		  ApopListpo1++;
  }
}



void MicroEnv::PCCellCheck(int Tstep, int drugindex)
{
  //list <PC>::iterator LivingListpo1;
  //list <PC>::iterator ApopListpo1;

  //for (LivingListpo1=PCLivingList.begin();LivingListpo1!=PCLivingList.end();)
  //{
	 // int TempTimer=LivingListpo1->PCCurrentTimer+1;//every loop, the CurrentTimer is one less than the TimeStep	  
	 // if (TempTimer==Tstep) // to guarantee that every cell is dealt with only once in every timestep.
	 // {
		//  LivingListpo1->PCCurrentTimer=LivingListpo1->PCCurrentTimer+1;// Here just plus 1 inside the loop for considering about the new PC cells that added in this timestep
		//  if(LivingListpo1->PCApopStatus==APOPY)// if the cell is beginning to die
		//  {
		//	  LivingListpo1->PCApopTimer=LivingListpo1->PCApopTimer+1;// after 5 timesteps, put it to HHAbsorbList
		//	  PCApopList.insert(PCApopList.end(),*LivingListpo1); // put it to being apoptosis cell list
		//	  LivingListpo1=PCLivingList.erase(LivingListpo1); // then delete it from the living cell list.
		//  }
		//  else // the cell is not apoptosis
		//  {
		//	  // To judge if apoptosis rate < threshold
		//	  double thresh=getRandom(); // generate threshold
		//	  double e0,kre,emax,apop,ka,AA,hill1,hill2,apopdrug;
		//	  apopdrug=PCapopdrug;  e0=E0;
		//	  emax=EMAX;  kre=(emax-e0)/2;
		//	  apop=PCapop0;  e0=E0;
		//	  ka=K_A;
		//	  if ( (Tstep<100) | (Tstep >= 300))
		//		  AA=0;
		//	  else
		//	  {
		//		  if (drugindex==-1)
		//			  AA=0;
		//		  else
		//			  AA=drug[drugindex]; // AA=0: no drug; AA=1: drug
		//	  }
		//	  hill1=AA/ka;  hill1=pow(hill1,2);
		//	  hill1=hill1/(1+hill1);
		//	  hill2=STIFF[LivingListpo1->PCLx][LivingListpo1->PCLy][LivingListpo1->PCLz];
		//	  hill2=hill2-e0;
		//	  hill2=hill2/kre;
		//	  hill2=pow(hill2,2);
		//	  hill2+=1;
		//	  hill2=1/hill2;
		//	  apopdrug=apopdrug*hill1;
		//	  apopdrug=apopdrug*hill2;	
		//	  apop=apop+apopdrug;	
		//	  if (apop>thresh) // if apoptosis rate < Threshold  ===> apoptosis
		//	  {
		//		  LivingListpo1->PCApopTimer=LivingListpo1->PCApopTimer+1;// after 5 timesteps, put it to MICAbsorbList
		//		  //LivingListpo1->MICApopStatus=APOPY;
		//		  PCApopList.insert(PCApopList.end(),*LivingListpo1); // put it to being apoptosis cell list
		//		  LivingListpo1=PCLivingList.erase(LivingListpo1); // then delete it from the living cell list.
		//	  }
		//	  else // apoptosis rate > Threshold, go on following the flowchart
		//	  {

		//	  if (LivingListpo1->PCCellStatus!=StaN)// enter the cell cycle
		//	  {
		//		  LivingListpo1->PCCellStatus=JudgeStatus(LivingListpo1->PCCycleTimer);
		//		  switch (LivingListpo1->PCCellStatus)// Judge which stage the cell is in cell cycle
		//		  {
		//			  case StaG01: // in the G0/G1 stage
		//				  {
		//					  // after five time steps-> one time step is 2 hours
		//						 //* migration1 *//
		//					  Location *a;  
		//					  double micd=DPC;
		//					  a=MSearch1(LivingListpo1->PCLx,LivingListpo1->PCLy,LivingListpo1->PCLz,micd);
		//					  if (a!=null)// There are free space to migrate
		//					  {
		//						  ECMHigh[LivingListpo1->PCLx][LivingListpo1->PCLy][LivingListpo1->PCLz]=NOCELL;
		//						  ECMCan[LivingListpo1->PCLx][LivingListpo1->PCLy][LivingListpo1->PCLz]=NOCELL;// begin to move, the current position will be empty

		//						  LivingListpo1->PCLx=a->PosX;
		//						  LivingListpo1->PCLy=a->PosY;
		//						  LivingListpo1->PCLz=a->PosZ;

		//						  ECMHigh[a->PosX][a->PosY][a->PosZ]=CTYPEPC;
		//						  ECMCan[a->PosX][a->PosY][a->PosZ]=CTYPEPC;

		//						  TGFb[a->PosX][a->PosY][a->PosZ]=TGFb[a->PosX][a->PosY][a->PosZ]+PCTGFb;
		//						  IL6[a->PosX][a->PosY][a->PosZ]=IL6[a->PosX][a->PosY][a->PosZ]+PCIL6;
		//						  //SDF1[a->PosX][a->PosY][a->PosZ]=SDF1[a->PosX][a->PosY][a->PosZ]+PCSDF1;
		//						  LivingListpo1->PCCycleTimer=LivingListpo1->PCCycleTimer+1;// after one time step, should plus 1
		//						  LivingListpo1++;
		//					  }
		//					  else // There are no free spaces for migration
		//					  {
		//						  LivingListpo1->PCCycleTimer=LivingListpo1->PCCycleTimer+1;
		//						  TGFb[LivingListpo1->PCLx][LivingListpo1->PCLy][LivingListpo1->PCLz]=TGFb[LivingListpo1->PCLx][LivingListpo1->PCLy][LivingListpo1->PCLz]+PCTGFb;
		//						  IL6[LivingListpo1->PCLx][LivingListpo1->PCLy][LivingListpo1->PCLz]=IL6[LivingListpo1->PCLx][LivingListpo1->PCLy][LivingListpo1->PCLz]+PCIL6;
		//						  //SDF1[LivingListpo1->PCLx][LivingListpo1->PCLy][LivingListpo1->PCLz]=SDF1[LivingListpo1->PCLx][LivingListpo1->PCLy][LivingListpo1->PCLz]+PCSDF1;
		//						  LivingListpo1++;
		//					  }

		//					 break;
 
		//				  }
		//			  case StaS: // in the S stage
		//				  {
		//					  //* migration1 *//
		//					  /* ===To set NEW LOCATION, migD is given*/
		//					  Location *a;
		//					  double micd=DPC;
		//					  a=MSearch1(LivingListpo1->PCLx,LivingListpo1->PCLy,LivingListpo1->PCLz,micd);
		//					  if (a!=null)
		//					  {
		//						  ECMHigh[LivingListpo1->PCLx][LivingListpo1->PCLy][LivingListpo1->PCLz]=NOCELL;
		//						  ECMCan[LivingListpo1->PCLx][LivingListpo1->PCLy][LivingListpo1->PCLz]=NOCELL;
		//						  LivingListpo1->PCLx=a->PosX;
		//						  LivingListpo1->PCLy=a->PosY;
		//						  LivingListpo1->PCLz=a->PosZ;
		//						  ECMHigh[a->PosX][a->PosY][a->PosZ]=CTYPEPC;
		//						  ECMCan[a->PosX][a->PosY][a->PosZ]=CTYPEPC;
		//						  TGFb[a->PosX][a->PosY][a->PosZ]=TGFb[a->PosX][a->PosY][a->PosZ]+PCTGFb;
		//						  IL6[a->PosX][a->PosY][a->PosZ]=IL6[a->PosX][a->PosY][a->PosZ]+PCIL6;
		//						  //SDF1[a->PosX][a->PosY][a->PosZ]=SDF1[a->PosX][a->PosY][a->PosZ]+PCSDF1;
		//						  LivingListpo1->PCCycleTimer=LivingListpo1->PCCycleTimer+1;
		//						  LivingListpo1++;
		//					  }
		//					  else
		//					  {
		//						  LivingListpo1->PCCycleTimer=LivingListpo1->PCCycleTimer+1;
		//						  TGFb[LivingListpo1->PCLx][LivingListpo1->PCLy][LivingListpo1->PCLz]=TGFb[LivingListpo1->PCLx][LivingListpo1->PCLy][LivingListpo1->PCLz]+PCTGFb;
		//						  IL6[LivingListpo1->PCLx][LivingListpo1->PCLy][LivingListpo1->PCLz]=IL6[LivingListpo1->PCLx][LivingListpo1->PCLy][LivingListpo1->PCLz]+PCIL6;
		//						  //SDF1[LivingListpo1->PCLx][LivingListpo1->PCLy][LivingListpo1->PCLz]=SDF1[LivingListpo1->PCLx][LivingListpo1->PCLy][LivingListpo1->PCLz]+PCSDF1;
		//						  LivingListpo1++;
		//					  }
		//					  break;


		//				  }
		//			  case StaG2:  // in the G2 stage
		//				  {
		//					  //* migration1 *//
		//					  Location *a;
		//					  double micd=DPC;
		//					  a=MSearch1(LivingListpo1->PCLx,LivingListpo1->PCLy,LivingListpo1->PCLz,micd);
		//					  if (a!=null)
		//					  {
		//						  ECMHigh[LivingListpo1->PCLx][LivingListpo1->PCLy][LivingListpo1->PCLz]=NOCELL;
		//						  ECMCan[LivingListpo1->PCLx][LivingListpo1->PCLy][LivingListpo1->PCLz]=NOCELL;
		//						  LivingListpo1->PCLx=a->PosX;
		//						  LivingListpo1->PCLy=a->PosY;
		//						  LivingListpo1->PCLz=a->PosZ;
		//						  ECMHigh[a->PosX][a->PosY][a->PosZ]=CTYPEPC;
		//						   ECMCan[a->PosX][a->PosY][a->PosZ]=CTYPEPC;
		//						  TGFb[a->PosX][a->PosY][a->PosZ]=TGFb[a->PosX][a->PosY][a->PosZ]+PCTGFb;
		//						  //SDF1[a->PosX][a->PosY][a->PosZ]=SDF1[a->PosX][a->PosY][a->PosZ]+PCSDF1;
		//						  IL6[a->PosX][a->PosY][a->PosZ]=IL6[a->PosX][a->PosY][a->PosZ]+PCIL6;
		//						  LivingListpo1->PCCycleTimer=LivingListpo1->PCCycleTimer+1;
		//						  LivingListpo1++;
		//					  }
		//					  else
		//					  {
		//						  LivingListpo1->PCCycleTimer=LivingListpo1->PCCycleTimer+1;
		//						  TGFb[LivingListpo1->PCLx][LivingListpo1->PCLy][LivingListpo1->PCLz]+=PCTGFb;
		//						  //SDF1[LivingListpo1->PCLx][LivingListpo1->PCLy][LivingListpo1->PCLz]+=PCSDF1;
		//						  IL6[LivingListpo1->PCLx][LivingListpo1->PCLy][LivingListpo1->PCLz]+=PCIL6;
		//						  LivingListpo1++;
		//					  }
		//					  break;

		//				  }
		//			  case StaM:
		//				  {
		//					  LivingListpo1->PCCycleTimer=CYCLEG01BEGIN; // first reset the Cell Cycle Timer to CYCLEG01BEGIN.
		//						/* ===Search whether there are free spaces */
		//					  Location *a;
		//					  int r=radx;
		//					  a=ProSearch(LivingListpo1->PCLx,LivingListpo1->PCLy,LivingListpo1->PCLz,r);
		//					  if (a!=null)
		//						  LivingListpo1->PCFreeSpace=FREESPAY;
		//					  else
		//						  LivingListpo1->PCFreeSpace=FREESPAN;

		//					  if(LivingListpo1->PCFreeSpace==FREESPAY)
		//					  {
		//						  if (LivingListpo1->PCGenNo>MAXPC)
		//						  {
		//						    MM c1(i1++,a->PosX,a->PosY,a->PosZ);
		//							c1.MMCurrentTimer=Tstep;// to ensure that the new cell c will not be checked in this time step
		//							c1.MMGenNo=LivingListpo1->PCGenNo+1;
		//							ECMHigh[a->PosX][a->PosY][a->PosZ]=CTYPEMM;
		//							ECMCan[a->PosX][a->PosY][a->PosZ]=CTYPEMM;
		//							TGFb[a->PosX][a->PosY][a->PosZ]+=MMTGFb;
		//							IL6[a->PosX][a->PosY][a->PosZ]+=MMIL6;
		//							//SDF1[a->PosX][a->PosY][a->PosZ]+=MMSDF1;
		//							
		//							MMLivingList.insert(MMLivingList.end(),c1);
		//							
		//							// this is added on 8/25/2011
		//							// original PC turns to MM also

		//							MM c2(LivingListpo1->ID,LivingListpo1->PCLx,LivingListpo1->PCLy,LivingListpo1->PCLz);
		//							c2.MMCurrentTimer=Tstep;// to ensure that the cell changed from mother cell will not be checked in this time step
		//							c2.MMGenNo=LivingListpo1->PCGenNo+1;
		//							ECMHigh[LivingListpo1->PCLx][LivingListpo1->PCLy][LivingListpo1->PCLz]=CTYPEMM;
		//							ECMCan[LivingListpo1->PCLx][LivingListpo1->PCLy][LivingListpo1->PCLz]=CTYPEMM;
		//							TGFb[LivingListpo1->PCLx][LivingListpo1->PCLy][LivingListpo1->PCLz]+=MMTGFb;
		//							IL6[LivingListpo1->PCLx][LivingListpo1->PCLy][LivingListpo1->PCLz]+=MMIL6;
		//							//SDF1[LivingListpo1->PCLx][LivingListpo1->PCLy][LivingListpo1->PCLz]+=MMSDF1;
		//							
		//							MMLivingList.insert(MMLivingList.end(),c2);
		//							// the mother cell changed to MM , so it has to be removed from the PC living list
		//							//LivingListpo1->PCGenNo=LivingListpo1->PCGenNo+1;
		//							LivingListpo1=PCLivingList.erase(LivingListpo1);//delete PC
		//						  }
		//						  else
		//						  {
		//						    PC c(i1++,a->PosX,a->PosY,a->PosZ);
		//							c.PCCurrentTimer=Tstep;// to ensure that the new cell c will not be checked in this time step
		//							c.PCGenNo=LivingListpo1->PCGenNo+1;
		//							ECMHigh[a->PosX][a->PosY][a->PosZ]=CTYPEPC;
		//							ECMCan[a->PosX][a->PosY][a->PosZ]=CTYPEPC;
		//							TGFb[a->PosX][a->PosY][a->PosZ]+=PCTGFb;
		//							IL6[a->PosX][a->PosY][a->PosZ]+=PCIL6;
		//							//SDF1[a->PosX][a->PosY][a->PosZ]+=PCSDF1;
		//							PCLivingList.insert(PCLivingList.end(),c);


		//							TGFb[LivingListpo1->PCLx][LivingListpo1->PCLy][LivingListpo1->PCLz]+=PCTGFb;
		//							IL6[LivingListpo1->PCLx][LivingListpo1->PCLy][LivingListpo1->PCLz]+=PCIL6;
		//							//SDF1[LivingListpo1->PCLx][LivingListpo1->PCLy][LivingListpo1->PCLz]+=PCSDF1;
		//							LivingListpo1->PCGenNo=LivingListpo1->PCGenNo+1;
		//							LivingListpo1++;

		//						  }

		//					  }
		//					  else // there are no free spaces for pro, 
		//					  {
		//						  //if(LivingListpo1->PCWaitProTimer>=T1) //T1 is the maximum pro time
		//						  //{
		//							 // LivingListpo1->PCApopTimer=1;// after 5 timesteps, put it to HHAbsorbList
		//							 // LivingListpo1->PCApopStatus=APOPY;
		//							 // PCApopList.insert(PCApopList.end(),*LivingListpo1); // put it to being apoptosis cell list
		//							 // LivingListpo1=PCLivingList.erase(LivingListpo1); // then delete it from the living cell list.
		//						  //}
		//						  //else
		//						  //{
		//							 // LivingListpo1->PCWaitProTimer=LivingListpo1->PCWaitProTimer+1;
		//							 // LivingListpo1->PCCellStatus=StaM;
		//							 // LivingListpo1->PCCycleTimer=CYCLEMBEGIN; // to ensure the cell is still in M Status for the next time step
		//							 // LivingListpo1++;
		//						  //}
		//						  LivingListpo1->PCCellStatus=StaM;
		//						  LivingListpo1->PCCycleTimer=CYCLEMBEGIN;
		//						  LivingListpo1++;
		//					  }
		//					  
		//				  }
		//		  }
		//	  }
		//	  else // Not cell cycle on // ======= THE PATHWAY FUNCTION SHOULD BE CHANGED!!!!!!!!
		//	  {
		//		  //Diff(D_TGFb,Kd_TGFb,theta_d,TGFb); // to get new values of TGFb at each grid point
		//		  
		//		  double pathdo=0;
		//		  double e0,emax,Ke,pathway,p0,pprol;
		//		  int count1=0;
		//		  e0=E0;
		//		  emax=EMAX;
		//		  Ke=(emax-e0)/2;
		//			for(int x=-1;x<=1;x++)
		//			{
		//			 for(int y=-1;y<=1;y++)
		//			  {
		//				for(int z=-1;z<=1;z++)
		//				{
		//					if( (BoundaryCheck(LivingListpo1->PCLx+x,LivingListpo1->PCLy+y,LivingListpo1->PCLz+z)==0)&& (ECMHigh[LivingListpo1->PCLx+x][LivingListpo1->PCLy+y][LivingListpo1->PCLz+z]==BMSCCELL) )	
		//				 {
		//				 //if (pathdo<STIFF[LivingListpo1->PCLx+x][LivingListpo1->PCLy+y][LivingListpo1->PCLz+z])
		//				 {
		//					 pathdo+=STIFF[LivingListpo1->PCLx+x][LivingListpo1->PCLy+y][LivingListpo1->PCLz+z];
		//					 count1++;

		//				 }
		//				 }
		//				}
		//			 }
		//		   }
	
		//		  if (count1!=0)
		//			  pathdo=pathdo/count1;

		//		  pathdo=(pathdo-e0)/Ke;
		//		  pathdo=pow(pathdo,2);
		//		  pathdo=pathdo/(1+pathdo);
		//		  p0=PCp0;
		//		  pathway=PCP_pathway;
		//		  pprol=p0+pathway*pathdo;

		//		  /*double WNTeff=(WNT2[LivingListpo1->PCLx][LivingListpo1->PCLy][LivingListpo1->PCLz]);
		//		  double wnt1=DKK2[LivingListpo1->PCLx][LivingListpo1->PCLy][LivingListpo1->PCLz]/WNTeff;
		//		  double c1=c_dkk;
		//		  wnt1=c1*wnt1+1;
		//		  WNTeff=WNTeff/wnt1;///(1+c_dkk*(DKK2[LivingListpo1->HHLx][LivingListpo1->HHLy])/(WNT2[LivingListpo1->HHLx][LivingListpo1->HHLy]));
		//		  double Kwnt=K_WNT;
		//		  wnt1=WNTeff/Kwnt;
		//		  wnt1=wnt1*wnt1;
		//		  double Ktgfb=K_TGFb;
		//		  double TGFb1=TGFb[LivingListpo1->PCLx][LivingListpo1->PCLy][LivingListpo1->PCLz]/Ktgfb;
		//		  TGFb1=TGFb1*TGFb1;
		//		  double pprol,pathway=OTP_pathway;
		//		  wnt1=wnt1/(1+wnt1+TGFb1);
		//		  pathway=pathway*wnt1;
		//		  pathway=pathway+OTp0;
		//		  pprol=pathway;
		//		  *///p_prol=0.5;
		//		  /*int t =(int) time(NULL) ; 
		//		  srand (t); 
		//		  double ran = rand(); */
		//		  double dr =getRandom();
		//		  //dr=0.4;
		//		  if (dr>=0 && dr<pprol)
		//		  {
		//			  LivingListpo1->PCCycleTimer=CYCLEG01BEGIN;
		//			  LivingListpo1->PCCellStatus=StaG01;
		//			  LivingListpo1++;
		//		  }
		//		  else
		//		  {
		//			  Location *a;
		//			  double micd=DPC;
		//			  a=MSearch1(LivingListpo1->PCLx,LivingListpo1->PCLy,LivingListpo1->PCLz,micd);
		//			  if (a!=null)
		//			  {
		//				  ECMHigh[LivingListpo1->PCLx][LivingListpo1->PCLy][LivingListpo1->PCLz]=NOCELL;
		//				  ECMCan[LivingListpo1->PCLx][LivingListpo1->PCLy][LivingListpo1->PCLz]=NOCELL;
		//				  LivingListpo1->PCLx=a->PosX;
		//				  LivingListpo1->PCLy=a->PosY;
		//				  LivingListpo1->PCLz=a->PosZ;
		//				  ECMHigh[a->PosX][a->PosY][a->PosZ]=CTYPEPC;
		//				  ECMCan[a->PosX][a->PosY][a->PosZ]=CTYPEPC;
		//				  LivingListpo1->PCCycleTimer=LivingListpo1->PCCycleTimer+1;
		//				  TGFb[a->PosX][a->PosY][a->PosZ]+=PCTGFb;
		//				  //SDF1[a->PosX][a->PosY][a->PosZ]+=PCSDF1;	
		//				  IL6[a->PosX][a->PosY][a->PosZ]+=PCIL6;
		//				  LivingListpo1++;
		//				  //LivingListpo1->HHCurrentTimer=LivingListpo1->HHCurrentTimer+1;
		//			  }
		//			  else
		//			  {
		//				  LivingListpo1->PCCycleTimer=LivingListpo1->PCCycleTimer+1;
		//				  TGFb[LivingListpo1->PCLx][LivingListpo1->PCLy][LivingListpo1->PCLz]+=PCTGFb;
		//				  IL6[LivingListpo1->PCLx][LivingListpo1->PCLy][LivingListpo1->PCLz]+=PCIL6;
		//				  //SDF1[LivingListpo1->PCLx][LivingListpo1->PCLy][LivingListpo1->PCLz]+=PCSDF1;
		//				  LivingListpo1++;
		//			  }

		//		  }
		//  
		//	  }
		//	  }
		//  }
	 // 
	 // }
	 // else
		//  LivingListpo1++;
  //}

  //for (ApopListpo1=PCApopList.begin();ApopListpo1!=PCApopList.end();)// to search Apoplist
  //{
	 // int Timer1=ApopListpo1->PCCurrentTimer+1;//every loop, the CurrentTimer is the same to the TimeStep
	 // if (Timer1==Tstep)// if the cell is inserted from the above process, then its PCCurrentTimer+1!=Tstep
	 //                                    // therefore it will not be considered in this process
	 // {
		//   ApopListpo1->PCApopTimer = ApopListpo1->PCApopTimer+1;
		//   ApopListpo1->PCCurrentTimer=ApopListpo1->PCCurrentTimer+1;
		//  if (ApopListpo1->PCApopTimer >=5)
		//  {
		//	  ApopListpo1->PCApopStatus=APOPA;
		//	  ECMHigh[ApopListpo1->PCLx][ApopListpo1->PCLy][ApopListpo1->PCLz]=NOCELL;
		//	  ECMCan[ApopListpo1->PCLx][ApopListpo1->PCLy][ApopListpo1->PCLz]=NOCELL;
		//	  PCAbsorbList.insert(PCAbsorbList.end(),*ApopListpo1);
		//	  ApopListpo1=PCApopList.erase(ApopListpo1);
		//  }
		//  else
		//	  ApopListpo1++;
	 // }
	 // else
		//  ApopListpo1++;
  //}
}

void MicroEnv::TMMCellCheck(int Tstep, int drugindex)
{
		//
  //list <TMM>::iterator LivingListpo1;
  //list <TMM>::iterator ApopListpo1;

  //for (LivingListpo1=TMMLivingList.begin();LivingListpo1!=TMMLivingList.end();)
  //{
	 // int TempTimer=LivingListpo1->TMMCurrentTimer+1;//every loop, the CurrentTimer is one less than the TimeStep	  
	 // if (TempTimer==Tstep) // to guarantee that every cell is dealt with only once in every timestep.
	 // {
		//  LivingListpo1->TMMCurrentTimer=LivingListpo1->TMMCurrentTimer+1;

		//  if(LivingListpo1->TMMApopStatus==APOPY)// if the cell is beginning to die
		//  {
		//	  LivingListpo1->TMMApopTimer=LivingListpo1->TMMApopTimer+1;// after 5 timesteps, put it to HHAbsorbList
		//	  TMMApopList.insert(TMMApopList.end(),*LivingListpo1); // put it to being apoptosis cell list
		//	  LivingListpo1=TMMLivingList.erase(LivingListpo1); // then delete it from the living cell list.
		//  }
		//  else // the cell is not apoptosis, then judge if apoptosis rate < threshold
		//  {
		//	  double thresh=getRandom(); // generate threshold
		//	  double e0,kre,emax,apop,ka,AA,hill1,hill2,apopdrug;
		//	  apopdrug=TMMapopdrug;  e0=E0;
		//	  emax=EMAX;  kre=(emax-e0)/2;
		//	  apop=TMMapop0;  e0=E0;
		//	  ka=K_A;
		//	  if ( (Tstep<100) | (Tstep >= 300))
		//		  AA=0;
		//	  else
		//	  {
		//		  if (drugindex==-1)
		//			  AA=0;
		//		  else
		//			  AA=drug[drugindex]; // AA=0: no drug; AA=1: drug
		//	  }
		//	  hill1=AA/ka;  hill1=pow(hill1,2);
		//	  hill1=hill1/(1+hill1);
		//	  hill2=STIFF[LivingListpo1->TMMLx][LivingListpo1->TMMLy][LivingListpo1->TMMLz];
		//	  hill2=hill2-e0;
		//	  hill2=hill2/kre;
		//	  hill2=pow(hill2,2);
		//	  hill2+=1;
		//	  hill2=1/hill2;
		//	  apopdrug=apopdrug*hill1;
		//	  apopdrug=apopdrug*hill2;	
		//	  apop=apop+apopdrug;	
		//	  if (apop>thresh) // if apoptosis rate < Threshold  ===> apoptosis
		//	  {
		//		  //LivingListpo1->MICApopStatus=APOPY;
		//		  LivingListpo1->TMMApopTimer=LivingListpo1->TMMApopTimer+1;// after 5 timesteps, put it to TMMAbsorbList
		//		  TMMApopList.insert(TMMApopList.end(),*LivingListpo1); // put it to being apoptosis cell list
		//		  LivingListpo1=TMMLivingList.erase(LivingListpo1); // then delete it from the living cell list.
		//	  }
		//	  else
		//	  {
		//	  	
		//  	  Location *a;
		//	  double micd=DTMM;
		//	  a=MSearch1(LivingListpo1->TMMLx,LivingListpo1->TMMLy,LivingListpo1->TMMLz,micd);
		//	  if (a!=null)
		//	  {
		//		  ECMHigh[LivingListpo1->TMMLx][LivingListpo1->TMMLy][LivingListpo1->TMMLz]=NOCELL;
		//		  ECMCan[LivingListpo1->TMMLx][LivingListpo1->TMMLy][LivingListpo1->TMMLz]=NOCELL;
		//		  LivingListpo1->TMMLx=a->PosX;
		//		  LivingListpo1->TMMLy=a->PosY;
		//		  LivingListpo1->TMMLz=a->PosZ;
		//		  ECMHigh[a->PosX][a->PosY][a->PosZ]=CTYPETMM;
		//		  ECMCan[a->PosX][a->PosY][a->PosZ]=CTYPETMM;
		//		  TGFb[a->PosX][a->PosY][a->PosZ]+=TMMTGFb;
		//		  IL6[a->PosX][a->PosY][a->PosZ]=IL6[a->PosX][a->PosY][a->PosZ]+TMMIL6;
		//		  //SDF1[a->PosX][a->PosY][a->PosZ]=SDF1[a->PosX][a->PosY][a->PosZ]+TMMSDF1;
		//		  LivingListpo1->TMMCycleTimer=LivingListpo1->TMMCycleTimer+1;
		//		  LivingListpo1++;
		//	  }
		//	  else
		//	  {
		//		  LivingListpo1->TMMCycleTimer=LivingListpo1->TMMCycleTimer+1;
		//		  TGFb[LivingListpo1->TMMLx][LivingListpo1->TMMLy][LivingListpo1->TMMLz]=TGFb[LivingListpo1->TMMLx][LivingListpo1->TMMLy][LivingListpo1->TMMLz]+TMMTGFb;
		//		  IL6[LivingListpo1->TMMLx][LivingListpo1->TMMLy][LivingListpo1->TMMLz]=IL6[LivingListpo1->TMMLx][LivingListpo1->TMMLy][LivingListpo1->TMMLz]+TMMIL6;
		//		  //SDF1[LivingListpo1->TMMLx][LivingListpo1->TMMLy][LivingListpo1->TMMLz]=SDF1[LivingListpo1->TMMLx][LivingListpo1->TMMLy][LivingListpo1->TMMLz]+TMMSDF1;
		//		  LivingListpo1++;
		//	  }
		//  }
		//  }
		//  }
	 // else
		//  LivingListpo1++;
	 // }

  //for (ApopListpo1=TMMApopList.begin();ApopListpo1!=TMMApopList.end();)// to search Apoplist
  //{
	 // int Timer1=ApopListpo1->TMMCurrentTimer+1;//every loop, the CurrentTimer is the same to the TimeStep
	 // if (Timer1==Tstep)// if the cell is inserted from the above process, then its MMCurrentTimer+1!=Tstep
		//                                     // therefore it will not be considered in this process
	 // {
		//   ApopListpo1->TMMApopTimer = ApopListpo1->TMMApopTimer+1;
		//   ApopListpo1->TMMCurrentTimer=ApopListpo1->TMMCurrentTimer+1;
		//  if (ApopListpo1->TMMApopTimer >=5)
		//  {
		//	  ApopListpo1->TMMApopStatus=APOPA;
		//	  ECMHigh[ApopListpo1->TMMLx][ApopListpo1->TMMLy][ApopListpo1->TMMLz]=NOCELL;
		//	  ECMCan[ApopListpo1->TMMLx][ApopListpo1->TMMLy][ApopListpo1->TMMLz]=NOCELL;
		//	  TMMAbsorbList.insert(TMMAbsorbList.end(),*ApopListpo1);
		//	  ApopListpo1=TMMApopList.erase(ApopListpo1);
		//  }
		//  else
		//	  ApopListpo1++;
	 // }
	 // else
		//  ApopListpo1++;
  //}

}
  
MicroEnv::~MicroEnv()
{
	delete []TGFb;
	delete []SDF1;
	delete []IL6;
	delete []STIFF;
	delete []ECMHigh;
	delete []ECMCan;
	delete []WNT2;
	delete []DKK2;
	//delete []A;
	delete []tmp;
	//delete []drug;
	//delete []codrug;
}
//-------------------------------------------------------------------------------------
void MicroEnv::PrintMIC(int Tstep, int drugindex, int coindex, string isM)
{
  string sFilename = "";
  sFilename = par::sOutDir + "/PTMIC/MIC_Thread_" + string(itostr(par::iT)) + "_Time_" + string(itostr(Tstep)) + ".txt";
  ofstream fout1;

  list <MIC>::iterator LivingListpo1;
  fout1.open(sFilename,ios::app);
  fout1<<"ID\tx\ty\tz\tstatus"<<endl;
  for (LivingListpo1=MICLivingList.begin();LivingListpo1!=MICLivingList.end();LivingListpo1++){
    fout1 << LivingListpo1->ID << "\t"
	  <<LivingListpo1->MICLx << "\t"
	  <<LivingListpo1->MICLy << "\t"
	  <<LivingListpo1->MICLz << "\t"
	  <<(LivingListpo1->MICCellStatus == 0?"1":"2")
	  <<endl;
  }
  for (LivingListpo1=MICApopList.begin();LivingListpo1!=MICApopList.end();LivingListpo1++) {
    fout1<< LivingListpo1->ID<<"\t"
	 <<LivingListpo1->MICLx<<"\t"
	 <<LivingListpo1->MICLy<<"\t"
	 <<LivingListpo1->MICLz<<"\t"
	 <<"3"
	 <<endl;
  }
  for (LivingListpo1=MICAbsorbList.begin();LivingListpo1!=MICAbsorbList.end();LivingListpo1++) {
    fout1<< LivingListpo1->ID<<"\t"
	 <<LivingListpo1->MICLx<<"\t"
	 <<LivingListpo1->MICLy<<"\t"
	 <<LivingListpo1->MICLz<<"\t"
	 <<"4"
	 <<endl;
  }
  fout1.close();
}

void MicroEnv::PrintPC(int Tstep, int drugindex, int coindex, string isM)
{
 /* string sFilename = "";
  sFilename = par::sOutDir + "/PC_" + isM + "_A_" + string(itostr(par::iA)) + "_B_" + string(itostr(par::iB)) + "_TP_" + string(itostr(Tstep)) + ".txt";
  ofstream fout1;

  list <PC>::iterator LivingListpo1;
  fout1.open(sFilename,ios::app);
  fout1<<"ID\tx\ty\tz\tstatus"<<endl;
  for (LivingListpo1=PCLivingList.begin();LivingListpo1!=PCLivingList.end();LivingListpo1++){
    fout1 << LivingListpo1->ID << "\t"
	  <<LivingListpo1->PCLx << "\t"
	  <<LivingListpo1->PCLy << "\t"
	  <<LivingListpo1->PCLz << "\t"
	  <<(LivingListpo1->PCCellStatus == 0?"1":"2")
	  <<endl;
  }
  for (LivingListpo1=PCApopList.begin();LivingListpo1!=PCApopList.end();LivingListpo1++) {
    fout1<< LivingListpo1->ID<<"\t"
	 <<LivingListpo1->PCLx<<"\t"
	 <<LivingListpo1->PCLy<<"\t"
	 <<LivingListpo1->PCLz<<"\t"
	 <<"3"
	 <<endl;
  }
  for (LivingListpo1=PCAbsorbList.begin();LivingListpo1!=PCAbsorbList.end();LivingListpo1++) {
    fout1<< LivingListpo1->ID<<"\t"
	 <<LivingListpo1->PCLx<<"\t"
	 <<LivingListpo1->PCLy<<"\t"
	 <<LivingListpo1->PCLz<<"\t"
	 <<"4"
	 <<endl;
  }
  fout1.close();*/
  
    // 	ofstream fout1;
    // 	char T[20],L[20]="LivingPC",B[20]="ApoptosisPC",A[20]="AbsorbPC",P[20]=".dat";
    // sprintf(T,"%d",Tstep);
    // strcat(L,T);strcat(B,T);strcat(A,T);
    // strcat(L,P);strcat(B,P);strcat(A,P);

    // 	fout1.open(L,ios::out);
    // 	//fout1<<"Living PC Cell list:"<<endl;
    // 	for (LivingListpo1=PCLivingList.begin();LivingListpo1!=PCLivingList.end();LivingListpo1++)
    // 	{
    // 		fout1<< LivingListpo1->ID<<" "<<LivingListpo1->PCLx<<" "<<LivingListpo1->PCLy<<" "<<LivingListpo1->PCLz<<endl;
    // 	}
    // 	fout1.close();
    // 	fout1.open(B,ios::out);
    // 	//fout1<<"Living MIC Cell list:"<<endl;
    // 	for (LivingListpo1=PCApopList.begin();LivingListpo1!=PCApopList.end();LivingListpo1++)
    // 	{
    // 		fout1<< LivingListpo1->ID<<" "<<LivingListpo1->PCLx<<" "<<LivingListpo1->PCLy<<" "<<LivingListpo1->PCLz<<endl;
    // 	}
    // 	fout1.close();
    // 	fout1.open(A,ios::out);

    // 	for (LivingListpo1=PCAbsorbList.begin();LivingListpo1!=PCAbsorbList.end();LivingListpo1++)
    // 	{
    // 		fout1<< LivingListpo1->ID<<" "<<LivingListpo1->PCLx<<" "<<LivingListpo1->PCLy<<" "<<LivingListpo1->PCLz<<endl;
    // 	}
    // 	fout1.close();
}

void MicroEnv::PrintMM(int Tstep, int drugindex, int coindex, string isM)
{
  string sFilename = "";
  sFilename = par::sOutDir + "/PTMM/MM_Thread_" + string(itostr(par::iT)) + "_Time_" + string(itostr(Tstep)) + ".txt";
  ofstream fout1;

  list <MM>::iterator LivingListpo1;
  fout1.open(sFilename,ios::app);
  fout1<<"ID\tx\ty\tz\tstatus"<<endl;
  for (LivingListpo1=MMLivingList.begin();LivingListpo1!=MMLivingList.end();LivingListpo1++){
    fout1 << LivingListpo1->ID << "\t"
	  <<LivingListpo1->MMLx << "\t"
	  <<LivingListpo1->MMLy << "\t"
	  <<LivingListpo1->MMLz << "\t"
	  <<(LivingListpo1->MMCellStatus == 0?"1":"2")
	  <<endl;
  }
  for (LivingListpo1=MMApopList.begin();LivingListpo1!=MMApopList.end();LivingListpo1++) {
    fout1<< LivingListpo1->ID<<"\t"
	 <<LivingListpo1->MMLx<<"\t"
	 <<LivingListpo1->MMLy<<"\t"
	 <<LivingListpo1->MMLz<<"\t"
	 <<"3"
	 <<endl;
  }
  for (LivingListpo1=MMAbsorbList.begin();LivingListpo1!=MMAbsorbList.end();LivingListpo1++) {
    fout1<< LivingListpo1->ID<<"\t"
	 <<LivingListpo1->MMLx<<"\t"
	 <<LivingListpo1->MMLy<<"\t"
	 <<LivingListpo1->MMLz<<"\t"
	 <<"4"
	 <<endl;
  }
  fout1.close();
    
}

void MicroEnv::PrintTMM(int Tstep, int drugindex, int coindex, string isM)
{
  /*string sFilename = "";
  sFilename = par::sOutDir + "/TMM_" + isM + "_A_" + string(itostr(par::iA))
    + "_B_" + string(itostr(par::iB))
    + "_TP_" + string(itostr(Tstep))
    + ".dat";
  ofstream fout1;

  list <TMM>::iterator LivingListpo1;
  fout1.open(sFilename,ios::out);
  fout1<<"ID\tx\ty\tz\tstatus"<<endl;
  for (LivingListpo1=TMMLivingList.begin();LivingListpo1!=TMMLivingList.end();LivingListpo1++){
    fout1 << LivingListpo1->ID << "\t"
	  <<LivingListpo1->TMMLx << "\t"
	  <<LivingListpo1->TMMLy << "\t"
	  <<LivingListpo1->TMMLz << "\t"
	  <<(LivingListpo1->TMMCellStatus == 0?"1":"2")
	  <<endl;
  }
  for (LivingListpo1=TMMApopList.begin();LivingListpo1!=TMMApopList.end();LivingListpo1++) {
    fout1<< LivingListpo1->ID<<"\t"
	 <<LivingListpo1->TMMLx<<"\t"
	 <<LivingListpo1->TMMLy<<"\t"
	 <<LivingListpo1->TMMLz<<"\t"
	 <<"3"
	 <<endl;
  }
  for (LivingListpo1=TMMAbsorbList.begin();LivingListpo1!=TMMAbsorbList.end();LivingListpo1++) {
    fout1<< LivingListpo1->ID<<"\t"
	 <<LivingListpo1->TMMLx<<"\t"
	 <<LivingListpo1->TMMLy<<"\t"
	 <<LivingListpo1->TMMLz<<"\t"
	 <<"4"
	 <<endl;
  }
  fout1.close();*/
    // 	list <TMM>::iterator LivingListpo1;
    // 	ofstream fout1;
    // 	    char T[20],L[20]="LivingTMM",B[20]="ApoptosisTMM",A[20]="AbsorbTMM",P[20]=".dat";
    // sprintf(T,"%d",Tstep);
    // strcat(L,T);strcat(B,T);strcat(A,T);
    // strcat(L,P);strcat(B,P);strcat(A,P);

    // 	fout1.open(L,ios::out);
    // 	//fout1<<"Living HH Cell list:"<<endl;
    // 	for (LivingListpo1=TMMLivingList.begin();LivingListpo1!=TMMLivingList.end();LivingListpo1++)
    // 	{
    // 		fout1<< LivingListpo1->ID<<" "<<LivingListpo1->TMMLx<<" "<<LivingListpo1->TMMLy<<" "<<LivingListpo1->TMMLz<<endl;
    // 	}
    // 	fout1.close();
    // 	fout1.open(B,ios::out);
    // 	//fout1<<"Living TMM Cell list:"<<endl;
    // 	for (LivingListpo1=TMMApopList.begin();LivingListpo1!=TMMApopList.end();LivingListpo1++)
    // 	{
    // 		fout1<< LivingListpo1->ID<<" "<<LivingListpo1->TMMLx<<" "<<LivingListpo1->TMMLy<<" "<<LivingListpo1->TMMLz<<endl;
    // 	}
    // 	fout1.close();
    // 	fout1.open(A,ios::out);

    // 	for (LivingListpo1=TMMAbsorbList.begin();LivingListpo1!=TMMAbsorbList.end();LivingListpo1++)
    // 	{
    // 		fout1<< LivingListpo1->ID<<" "<<LivingListpo1->TMMLx<<" "<<LivingListpo1->TMMLy<<" "<<LivingListpo1->TMMLz<<endl;
    // 	}
    // 	fout1.close();
}

void MicroEnv::PrintBMSC(int Tstep, int drugindex, int coindex)
{
	/*list <BMSC>::iterator LivingListpo1;
	ofstream fout1;
	char T[20],L[20]="BMSC",P[20]=".dat";
    sprintf_s(T,"%d",Tstep);
    strcat_s(L,T);//strcat(B,T);strcat(A,T);
    strcat_s(L,P);//strcat(B,P);strcat(A,P);

	fout1.open(L,ios::out);
	
	for (LivingListpo1=BMSCList.begin();LivingListpo1!=BMSCList.end();LivingListpo1++)
	{
		if ( STIFF[LivingListpo1->BMSCLx][LivingListpo1->BMSCLy][LivingListpo1->BMSCLz]>=400)
		{
			fout1<< LivingListpo1->ID<<" "<<LivingListpo1->BMSCLx<<" "<<LivingListpo1->BMSCLy<<" "<<LivingListpo1->BMSCLz<<endl;
		}

	}

	fout1.close();	*/
}

void MicroEnv::PrintSDF1(int Tstep, int drugindex, int coindex)
{
	//ofstream fout1;
	int i,j,k;
	//char T[20],L[20]="SDF1_",P[20]=".dat";
    //sprintf(T,"%d",Tstep);
    //strcat(L,T);
    //strcat(L,P);
    double sumsdf=0;
	//fout1.open(L,ios::out);
	//fout1<<"Living HH Cell list:"<<endl;
	for(i=0;i<GridWidth;i++)
	{
	  for(j=0;j<GridHeight;j++)
	  {
		  for(k=0;k<GridLength;k++)
		  {
			  //fout1<<i<<" "<<j<<" "<<k<<" "<<SDF1[i][j][k]<<endl;
			  sumsdf=sumsdf+SDF1[i][j][k];

		  }
	  }
	}
	//fout1.close();
	foutsdf<<Tstep<<" 	"<<sumsdf<<endl;

}

void MicroEnv::PrintSTIFF(int Tstep, int drugindex, int coindex, string isM)
{
  string sFilename = "";
  sFilename = par::sOutDir + "/STIFF_" + isM + "_A_" + string(itostr(par::iA))
    + "_B_" + string(itostr(par::iB))
    + "_TP_" + string(itostr(Tstep))
    + ".dat";
  ofstream fout1;
  int i,j,k;
  
  fout1.open(sFilename,ios::out);
  fout1<<"i\tj\tk\tstiff"<<endl;
  for(i=0;i<GridWidth;i++){
    for(j=0;j<GridHeight;j++) {
      for (k=0;k<GridLength;k++) {
	fout1<<i<<"\t"<<j<<"\t"<<k<<"\t"<<STIFF[i][j][k]<<endl;
      }
    }
  }
  fout1.close();
}

void MicroEnv::PrintTGFb(int Tstep, int drugindex, int coindex)
{
	/*ofstream fout1;
	int i,j,k;
	char T[20],L[20]="TGFb_",P[20]=".dat";
    sprintf_s(T,"%d",Tstep);
    strcat_s(L,T);
    strcat_s(L,P);

	fout1.open(L,ios::out);
	//fout1<<"Living HH Cell list:"<<endl;
	for(i=0;i<GridWidth;i++)
	{
	  for(j=0;j<GridHeight;j++)
	  {
		  for(k=0;k<GridLength;k++)
		  {
			  fout1<<i<<" "<<j<<" "<<k<<" "<<TGFb[i][j][k]<<endl;
		  }
	  }
	}
	fout1.close();*/


}

void MicroEnv::PrintECM(int Tstep, int drugindex, int coindex)
{
	/*ofstream fout1;
	int i,j,k;
	char T[20],L[20]="ECM_",P[20]=".dat";
    sprintf_s(T,"%d",Tstep);
    strcat_s(L,T);
    strcat_s(L,P);

	fout1.open(L,ios::out);
	//fout1<<"Living HH Cell list:"<<endl;
	for(i=0;i<GridWidth;i++)
	{
	  for(j=0;j<GridHeight;j++)
	  {
		  for(k=0;k<GridLength;k++)
		  {
			  fout1<<i<<" "<<j<<" "<<k<<" "<<ECMHigh[i][j][k]<<endl;
		  }
	  }
	  
	}
	fout1.close();*/

}
//-------------------------------------------------------------------------------------
 
double MicroEnv::Distance(int x1,int y1,int z1,int x2,int y2,int z2)
{
  double x,y,z,w;
  x = (double)(x1-x2);
  y = (double)(y1-y2);
  w = (double)(z1-z2);
  x = pow(x,2);
  y = pow(y,2);
  w = pow(w,2);
  z = sqrt(x+y+w);
  return z;
}

int MicroEnv::BoundaryCheck(int PosX,int PosY,int PosZ)
{
   int CheckFlag = 0;
   if ((PosX>GridWidth)||(PosX<0))
	   CheckFlag = 1;
   if((PosY>GridHeight)||(PosY<0))
	   CheckFlag = 1;
   if((PosZ>GridLength)||(PosZ<0))  //first must be in the cube
	   CheckFlag = 1;
   if ( ((PosX-20)*(PosX-20)+(PosZ-20)*(PosZ-20))>=140 )//( ((PosX-50)*(PosX-50)+(PosZ-50)*(PosZ-50))>=1000 ) // then x-y interception must be in the circle, 
	   CheckFlag = 1;        							 //	then ==> guarantee in the cylinder
   return CheckFlag;
}


void MicroEnv::PrintNumber(int Tstep, int drugindex, int coindex)
{
  //ofstream fout1;
  //int i,j,k;

  //fout4.open("CellInfo.txt",ios::out);
  // fout41.open("DeadCell.txt",ios::out);
  // fout4<<"T MIC PC   MM   TMM---Living"<<endl;
  // fout41<<"T MIC PC   MM   TMM---Dead"<<endl;
  // foutsdf.open("SdfInfo.txt",ios::out);
  // foutsdf<<"T	 SDF1total"<<endl;
  
  // fout4<<Tstep<<" "<<MICLivingList.size()<<"
  // "<<PCLivingList.size()<<" "<<MMLivingList.size()<<"
  // "<<TMMLivingList.size()<<" "<<endl;
  fout4 << Tstep
	<< "\t" << MICLivingList.size()
	<< "\t" << PCLivingList.size()
	<< "\t" << MMLivingList.size()
	<< "\t" << TMMLivingList.size()
	<< endl;
//	fout41<<Tstep<<" "<<MICAbsorbList.size()<<"
//	"<<PCAbsorbList.size()<<" "<<MMAbsorbList.size()<<"
//	"<<TMMAbsorbList.size()<<" "<<endl; fout1.close();
	
}

void MicroEnv::CloseFile()
{
	fout4.close();
	fout41.close();
	foutsdf.close();
	foutpop.close();
	foutsdf1.close();
	foutstiff.close();
	fouttgfb.close();
}

//--------------------------------------------------------------------------------------------------
double MicroEnv::Adh_MIC(double sti,double btz)  //calculation the adhesion rate of MIC to BMSC
{
	    /*double ADMIC;
	    mwArray mwsti(1,1,mxDOUBLE_CLASS);
	    mwArray mwbtz(1,1,mxDOUBLE_CLASS);
	    mwArray mwsuv(1,1,mxDOUBLE_CLASS);
	    mwArray mwad(1,1,mxDOUBLE_CLASS);
	 							 	  
	    if(sti<300)
		     sti=0.001+0.009*(sti-200)/100;
	    else if(sti<400)
		     sti=0.01+0.09*(sti-300)/100;
	    else if (sti>400)
		     sti=0.1+0.9*(sti-400)/130;
	    else
		     sti=1;
	  	  
	    mwbtz(1,1)=0;
	    mwsti(1,1)=sti;

	    Suv_Adh_MIC(2,mwsuv,mwad,mwsti,mwbtz);
	    ADMIC=mwad(1,1);  //adhesion rate

		return ADMIC;*/
	    return 1;
}

double MicroEnv::Adh_MIC2(double sti,double btz)
{
	double par[N_PAR]={0,0.5932, 0.2164, 0.1284, 0.3480, 0.8460, 0.0310, 0.9893, 0.5240, 0.8107, 0.8363, 0.6183, 0.3256, 0.0447, 0.1347, 0.4187, 0.4953, 0.2110, 0.6662, 0.1271, 0.6820, 0.7570, 0.3366, 0.5238, 0.1182, 0.5872, 0.3867, 0.5586, 0.3591, 0.8945, 0.3076, 0.8035,0.7326,0.1023,0};
	
	if (sti<300)
	    sti=0.001+0.009*(sti-200)/100;
	else if(sti<400)
		sti=0.01+0.09*(sti-300)/100;
	else if(sti>400)
		sti=0.1+0.9*(sti-400)/130;
	
	par[0]=sti;
	par[34]=btz;
	  
    double iniValue[N_VAR]={1,1,1,1,1,1,1,0.1749,0.6897};  //ok
    double iniTime = 0;
    int ntimepoints = 3;
    double timescale[3]={0,48,192};
	double *sol;
	int nsol = N_VAR * ntimepoints;
	int i,j;
	double rst;
	//--------------------------------------------------------
	sol = (double*)malloc(nsol * sizeof(double));
	if (!ODEmodel2(par, iniValue, iniTime, ntimepoints, timescale, sol))
	 {
	    printf("Something wrong with ODEs solution!\n");
	    return 0;
	 }
	   
	rst=sol[2*N_VAR-2]; //adhesion rate
	
		
	return rst;
	
}

Location* MicroEnv::MSearch(int Lx,int Ly,int Lz, double Di)  //MIC
{
  CpList.erase(CpList.begin(),CpList.end());
  int radius;
  double radius1;
  radius1=2*(1+4*Di*2)/10;
  radius=(int)radius1+1;  //radius=4;
  int ismoved;
  double ADMIC,ctmp;
  ADMIC=0;	
  if (ECMCan[Lx][Ly][Lz]!=BMSCCELL)  //current position of MIC is on BMSC 
   {
		ismoved=1;
   }
  else
   {
		//printf("%d  %d  %d  ",Lx,Ly,Lz);
		//printf("%f  %f \n",ECMHigh[Lx][Ly][Lz],ECMCan[Lx][Ly][Lz]);
	   //----------------------------------------------
        double mthresth=getRandom();	
	 
					
		//ADMIC=Adh_MIC(STIFF[Lx][Ly][Lz],0);
		ADMIC=Adh_MIC2(STIFF[Lx][Ly][Lz],0);
               // printf("Adhesion:%f    %f",ADMIC, Adh_MIC2(530,0));             	                	
  		if (mthresth>=ADMIC)
			ismoved=1;
		else
			ismoved=0;
              
   }
   //printf("   Adhesion : %d %d %d %f  %f  moved? %d\n",Lx,Ly,Lz,STIFF[Lx][Ly][Lz],ADMIC,ismoved);
  // getchar();

   if (ismoved==1)
   {
     for(int i=Lx-radius;i<=Lx+radius;i++)
       {
         for(int j=Ly-radius;j<=Ly+radius;j++)
          {
		     for(int k=Lz-radius;k<=Lz+radius;k++)
		      {
			     if(Distance(Lx, Ly,Lz, i, j,k)<=radius)
			  	  {				
					 if (BoundaryCheck(i,j,k)==0)						
				      {						
						//MIC should bind to BMSC, so that if ECMHigh=BMSC, MIC can move on BMSC
						if ((ECMHigh[i][j][k]==NOCELL)||(ECMHigh[i][j][k]==BMSCCELL)) // in migration search, the original point can not be candidate position
						{							
							if((i==Lx)&&(j==Ly)&&(k==Lz))
								continue;
							Location Loc(i,j,k);
							CpList.insert(CpList.end(),Loc);
							//printf("\n %d  %d  %d ",Loc.PosX,Loc.PosY,Loc.PosZ);
							
						}
					 }
			      }
              }
          }
       }
   }
 
   Location *LocReturn;
   LocReturn = null;
   if(CpList.empty()!=true)
   {
	   LocReturn = CompareLocation(Lx,Ly,Lz,Di);
   }
   else
   {       
	   LocReturn = null;	  
   }
   return LocReturn;
}

Location* MicroEnv::CompareLocation(int Lx, int Ly,int Lz,double Di)
{
   int counter,index;
   Location *TempLoc;
   double dr,ri,p,ctmp2,c,*r,v,sum=0;
   
   index=CpList.size();
   r=new double [index];
   index=0;
   TempLoc = null;
   list <Location> :: iterator LocCan;
   //printf("%d  %d  %d   %f  %f\n",Lx,Ly,Lz,ECMHigh[Lx][Ly][Lz],ECMCan[Lx][Ly][Lz]);
   //getchar();
   for(LocCan=CpList.begin();LocCan!=CpList.end();LocCan++)
    {
		ctmp2=1;
		counter=0;
		ri=Distance(LocCan->PosX,LocCan->PosY,LocCan->PosZ, Lx,Ly,Lz);
		p=1.0/(4*PI*Di)*exp(-pow(ri,2)/(4*PI*Di));

		if (ECMHigh[LocCan->PosX][LocCan->PosY][LocCan->PosZ]==BMSCCELL)
		{
			c=1.0;
			for(int x=-1;x<=1;x++)
		    {
		      for(int y=-1;y<=1;y++)
		      {
			    for(int z=-1;z<=1;z++)
			    {
			     if( (BoundaryCheck(LocCan->PosX+x,LocCan->PosY+y,LocCan->PosZ+z)==0)&& (Distance(LocCan->PosX+x,LocCan->PosY+y,LocCan->PosZ+z,LocCan->PosX,LocCan->PosY,LocCan->PosZ)<=1) )	
			      {
			        if (ECMHigh[LocCan->PosX+x][LocCan->PosY+y][LocCan->PosZ+z]!=NOCELL)
			        {
				     counter++;
					}
					if (ECMHigh[LocCan->PosX+x][LocCan->PosY+y][LocCan->PosZ+z]==CTYPECTL)
				    {
					   ctmp2=0.25;
				    }
				 }
				}
			  }
			}
		}
		else
		{
		  c=1.0/4;
		  for(int x=-1;x<=1;x++)
		 {
		 for(int y=-1;y<=1;y++)
		  {
			for(int z=-1;z<=1;z++)
			{
			 if( (BoundaryCheck(LocCan->PosX+x,LocCan->PosY+y,LocCan->PosZ+z)==0)&& (Distance(LocCan->PosX+x,LocCan->PosY+y,LocCan->PosZ+z,LocCan->PosX,LocCan->PosY,LocCan->PosZ)<=1) )	
			  {
			   if (ECMHigh[LocCan->PosX+x][LocCan->PosY+y][LocCan->PosZ+z]!=NOCELL)
			   {
				  counter++;				 
				  if (ECMCan[LocCan->PosX+x][LocCan->PosY+y][LocCan->PosZ+z]==BMSCCELL)
				  {					 
					 double e0,emax,Kfe,Fe0,Femax,edomina,ctmp;
					 e0=par::basicstiffness; 
					 emax=par::maxstiffness; 
					 Fe0=FE0; Femax=FEMAX;
					 Kfe=(emax-e0)/2;
					 edomina=STIFF[LocCan->PosX+x][LocCan->PosY+y][LocCan->PosZ+z];
					 edomina=(edomina-e0)/Kfe;
					 edomina=pow(edomina,2);
					 edomina=edomina/(1+edomina);
					 ctmp=edomina*Femax+Fe0;	 //  ctmp>0.5
					 
					 if (c<ctmp)
						 c=ctmp;
				  }
				  if (ECMHigh[LocCan->PosX+x][LocCan->PosY+y][LocCan->PosZ+z]==CTYPECTL)
				    {
					   ctmp2=0.25;
				    }
		       }
			  }
		  }
		 }
	     }
		} 
	   
	//printf("%d %d %d   %d  %f %f\n",LocCan->PosX,LocCan->PosY,LocCan->PosZ,counter,c,ctmp2);	
	  
	   v=1.0/16; // has no neighbors	
	   if (counter>=5&&counter<=6)
		   v=1.0/8;
	   if (counter>=3&&counter<=4)
		   v=1.0/4;
	   if (counter>=1&&counter<=2)
		   v=1.0;
	   
	   r[index++]=p*c*v*ctmp2;
	   
    }
   
   index=0;
   for (LocCan=CpList.begin();LocCan!=CpList.end();LocCan++)
	   sum=sum+r[index++];
   index=0;
   for (LocCan=CpList.begin();LocCan!=CpList.end();LocCan++)
	   r[index]=r[index++]/sum;

	sum=0;
   
	dr=getRandom();
	index=0;
	for (LocCan=CpList.begin();LocCan!=CpList.end();LocCan++)
	{
		sum=sum+r[index];
		index=index+1;
		if(dr<sum)
		{
			TempLoc=&(*(LocCan));
			break;
		}		
	}
	//printf("%d  %d  %d   %f  %f\n",TempLoc->PosX, TempLoc->PosY, TempLoc->PosZ,ECMHigh[TempLoc->PosX][TempLoc->PosY][TempLoc->PosZ],ECMCan[TempLoc->PosX][TempLoc->PosY][TempLoc->PosZ]);
   //  getchar();
	
	if (TempLoc != null)
		return TempLoc;
    else
	    return null; 
	   
}

//-----------------------------------------------------------------------------------------------------
Location* MicroEnv::MSearch_mm(int Lx,int Ly,int Lz, double Di)  //MM
{
  CpList.erase(CpList.begin(),CpList.end());
  int radius;
  double radius1;
  radius1=2*(1+4*Di*2)/10;
  radius=(int)radius1+1;
  
  for(int i=Lx-radius;i<=Lx+radius;i++)
   {
     for(int j=Ly-radius;j<=Ly+radius;j++)
      {
		  for(int k=Lz-radius;k<=Lz+radius;k++)
		  {
			  if(Distance(Lx, Ly,Lz, i, j,k)<=radius)			 
			  {
				 if (BoundaryCheck(i,j,k)==0)
				 {
					if ((ECMHigh[i][j][k]==NOCELL)||(ECMHigh[i][j][k]==BMSCCELL))
						{
							if((i==Lx)&&(j==Ly)&&(k==Lz))
								continue;
							Location Loc(i,j,k);							
							CpList.insert(CpList.end(),Loc);
							
						}
				 }
			 }
        }
      }
    }
   Location *LocReturn;
   LocReturn = null;
   if(CpList.empty()!=true)
    { 
		LocReturn = CompareLocation_migmm(Lx,Ly,Lz,Di);     
    }
   else
    {
        LocReturn = null;	
    }
   int i,j,k;
   i=LocReturn->PosX;
   j=LocReturn->PosY;
   k=LocReturn->PosZ;
   
   return LocReturn;
}  


Location* MicroEnv::CompareLocation_migmm(int Lx, int Ly,int Lz,double Di)
{
   int counter,index;
   Location *TempLoc;
   double dr,ri,p,c,*r,v,sum=0;
   double ctmp,ctmp2;
   index=CpList.size();
   r=new double [index];
   index=0;
   TempLoc = null;
   list <Location> :: iterator LocCan;
   //printf("%d  %d  %d   %f  %f\n",Lx,Ly,Lz,ECMHigh[Lx][Ly][Lz],ECMCan[Lx][Ly][Lz]);
   //getchar();
   for(LocCan=CpList.begin();LocCan!=CpList.end();LocCan++)
    {
		ctmp2=1;
		counter=0;
		ri=Distance(LocCan->PosX,LocCan->PosY,LocCan->PosZ, Lx,Ly,Lz);
		p=1.0/(4*PI*Di)*exp(-pow(ri,2)/(4*PI*Di));
		if (ECMHigh[LocCan->PosX][LocCan->PosY][LocCan->PosZ]==NOCELL)
		{
		   c=0.25;
		   for(int x=-1;x<=1;x++)
		   {
		     for(int y=-1;y<=1;y++)
		     {
			  for(int z=-1;z<=1;z++)
			  {
			   if( (BoundaryCheck(LocCan->PosX+x,LocCan->PosY+y,LocCan->PosZ+z)==0) && (Distance(LocCan->PosX+x,LocCan->PosY+y,LocCan->PosZ+z,LocCan->PosX,LocCan->PosY,LocCan->PosZ)<=1) )	
			   {
			     if ((ECMHigh[LocCan->PosX+x][LocCan->PosY+y][LocCan->PosZ+z]!=NOCELL) && (x*x+y*y+z*z!=0))
			      {
				    counter++;				 
				    if (ECMCan[LocCan->PosX+x][LocCan->PosY+y][LocCan->PosZ+z]==BMSCCELL)
				    {
					   ctmp=1.0;			
					   if (ctmp>c)
						   c=ctmp;
				    }
					if (ECMHigh[LocCan->PosX+x][LocCan->PosY+y][LocCan->PosZ+z]==CTYPECTL)
				    {
					   ctmp2=0.25;
				    }
		          }
			   }
			  }
		     }
	       }
		}
		else
		{
		  c=0.5;
		  for(int x=-1;x<=1;x++)
		  {
		     for(int y=-1;y<=1;y++)
		     {
			  for(int z=-1;z<=1;z++)
			  {
			   if( (BoundaryCheck(LocCan->PosX+x,LocCan->PosY+y,LocCan->PosZ+z)==0)&& (Distance(LocCan->PosX+x,LocCan->PosY+y,LocCan->PosZ+z,LocCan->PosX,LocCan->PosY,LocCan->PosZ)<=1) )	
			   {
			     if (ECMHigh[LocCan->PosX+x][LocCan->PosY+y][LocCan->PosZ+z]!=NOCELL)
			      {
				    counter++;				 
				    if (ECMCan[LocCan->PosX+x][LocCan->PosY+y][LocCan->PosZ+z]==BMSCCELL)
				    {
					   ctmp=1.0;
					   if (c<ctmp)
					     	 c=ctmp;
				    }
					if (ECMHigh[LocCan->PosX+x][LocCan->PosY+y][LocCan->PosZ+z]==CTYPECTL)
				    {
					   ctmp2=0.25;
				    }
		          }
			   }
			  }
		     }
	      }
		}

	   // printf("%d  %f  %f\n",counter,c,ctmp2);	  
	   	   v=1.0/16; // has no neighbors	
	   if (counter>=5&&counter<=6)
		   v=1.0/8;
	   if (counter>=3&&counter<=4)
		   v=1.0/4;
	   if (counter>=1&&counter<=2)
		   v=1.0;
	   
	   r[index++]=p*c*v*ctmp2;

    }
   index=0;
   for (LocCan=CpList.begin();LocCan!=CpList.end();LocCan++)
	   sum=sum+r[index++];
   index=0;
   for (LocCan=CpList.begin();LocCan!=CpList.end();LocCan++)
	   r[index]=r[index++]/sum;

	sum=0;
   
	dr=getRandom();
	index=0;
	for (LocCan=CpList.begin();LocCan!=CpList.end();LocCan++)
	{
		sum=sum+r[index];
		index=index+1;
		if(dr<sum)
		{
			TempLoc=&(*(LocCan));
			break;
		}		
	}

	//printf("%d  %d  %d   %f  %f\n",TempLoc->PosX, TempLoc->PosY, TempLoc->PosZ,ECMHigh[TempLoc->PosX][TempLoc->PosY][TempLoc->PosZ],ECMCan[TempLoc->PosX][TempLoc->PosY][TempLoc->PosZ]);
    //getchar();
	if (TempLoc != null)
		return TempLoc;
    else
	    return null; 
   
}

//------------------------------------------------------------------------------------------------------
Location* MicroEnv::MSearch_ctl(int Lx,int Ly,int Lz, double Di)  //CTL
{
  CpList.erase(CpList.begin(),CpList.end());
  int radius;
  double radius1;
  radius1=2*(1+4*Di*2)/10;
  radius=(int)radius1+1;
  
  for(int i=Lx-radius;i<=Lx+radius;i++)
   {
     for(int j=Ly-radius;j<=Ly+radius;j++)
      {
		  for(int k=Lz-radius;k<=Lz+radius;k++)
		  {
			  if(Distance(Lx, Ly,Lz, i, j,k)<=radius)			 
			  {
				 if (BoundaryCheck(i,j,k)==0)
				 {
					if ((ECMHigh[i][j][k]==NOCELL)||(((ECM3[i][j][k]==CTYPEMM)||(ECM3[i][j][k]==CTYPEMIC))&&(ECMCan[i][j][k]!=BMSCCELL)))
						{
							if((i==Lx)&&(j==Ly)&&(k==Lz))
								continue;
							Location Loc(i,j,k);							
							CpList.insert(CpList.end(),Loc);
							//printf("%d  %d  %d  %f\n",i,j,k,ECMHigh[i][j][k]);
							
						}
				 }
			 }
        }
      }
    }
   Location *LocReturn;
   LocReturn = null;
   if(CpList.empty()!=true)
    { 
		LocReturn = CompareLocation_migctl(Lx,Ly,Lz,Di);     
    }
   else
    {
        LocReturn = null;	
    }
   // getchar();
   //int i,j,k;
   //i=LocReturn->PosX;
   //j=LocReturn->PosY;
   //k=LocReturn->PosZ;
   ////printf("%d  %d  %d  %f",i,j,k,ECMHigh[i][j][k]);
   //getchar();
   return LocReturn;
}


Location* MicroEnv::CompareLocation_migctl(int Lx, int Ly,int Lz,double Di)  //CTL
{
   int counter,index;
   Location *TempLoc;
   double dr,ri,p,c,*r,v,sum=0;
   
   index=CpList.size();
   r=new double [index];
   index=0;
   TempLoc = null;
   list <Location> :: iterator LocCan;
   //printf("\n%d  %d  %d  %f\n",Lx,Ly,Lz,ECMHigh[Lx][Ly][Lz]);

   for(LocCan=CpList.begin();LocCan!=CpList.end();LocCan++)
    {
		counter=0;
		ri=Distance(LocCan->PosX,LocCan->PosY,LocCan->PosZ, Lx,Ly,Lz);
		p=1.0/(4*PI*Di)*exp(-pow(ri,2)/(4*PI*Di));		

		if ((ECMHigh[LocCan->PosX][LocCan->PosY][LocCan->PosZ]==CTYPEMIC) || (ECMHigh[LocCan->PosX][LocCan->PosY][LocCan->PosZ]==CTYPEMM))
			 c=1;		
		else		
		     c=0.75;
		//printf("%d  %d  %d  %f  ",LocCan->PosX,LocCan->PosY,LocCan->PosZ,ECMHigh[LocCan->PosX][LocCan->PosY][LocCan->PosZ]);
		//printf("old_c:%f  ",c);
		double tmp,tgfbval;
		tgfbval=TGFb[LocCan->PosX][LocCan->PosY][LocCan->PosZ];
		tmp=(tgfbval/CTLTGFB0)*(tgfbval/CTLTGFB0);
		c=c*(1/(1+tmp));  //prefer to migrate to the position with low TGFB concentration
		//printf("new_c:%f  ",c);
		for(int x=-1;x<=1;x++)
		{
		  for(int y=-1;y<=1;y++)
		   {
			for(int z=-1;z<=1;z++)
			 {
			  if( (BoundaryCheck(LocCan->PosX+x,LocCan->PosY+y,LocCan->PosZ+z)==0)&& (Distance(LocCan->PosX+x,LocCan->PosY+y,LocCan->PosZ+z,LocCan->PosX,LocCan->PosY,LocCan->PosZ)<=1) )	
			  {
			   if ((ECMHigh[LocCan->PosX+x][LocCan->PosY+y][LocCan->PosZ+z]!=NOCELL)&& (x*x+y*y+z*z !=0) && (((ECMHigh[LocCan->PosX+x][LocCan->PosY+y][LocCan->PosZ+z]==CTYPEMIC)||(ECMHigh[LocCan->PosX+x][LocCan->PosY+y][LocCan->PosZ+z]==CTYPEMM))))
			   {
				 counter++;	
		       }
			  }
			 }
		    }
	     }
		
	   //printf("counter: %d  ",counter)	;  
	   
	   v=1.0/16; // has no neighbors	
	   if (counter>=5&&counter<=6)
		   v=1/8;
	   if (counter>=3&&counter<=4)
		   v=1.0/4;
	   if (counter>=1&&counter<=2)
		   v=1;
	   
	   r[index++]=p*c*v;
	   //printf("  %f\n",p*c*v);

    }
   index=0;
   for (LocCan=CpList.begin();LocCan!=CpList.end();LocCan++)
	   sum=sum+r[index++];
   index=0;
   for (LocCan=CpList.begin();LocCan!=CpList.end();LocCan++)
	   r[index]=r[index++]/sum;

	sum=0;
  
	dr=getRandom();
	index=0;
	for (LocCan=CpList.begin();LocCan!=CpList.end();LocCan++)
	{
		sum=sum+r[index];
		index=index+1;
		if(dr<sum)
		{
			TempLoc=&(*(LocCan));
			break;
		}
		
	}
	//printf("\n%d  %d  %d  %f\n",TempLoc->PosX,TempLoc->PosY,TempLoc->PosZ,ECMHigh[TempLoc->PosX][TempLoc->PosY][TempLoc->PosZ]);
	//getchar();
	if (TempLoc != null)
		return TempLoc;
    else
	    return null; 
   
}


Location* MicroEnv::ProSearch_ctl(int Lx,int Ly,int Lz, double rad)
{
  CpList.erase(CpList.begin(),CpList.end());
  int radius;
  radius=(int)rad;
  

  for(int i=Lx-radius;i<=Lx+radius;i++)
    {
     for(int j=Ly-radius;j<=Ly+radius;j++)
      {
		  for (int k=Lz-radius;k<=Lz+radius;k++)
		  {
			  if(Distance(Lx, Ly, Lz, i, j,k)<=radius)
				  {
					  if (BoundaryCheck(i,j,k)==0)
						  {
							  if (ECMHigh[i][j][k]==NOCELL)
								  {
									  if((i==Lx)&&(j==Ly)&&(k==Lz))
										   continue;
									  Location Loc(i,j,k);
									  CpList.insert(CpList.end(),Loc);
									  //printf("%d  %d  %d %f\n",Loc.PosX,Loc.PosY,Loc.PosZ,ECMHigh[Loc.PosX][Loc.PosY][Loc.PosZ]);
									  //---------------------------------------									  
								   }
						  }
			       }
          }
      }
    }
   
   Location *LocReturn;
   LocReturn = null;
   if(CpList.empty()!=true)
     LocReturn = CompareLocation_proctl(Lx,Ly,Lz);
   else
     LocReturn = null;
   
   return LocReturn;
}


Location* MicroEnv::CompareLocation_proctl(int Lx, int Ly,int Lz)
{
   int counter,index = 0;
   Location *TempLoc;
   double dr,ri,p,c,v,sum,*r,ctmp;
   sum=0;
   r=new double [6];
   for (int i=0;i<6;i++)  //6 immediate neighbours	
	   r[i]=0;
  
   TempLoc = null;
   list <Location> :: iterator LocCan;
   for(LocCan=CpList.begin();LocCan!=CpList.end();LocCan++)
    {
		c=1.0/4;counter=0;
		ri=Distance(LocCan->PosX,LocCan->PosY,LocCan->PosZ, Lx,Ly,Lz);
		p=1.0;
	   for(int x=-1;x<=1;x++)
		{
		 for(int y=-1;y<=1;y++)
		  {
			  for (int z=-1;z<=1;z++)
			  {
				  if ( (BoundaryCheck(LocCan->PosX+x,LocCan->PosY+y,LocCan->PosZ+z)==0) && (Distance(LocCan->PosX+x,LocCan->PosY+y,LocCan->PosZ+z,LocCan->PosX,LocCan->PosY,LocCan->PosZ)<=1)   )
					  {
						  if (ECMHigh[LocCan->PosX+x][LocCan->PosY+y][LocCan->PosZ+z]!=NOCELL)
							  {
								  counter++;
								  if ((ECMHigh[LocCan->PosX+x][LocCan->PosY+y][LocCan->PosZ+z]==CTYPEMIC) || (ECMHigh[LocCan->PosX+x][LocCan->PosY+y][LocCan->PosZ+z]==CTYPEMM))
								  {
									  ctmp=0.5;
									  if(ctmp>c)
										  c=ctmp;
									 
								  }
								
							   }
				       }
			  }
		  }
	    }

	   //if (counter>=7||counter<1)
		   v=1.0/16;
	   if (counter>=5&&counter<=6)
		   v=1.0/8;
	   if (counter>=3&&counter<=4)
		   v=1.0/4;
	   if (counter>=1&&counter<=2)
		   v=1.0;
	   r[index]=p*c*v;
	   index++;	   	  
    }
   index=0;
    for(LocCan=CpList.begin();LocCan!=CpList.end();LocCan++)
	   sum=sum+r[index++];
	index=0;
   for (LocCan=CpList.begin();LocCan!=CpList.end();LocCan++)
	   r[index]=r[index++]/sum;
   sum=0;
   
    dr=getRandom();
	index=0;
	for (LocCan=CpList.begin();LocCan!=CpList.end();LocCan++)
	{
		sum=r[index]+sum;
		index=index+1;
		if(dr<sum)
		{
			TempLoc=&(*(LocCan));
			break;
		}
		
	}
	
	if (TempLoc != null)
		return TempLoc;
    else
	    return null; 
   
}

//---------------------------------------------------------------------------------------------------------------
//MIC/MM proliferation
Location* MicroEnv::ProSearch(int Lx,int Ly,int Lz, double rad)   //MIC/MM proliferation
{
  CpList.erase(CpList.begin(),CpList.end());
  int radius;
  radius=(int)rad;
  //radius++;
  //if (ECMCan[Lx][Ly][Lz]!=BMSCCELL)
//	ECMCan[Lx][Ly][Lz]=NOCELL;

  for(int i=Lx-radius;i<=Lx+radius;i++)
    {
     for(int j=Ly-radius;j<=Ly+radius;j++)
      {
		  for (int k=Lz-radius;k<=Lz+radius;k++)
		  {
			  if(Distance(Lx, Ly, Lz, i, j,k)<=radius)
				  {
					  if (BoundaryCheck(i,j,k)==0)
						  {
							  if (ECMHigh[i][j][k]==NOCELL||ECMHigh[i][j][k]==BMSCCELL)
								  {
									  if((i==Lx)&&(j==Ly)&&(k==Lz))
										   continue;
									  Location Loc(i,j,k);
									  CpList.insert(CpList.end(),Loc);
									  //---------------------------------------
									  
								   }
						  }
			       }
          }
      }
    }
   
   Location *LocReturn;
   LocReturn = null;
   if(CpList.empty()!=true)
     LocReturn = CompareLocation1(Lx,Ly,Lz);
   else
     LocReturn = null;
   return LocReturn;
}

//MIC/MM proliferation
Location* MicroEnv::CompareLocation1(int Lx, int Ly,int Lz)  //MIC/MM proliferation
{
   int counter,index = 0;
   Location *TempLoc;
   double dr,ri,p,c,v,sum,*r;
   sum=0;
   r=new double [6];
   for (int i=0;i<6;i++)   //6 immediate neighbours
	   r[i]=0;
   //printf("%d  %d  %d   %f\n",Lx,Ly,Lz,ECMHigh[Lx][Ly][Lz]);
   TempLoc = null;
   list <Location> :: iterator LocCan;
   for(LocCan=CpList.begin();LocCan!=CpList.end();LocCan++)
    {
		c=1.0/4;counter=0;
		ri=Distance(LocCan->PosX,LocCan->PosY,LocCan->PosZ, Lx,Ly,Lz);
		p=1.0;
		if (ECMHigh[LocCan->PosX][LocCan->PosY][LocCan->PosZ]==NOCELL)  //candidate free position is NOCELL
		{
		   for(int x=-1;x<=1;x++)
		   {
		       for(int y=-1;y<=1;y++)
		       {
			      for (int z=-1;z<=1;z++)
			      {
				         if ( (BoundaryCheck(LocCan->PosX+x,LocCan->PosY+y,LocCan->PosZ+z)==0) && (Distance(LocCan->PosX+x,LocCan->PosY+y,LocCan->PosZ+z,LocCan->PosX,LocCan->PosY,LocCan->PosZ)<=1)   )
					       {
						     if ((ECMHigh[LocCan->PosX+x][LocCan->PosY+y][LocCan->PosZ+z]!=NOCELL) && (x*x+y*y+z*z!=0))
							  {
								  counter++;
							   }
				            }
			      }
		       }
	       }
		}
		else
		{
			c=0.5;
		    for(int x=-1;x<=1;x++)
		    {
		        for(int y=-1;y<=1;y++)
		         {
			         for (int z=-1;z<=1;z++)
			          {
				           if ( (BoundaryCheck(LocCan->PosX+x,LocCan->PosY+y,LocCan->PosZ+z)==0) && (Distance(LocCan->PosX+x,LocCan->PosY+y,LocCan->PosZ+z,LocCan->PosX,LocCan->PosY,LocCan->PosZ)<=1)   )
					        {
						      if ((ECMHigh[LocCan->PosX+x][LocCan->PosY+y][LocCan->PosZ+z]!=NOCELL) && (x*x+y*y+z*z!=0))
							   {
								  counter++;
								  if (ECMCan[LocCan->PosX+x][LocCan->PosY+y][LocCan->PosZ+z]==BMSCCELL )
								  {
									  double e0,emax,Kfe,Fe0,Femax,edomina,ctmp;
									  e0=par::basicstiffness;
									  emax=par::maxstiffness;
									 									  
									  Fe0=FE0; Femax=FEMAX;
									  Kfe=(emax-e0)/2;
									  edomina=STIFF[LocCan->PosX+x][LocCan->PosY+y][LocCan->PosZ+z];
									  edomina=(edomina-e0)/Kfe;
									  edomina=pow(edomina,2);
									  edomina=edomina/(1+edomina);
									  ctmp=edomina*Femax+Fe0;
									  if(ctmp>c)
										  c=ctmp;
									 
								  }
								
							   }
				           }
			         }
		         }
	        }
		}

	   
	   
		   v=1.0/16;
	   if (counter>=5&&counter<=6)
		   v=1.0/8;
	   if (counter>=3&&counter<=4)
		   v=1.0/4;
	   if (counter>=1&&counter<=2)
		   v=1.0;
	   r[index]=p*c*v;
	   index++;
	   
	   //printf("%d  %d  %d  %d  %f  %f  %f  %f\n",LocCan->PosX,LocCan->PosY,LocCan->PosZ, counter,ECMHigh[LocCan->PosX][LocCan->PosY][LocCan->PosZ], p,c,v);
    }

   index=0;
    for(LocCan=CpList.begin();LocCan!=CpList.end();LocCan++)
	   sum=sum+r[index++];
	index=0;
   for (LocCan=CpList.begin();LocCan!=CpList.end();LocCan++)
	   r[index]=r[index++]/sum;
   sum=0;
    
    dr=getRandom();
	index=0;
	for (LocCan=CpList.begin();LocCan!=CpList.end();LocCan++)
	{
		sum=r[index]+sum;
		index=index+1;
		if(dr<sum)
		{
			TempLoc=&(*(LocCan));
			break;
		}
		
	}
	//printf("%d  %d  %d   %f\n",TempLoc->PosX,TempLoc->PosY,TempLoc->PosZ,ECMHigh[TempLoc->PosX][TempLoc->PosY][TempLoc->PosZ]);
	//getchar();
	if (TempLoc != null)
		return TempLoc;
    else
	    return null; 
   
}

//------------------------------------------------------------------------------

Location* MicroEnv::ProSearch_treg(int Lx,int Ly,int Lz, double rad)  //Treg
{
  CpList.erase(CpList.begin(),CpList.end());  
  int radius;
  radius=(int)rad;
 
  //printf("original: %d  %d  %d  %f\n",Lx,Ly,Lz,ECMHigh[Lx][Ly][Lz]);

  for(int i=Lx-radius;i<=Lx+radius;i++)
    {
     for(int j=Ly-radius;j<=Ly+radius;j++)
      {
		  for (int k=Lz-radius;k<=Lz+radius;k++)
		  {
			  if(Distance(Lx, Ly, Lz, i, j,k)<=radius)
				  {
					  if (BoundaryCheck(i,j,k)==0)
						  {
							  if (ECMHigh[i][j][k]==NOCELL)
								  {
									  if((i==Lx)&&(j==Ly)&&(k==Lz))
										   continue;
									  Location Loc(i,j,k);
									  CpList.insert(CpList.end(),Loc);
									  
									  //---------------------------------------
									  
								   }
						  }
			       }
          }
      }
    }
   //getchar();
   Location *LocReturn;
   LocReturn = null;
   if(CpList.empty()!=true)
     LocReturn = CompareLocation_protreg(Lx,Ly,Lz);
   else
     LocReturn = null;
   return LocReturn;
}

Location* MicroEnv::CompareLocation_protreg(int Lx, int Ly,int Lz)
{
   int counter,index = 0;
   Location *TempLoc;
   double dr,ri,p,c,ctmp,v,sum,*r;
   sum=0;
   r=new double [6];
   for (int i=0;i<6;i++)   //6 immediate neighbours
	   r[i]=0;
   //printf("%d  %d  %d  %f\n",Lx,Ly,Lz,ECMHigh[Lx][Ly][Lz]);
   TempLoc = null;
   list <Location> :: iterator LocCan;
   for(LocCan=CpList.begin();LocCan!=CpList.end();LocCan++)
    {
		//printf("%d  %d  %d  %f",LocCan->PosX,LocCan->PosY,LocCan->PosZ,ECMHigh[LocCan->PosX][LocCan->PosY][LocCan->PosZ]);
		c=1.0/4;counter=0;
		ri=Distance(LocCan->PosX,LocCan->PosY,LocCan->PosZ, Lx,Ly,Lz);
		p=1.0;
	   for(int x=-1;x<=1;x++)
		{
		 for(int y=-1;y<=1;y++)
		  {
			  for (int z=-1;z<=1;z++)
			  {
				  if ( (BoundaryCheck(LocCan->PosX+x,LocCan->PosY+y,LocCan->PosZ+z)==0) && (Distance(LocCan->PosX+x,LocCan->PosY+y,LocCan->PosZ+z,LocCan->PosX,LocCan->PosY,LocCan->PosZ)<=1)   )
					  {
						  if (ECMHigh[LocCan->PosX+x][LocCan->PosY+y][LocCan->PosZ+z]!=NOCELL)
							{
								 counter++;
								 if (ECMHigh[LocCan->PosX+x][LocCan->PosY+y][LocCan->PosZ+z]==CTYPECTL )
								 {									  
									ctmp=1;
									if(ctmp>c)
										c=ctmp;									 
								 }								
							 }
				       }
			  }
		  }
	    }
	   
	   //printf("  %d  ",counter);

		   v=1.0/16;
	   if (counter>=5&&counter<=6)
		   v=1.0/8;
	   if (counter>=3&&counter<=4)
		   v=1.0/4;
	   if (counter>=1&&counter<=2)
		   v=1.0;
	   r[index]=p*c*v;
	   //printf(" %f\n ",p*c*v);
	   index++;
	   

    }
   index=0;
    for(LocCan=CpList.begin();LocCan!=CpList.end();LocCan++)
	   sum=sum+r[index++];
	index=0;
   for (LocCan=CpList.begin();LocCan!=CpList.end();LocCan++)
	   r[index]=r[index++]/sum;
   sum=0;
   
    dr=getRandom();
	index=0;
	for (LocCan=CpList.begin();LocCan!=CpList.end();LocCan++)
	{
		sum=r[index]+sum;
		index=index+1;
		if(dr<sum)
		{
			TempLoc=&(*(LocCan));
			break;
		}		
	}

	//printf("%d  %d  %d  %f\n",TempLoc->PosX,TempLoc->PosY,TempLoc->PosZ,ECMHigh[TempLoc->PosX][TempLoc->PosY][TempLoc->PosZ]);
	//getchar();
	if (TempLoc != null)
		return TempLoc;
    else
	    return null;    
}


Location* MicroEnv::MSearch_treg(int Lx,int Ly,int Lz, double Di)  //Treg
{
  CpList.erase(CpList.begin(),CpList.end());
  int radius;
  double radius1;
  radius1=2*(1+4*Di*2)/10;
  radius=(int)radius1+1;
  //printf("original: %d  %d  %d  %f\n",Lx,Ly,Lz,ECMHigh[Lx][Ly][Lz]);
  for(int i=Lx-radius;i<=Lx+radius;i++)
   {
     for(int j=Ly-radius;j<=Ly+radius;j++)
      {
		  for(int k=Lz-radius;k<=Lz+radius;k++)
		  {
			  if(Distance(Lx, Ly,Lz, i, j,k)<=radius)			 
			  {
				 if (BoundaryCheck(i,j,k)==0)
				 {
					if ((ECMHigh[i][j][k]==NOCELL)||(ECMHigh[i][j][k]==CTYPECTL))
						{
							if((i==Lx)&&(j==Ly)&&(k==Lz))
								continue;
							Location Loc(i,j,k);							
							CpList.insert(CpList.end(),Loc);
							//printf("%d  %d  %d  %f\n",i,j,k,ECMHigh[i][j][k]);
							
						}
				 }
			 }
        }
      }
    }
   Location *LocReturn;
   LocReturn = null;
   
   if(CpList.empty()!=true)
    { 
		LocReturn = CompareLocation_migtreg(Lx,Ly,Lz,Di);     
    }
   else
    {
        LocReturn = null;	
    }
   
   /*int i,j,k;
   i=LocReturn->PosX;
   j=LocReturn->PosY;
   k=LocReturn->PosZ;
   printf("%d  %d  %d  %f\n",i,j,k,ECMHigh[i][j][k]);
   getchar();*/
   return LocReturn;
}

Location* MicroEnv::CompareLocation_migtreg(int Lx, int Ly,int Lz,double Di)
{
   int counter,index;
   Location *TempLoc;
   double dr,ri,p,c,c2,*r,sum=0;
   double v;
   
   index=CpList.size();
   r=new double [index];
   index=0;
   TempLoc = null;
   list <Location> :: iterator LocCan;
   //printf("selection: %d  %d  %d  %f\n",Lx,Ly,Lz,ECMHigh[Lx][Ly][Lz]);
   for(LocCan=CpList.begin();LocCan!=CpList.end();LocCan++)
    {
		//printf("candidate: %d  %d  %d  %f   ",LocCan->PosX,LocCan->PosY,LocCan->PosZ,ECMHigh[LocCan->PosX][LocCan->PosY][LocCan->PosZ]);
		c=0.25;counter=0;
		c2=0.5;
		ri=Distance(LocCan->PosX,LocCan->PosY,LocCan->PosZ, Lx,Ly,Lz);
		p=1.0/(4*PI*Di)*exp(-pow(ri,2)/(4*PI*Di));
		
		if (ECMHigh[LocCan->PosX][LocCan->PosY][LocCan->PosZ]==NOCELL)  //prefer to select the CD8+ candidate 
		    c=0.25;		
		else		
		    c=1;
		
	    for(int x=-1;x<=1;x++)
		{
		 for(int y=-1;y<=1;y++)
		  {
			for(int z=-1;z<=1;z++)
			{
			 if( (BoundaryCheck(LocCan->PosX+x,LocCan->PosY+y,LocCan->PosZ+z)==0)&& (Distance(LocCan->PosX+x,LocCan->PosY+y,LocCan->PosZ+z,LocCan->PosX,LocCan->PosY,LocCan->PosZ)<=1) )	
			 {
			 if ((ECMHigh[LocCan->PosX+x][LocCan->PosY+y][LocCan->PosZ+z]!=NOCELL) && (x*x+y*y+z*z!=0))
			 {
				 counter++;
				 
				 if (ECMHigh[LocCan->PosX+x][LocCan->PosY+y][LocCan->PosZ+z]==CTYPECTL)  //prefer to select the candidate who has CD8+ neighbour
				 {
					 double ctmp=1;
					 if (c2<ctmp)
						 c2=ctmp;
				 }

		     }
			 }
			}
		 }
	   }
	  
	   //printf("p:%f %d ",p,counter);
	   v=1.0/16; 	
	   
	   if (counter>=5&&counter<=6)
		   v=1.0/8;
	   if (counter>=3&&counter<=4)
		   v=1.0/2;
	   if (counter>=1&&counter<=2)
		   v=1.0;
	   //printf("%f    %f    %f\n",c,c2,p*c*v*c2);
	   r[index++]=p*c*v*c2;

    }
   //getchar();
   index=0;
   for (LocCan=CpList.begin();LocCan!=CpList.end();LocCan++)
	   sum=sum+r[index++];
   index=0;
   for (LocCan=CpList.begin();LocCan!=CpList.end();LocCan++)
	   r[index]=r[index++]/sum;

	sum=0;
   
	dr=getRandom();
	index=0;
	for (LocCan=CpList.begin();LocCan!=CpList.end();LocCan++)
	{
		sum=sum+r[index];
		index=index+1;
		if(dr<sum)
		{
			TempLoc=&(*(LocCan));
			break;
		}
		
	}
	//printf("%d  %d  %d  %f\n",TempLoc->PosX,TempLoc->PosY,TempLoc->PosZ,ECMHigh[TempLoc->PosX][TempLoc->PosY][TempLoc->PosZ]);
	//getchar();
	if (TempLoc != null)
		return TempLoc;
    else
	    return null; 
   
}

//----------------------------------------------------
//----------------------------------------------------------------

void MicroEnv::OgleBoundary( int TStep)
{

	/*char T[20],R[20]="r",B[20]="b",G[20]="g",O[20]="o",P[20]="cuon.dat";

        sprintf_s(T,"%d",TStep);
        strcat_s(R,T);strcat_s(B,T);strcat_s(G,T);strcat_s(O,T);
        strcat_s(R,P);strcat_s(B,P);strcat_s(G,P);strcat_s(O,P);

  FILE* outfileRed = fopen(R, "w");
  FILE* outfileGreen = fopen(G, "w");
  FILE* outfileBlue = fopen(B, "w");
  FILE* outfileOpac = fopen(O, "w");

  list <MIC>::iterator LivingListpo1;
  list <PC>::iterator LivingListpo2;
  list <MM>::iterator LivingListpo3;
  list <TMM>::iterator LivingListpo4;

  //
  //for (int i = 0; i < GridHeight; i++)
  // {
  // 	for (int j = 0; j < GridLength; j++)
  //   {
  //    for (int k = 0; k < GridWidth; k++)
  //     {
  //     	MatrixGrid[i][j][k] = 0;
  //     }
  //   }
  // }

  //for(LivingListpo1=MICLivingList.begin();LivingListpo1!=MICLivingList.end();LivingListpo1++)
  // {
	 //  MatrixGrid[LivingListpo1->MICLx][LivingListpo1->MICLy][LivingListpo1->MICLy] = CTYPEMIC;
  // }
	 
   
   
   
  for (int i = 0; i < GridWidth; i++)
   {
   	for (int j = 0; j < GridHeight; j++)
     {
      for (int k = 0; k < GridLength; k++)
       {
		  if ( BoundaryCheck(i,j,k)==1 )
		  {
       		  fputc((unsigned char) 200, outfileGreen);
    		  fputc((unsigned char) 200, outfileRed);
			  fputc((unsigned char) 200, outfileBlue);
			  fputc((unsigned char) 0, outfileOpac);
         }


        if (ECMHigh[i][j][k] ==CTYPEMIC) // MIC : red color
         {
       	  fputc((unsigned char) 0, outfileGreen);
    	  fputc((unsigned char) 0, outfileRed);
          fputc((unsigned char) 0, outfileBlue);
          fputc((unsigned char) 0, outfileOpac);
         }

        if(ECMHigh[i][j][k] ==CTYPEPC) // PC: blue
         {
          fputc((unsigned char) 0, outfileRed);
          fputc((unsigned char) 0, outfileGreen);
          fputc((unsigned char) 0, outfileBlue);
          fputc((unsigned char) 0, outfileOpac);
         }

       if (ECMHigh[i][j][k] ==CTYPEMM) // MM: green;
         {
       	  fputc((unsigned char) 0, outfileGreen);
    	  fputc((unsigned char) 0, outfileRed);
          fputc((unsigned char) 0, outfileBlue);
          fputc((unsigned char) 0, outfileOpac);
         }

        if (ECMHigh[i][j][k] ==CTYPETMM) // TMM: gray;
         {
       	  fputc((unsigned char) 0, outfileGreen);
    	  fputc((unsigned char) 0, outfileRed);
          fputc((unsigned char) 0, outfileBlue);
          fputc((unsigned char) 0, outfileOpac);
         }

        if (ECMHigh[i][j][k] ==BMSCCELL)// BMSC: yellow
         {
       	  fputc((unsigned char) 0, outfileGreen);
    	  fputc((unsigned char) 0, outfileRed);
          fputc((unsigned char) 0, outfileBlue);
          fputc((unsigned char) 0, outfileOpac);
         }

         if (ECMHigh[i][j][k] ==NOCELL) // empty;
         {
       	  fputc((unsigned char) 0, outfileGreen);
    	  fputc((unsigned char) 0, outfileRed);
          fputc((unsigned char) 0, outfileBlue);
          fputc((unsigned char) 0, outfileOpac);
         }

         // Reassign value
       }
     }
   }
  fputs("MIN 0 MAX 255\n", outfileRed);
  fputs("MIN 0 MAX 255\n", outfileGreen);
  fputs("MIN 0 MAX 255\n", outfileBlue);
  fputs("MIN 0 MAX 255\n", outfileOpac);

  fclose(outfileRed);
  fclose(outfileGreen);
  fclose(outfileBlue);
  fclose(outfileOpac);*/

}

void MicroEnv::OglePic(int TStep)
{

	/*char T[20],R[20]="r",B[20]="b",G[20]="g",O[20]="o",P[20]="cuon.dat";

        sprintf(T,"%d",TStep);
        strcat(R,T);strcat(B,T);strcat(G,T);strcat(O,T);
        strcat(R,P);strcat(B,P);strcat(G,P);strcat(O,P);

  FILE* outfileRed = fopen(R, "w");
  FILE* outfileGreen = fopen(G, "w");
  FILE* outfileBlue = fopen(B, "w");
  FILE* outfileOpac = fopen(O, "w");

   
   
   
  for (int i = 0; i < GridWidth; i++)
   {
   	for (int j = 0; j < GridHeight; j++)
     {
      for (int k = 0; k < GridLength; k++)
       {

        if (ECMHigh[i][j][k] ==CTYPEMIC) // MIC : red color
         {
       	  fputc((unsigned char) 0, outfileGreen);
    	  fputc((unsigned char) 255, outfileRed);
          fputc((unsigned char) 0, outfileBlue);
          fputc((unsigned char) 0, outfileOpac);
         }

        if(ECMHigh[i][j][k] ==CTYPEPC) // PC: red
         {
          fputc((unsigned char) 0, outfileRed);
          fputc((unsigned char) 0, outfileGreen);
          fputc((unsigned char) 255, outfileBlue);
          fputc((unsigned char) 0, outfileOpac);
         }

       if (ECMHigh[i][j][k] ==CTYPEMM) // MM: green;
         {
       	  fputc((unsigned char) 255, outfileGreen);
    	  fputc((unsigned char) 0, outfileRed);
          fputc((unsigned char) 0, outfileBlue);
          fputc((unsigned char) 0, outfileOpac);
         }

        if (ECMHigh[i][j][k] ==CTYPETMM) // TMM: gray;
         {
       	  fputc((unsigned char) 160, outfileGreen);
    	  fputc((unsigned char) 160, outfileRed);
          fputc((unsigned char) 160, outfileBlue);
          fputc((unsigned char) 0, outfileOpac);
         }

        if (ECMHigh[i][j][k] ==BMSCCELL)// BMSC: yellow
         {
			 if(j>50)//( k<50)
			 {
       		  fputc((unsigned char) 0, outfileGreen);
    		  fputc((unsigned char) 0, outfileRed);
			  fputc((unsigned char) 0, outfileBlue);
			  fputc((unsigned char) 0, outfileOpac);
			 }
			 else
			 {
       		  fputc((unsigned char) 0, outfileGreen);
    		  fputc((unsigned char) 0, outfileRed);
			  fputc((unsigned char) 0, outfileBlue);
			  fputc((unsigned char) 0, outfileOpac);
			 }

         }

         if (ECMHigh[i][j][k] ==NOCELL) // empty;
         {
       	  fputc((unsigned char) 0, outfileGreen);
    	  fputc((unsigned char) 0, outfileRed);
          fputc((unsigned char) 0, outfileBlue);
          fputc((unsigned char) 0, outfileOpac);
         }

         // Reassign value
       }
     }
   }
  fputs("MIN 0 MAX 255\n", outfileRed);
  fputs("MIN 0 MAX 255\n", outfileGreen);
  fputs("MIN 0 MAX 255\n", outfileBlue);
  fputs("MIN 0 MAX 255\n", outfileOpac);

  fclose(outfileRed);
  fclose(outfileGreen);
  fclose(outfileBlue);
  fclose(outfileOpac);*/

}

void MicroEnv::OglePic1(int TStep)
{

	/*char T[20],R[20]="r",B[20]="b",G[20]="g",O[20]="o",P[20]="cuon.dat";

        sprintf(T,"%d",TStep);
        strcat(R,T);strcat(B,T);strcat(G,T);strcat(O,T);
        strcat(R,P);strcat(B,P);strcat(G,P);strcat(O,P);

  FILE* outfileRed = fopen(R, "w");
  FILE* outfileGreen = fopen(G, "w");
  FILE* outfileBlue = fopen(B, "w");
  FILE* outfileOpac = fopen(O, "w");

	list <MIC>::iterator LivingListpo1;
	list <BMSC>::iterator LivingListpo2;


	  for (LivingListpo1=MICLivingList.begin();LivingListpo1!=MICLivingList.end();)
	  {
		  //int TempTimer=LivingListpo1->MICCurrentTimer+1;// every loop, the CurrentTimer is one less than the TimeStep	  
		  //if (TempTimer == TStep) // to guarantee that every cell is handled only once in each timestep.
		  {
			  //LivingListpo1->MICCurrentTimer=LivingListpo1->MICCurrentTimer+1;// Here just plus 1 inside the loop for considering about the new generated cells that added in each timestep
			  if (LivingListpo1->MICCellStatus!=StaN) //enter the cell cycle
				  MatrixGrid[LivingListpo1->MICLx][LivingListpo1->MICLy][LivingListpo1->MICLz]=2;
			  else
				  MatrixGrid[LivingListpo1->MICLx][LivingListpo1->MICLy][LivingListpo1->MICLz]=1;
			  LivingListpo1++;
		  }
		  //else
			 // LivingListpo1++;
	  }

	  for (LivingListpo2=BMSCList.begin();LivingListpo2!=BMSCList.end();)
	  {
		  //int TempTimer1=LivingListpo2->BMSCCurrentTimer+1;// every loop, the CurrentTimer is one less than the TimeStep	  
		  //if (TempTimer1 == TStep) // to guarantee that every cell is handled only once in each timestep.
		  {
			  //LivingListpo2->BMSCCurrentTimer=LivingListpo2->BMSCCurrentTimer+1;// Here just plus 1 inside the loop for considering about the new generated cells that added in each timestep
			  if (STIFF[LivingListpo2->BMSCLx][LivingListpo2->BMSCLy][LivingListpo2->BMSCLz]>=110) // BMSC is activated
				  MatrixGrid[LivingListpo2->BMSCLx][LivingListpo2->BMSCLy][LivingListpo2->BMSCLz]=3;
			  else
				  MatrixGrid[LivingListpo2->BMSCLx][LivingListpo2->BMSCLy][LivingListpo2->BMSCLz]=4;
			  LivingListpo2++;
		  }
		  else
			  LivingListpo2++;
	  }

   
   
   
  for (int i = 0; i < GridWidth; i++)
   {
   	for (int j = 0; j < GridHeight; j++)
     {
      for (int k = 0; k < GridLength; k++)
       {

        if (ECMHigh[i][j][k] ==CTYPEMIC) // MIC : red color
         {
			 if (MatrixGrid[i][j][k] ==2)//enter cell cycle MIC: light red
			 {
       			  fputc((unsigned char) 0, outfileGreen);
    			  fputc((unsigned char) 255, outfileRed);
				  fputc((unsigned char) 0, outfileBlue);
				  fputc((unsigned char) 0, outfileOpac);
			 }
			 else 
			 {
				 if (MatrixGrid[i][j][k] ==1)// not cell cycle MIC: gray
			 {
       			  fputc((unsigned char) 160, outfileGreen);
    			  fputc((unsigned char) 160, outfileRed);
				  fputc((unsigned char) 160, outfileBlue);
				  fputc((unsigned char) 0, outfileOpac);
			 }
				 else
				 {
					fputc((unsigned char) 0, outfileGreen);
    			  fputc((unsigned char) 0, outfileRed);
				  fputc((unsigned char) 0, outfileBlue);
				  fputc((unsigned char) 0, outfileOpac);

				 }
			 }	
         }

        if(ECMHigh[i][j][k] ==CTYPEPC) // PC: empty
         {
          fputc((unsigned char) 0, outfileRed);
          fputc((unsigned char) 0, outfileGreen);
          fputc((unsigned char) 0, outfileBlue);
          fputc((unsigned char) 0, outfileOpac);
         }

       if (ECMHigh[i][j][k] ==CTYPEMM) // MM: empty;
         {
       	  fputc((unsigned char) 0, outfileGreen);
    	  fputc((unsigned char) 0, outfileRed);
          fputc((unsigned char) 0, outfileBlue);
          fputc((unsigned char) 0, outfileOpac);
         }

        if (ECMHigh[i][j][k] ==CTYPETMM) // TMM: empty;
         {
       	  fputc((unsigned char) 0, outfileGreen);
    	  fputc((unsigned char) 0, outfileRed);
          fputc((unsigned char) 0, outfileBlue);
          fputc((unsigned char) 0, outfileOpac);
         }

        if (ECMHigh[i][j][k] ==BMSCCELL)// BMSC: green or grey
         {
			 if (MatrixGrid[i][j][k] ==3)// BMSC is activated: yellow
			 {
       			  fputc((unsigned char) 255, outfileGreen);
    			  fputc((unsigned char) 255, outfileRed);
				  fputc((unsigned char) 0, outfileBlue);
				  fputc((unsigned char) 0, outfileOpac);
			 }
			 if (MatrixGrid[i][j][k] ==4)// BMSC is not activated: empty
			 {
       			  fputc((unsigned char) 0, outfileGreen);
    			  fputc((unsigned char) 0, outfileRed);
				  fputc((unsigned char) 0, outfileBlue);
				  fputc((unsigned char) 0, outfileOpac);
			 }
         }

         if (ECMHigh[i][j][k] == NOCELL) // empty;
         {
       	  fputc((unsigned char) 0, outfileGreen);
    	  fputc((unsigned char) 0, outfileRed);
          fputc((unsigned char) 0, outfileBlue);
          fputc((unsigned char) 0, outfileOpac);
         }

         // Reassign value
       }
     }
   }


   
  fputs("MIN 0 MAX 255\n", outfileRed);
  fputs("MIN 0 MAX 255\n", outfileGreen);
  fputs("MIN 0 MAX 255\n", outfileBlue);
  fputs("MIN 0 MAX 255\n", outfileOpac);

  fclose(outfileRed);
  fclose(outfileGreen);
  fclose(outfileBlue);
  fclose(outfileOpac);*/


}
