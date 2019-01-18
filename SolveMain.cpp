//#include ".\stdafx.h"
//#include "StdAfx.h"
#include <stdio.h>
#include <stdlib.h>
#include <list>
#include "MIC.h"
#include "MM.h"
#include "PC.h"
#include "TMM.h"
#include "MicroEnv.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <fstream>
#include <string.h>
#include <iostream>
#include <getopt.h>

#include "par.hpp"
#include "Utilities.hpp"

using namespace std;


void get_opts(int argc, char *argv[]){
  int opt;
  while ((opt = getopt( argc, argv, "s:t:I:D:C:T:f:o:P:i:e:b:l:n:M:?h")) != -1) {
    switch (opt) {
    case 's':
      par::iT=atoi(optarg);  //Thread ID
      break;
    case 't':
      par::iMaxTime = atoi(optarg);
      // printf ("I got parameter t as: %d.\n",par::iMaxTime);
      break;
    case 'I':
      par::sw_mic = atoi(optarg);  //MIC Switch
      break;
    case 'D':
      par::sw_mm = atoi(optarg);   //MM Switch
      break;
    case 'C':
      par::sw_cd8 = atoi(optarg);  //CD8+ Switch
      break;
    case 'T':
      par::sw_treg = atoi(optarg); //Treg Switch 
      break;
    case 'f':
      par::sParFile = optarg;
      // printf ("I got parameter f as: %s.\n",par::sParFile.c_str());
      break;
    case 'o':
      par::sOutDir = optarg;
      // printf ("I got paramter o as: %s.\n",par::par::sOutDir.c_str());
      break;
    case 'P':
      par::bP = true;
      break;
    case 'i':
      par::iP_Start = atoi(optarg);
      break;
    case 'e':
      par::iP_End = atoi(optarg);
      break;
    case 'b': 
      par::iA = atoi(optarg);  //BTZ 
      break;
    case 'l':
      par::iB = atoi(optarg);  //LEN
      break;
    case 'n':
      par::iC = atoi(optarg);  //Third drug
      break;
    case 'M':
      if (("M" == par::sM) | ("N" == par::sM))
	        par::sM = optarg;
      else {
	          fprintf (stderr, "ERROR: %s is wrong BMSC type. \n",
		      par::sM.c_str());
	          exit(1);
           }
      break;
    case 'h':
      printf("Usage: %s\n", argv[0]);
      printf("\t[-s Thread ID]\n");
      printf("\t[-t timesteps]\n");
      printf("\t[-I The switch of MIC]; 1 for on, 0 for off. \n"); 
      printf("\t[-D The switch of MM]; 1 for on, 0 for off. \n");
      printf("\t[-C The switch of CD8+]; 1 for on, 0 for off. \n");
      printf("\t[-T The switch of Treg]; 1 for on, 0 for off. \n");
      printf("\t[-f parameter_file_name]\n");
      printf("\t[-o output_file_dir]\n");
      printf("\t[-P print details\n");
      printf("\t[-i timestep]; start time point for printing. \n");
      printf("\t[-e timestep]; end time point for printing. \n");
      printf("\t[-b level of BTZ, 0 to 11]; 0 means no drug. \n");
      printf("\t[-l level of LEN, 0 to 11]; 0 means no drug. \n");
      printf("\t[-n level of the third drug, 0 to 1]; 0 means no drug.\n");
      printf("\t[-M type of BMSC, M or N]\n");
      exit(1);
      break;
    default:
      abort();
    }
  }
}  


int main(int argc, char *argv[])
{
        get_opts(argc, argv);

       par::iA=par::iA/10;
       par::iB=par::iB/10;
       par::iC=par::iC/10;

       printf("thread: %d\n",par::iT);
       printf("time: %d\n",par::iMaxTime); 
       printf("MIC on or off? %d\n",par::sw_mic);
       printf("MM on or off? %d\n",par::sw_mm);
       printf("CD8 on or off? %d\n",par::sw_cd8);
       printf("Treg on or off? %d\n",par::sw_treg); 
       printf("BTZ: %f\n",par::iA);
       printf("LEN: %f\n",par::iB);
       printf("Third drug: %f\n",par::iC);
       printf("BMSC: %s\n",par::sM.c_str());
       //getchar();

       //srand(par::iT+1);

        int t_start = clock(); 

	ofstream foutpop1;
        string sFilename = "";
        sFilename = par::sOutDir + "/POP/T" + string(itostr(par::iT))+ "_allpop.txt";
        foutpop1.open(sFilename,ios::out);
	for(int kk=1;kk<=1;kk++)
	{
         int t =(int) time(NULL) ; 
         srand (t);
         int t_start = clock();
	
         if (par::sM=="N")
              par::inistiffness=200;
    
         int TimeStep;
         int counter1=0;
         int counter2=0;
         int counter3=0;
  
         list <MIC>::iterator Listpo;
         list <PC>::iterator Listpo1;
         list <MM>::iterator Listpo2;
         list <TMM>::iterator Listpo3;

         TimeStep=1;
         MicroEnv ME;
         ME.MicroEnvInit();

        
         ME.Microdrug(par::iB,par::iA,par::sM);
         ME.Diff(0,lamada_SDF1,ME.tmp,ME.SDF1);     //SDF1 diffustion
         ME.DiffTGFb(0,lamada_TGFb,ME.tmp,ME.TGFb); //TGFb diffustion 
         ME.DiffTrd(0,lamada_SDF1,ME.tmp,ME.SDF1);  //SDF1 degration by the third drug
		 
         // cout<<"+++++++++++++++++++"<<endl;
         // cout<< "Thread " << par::iT << " Time step: "<<TimeStep-1<<endl;
         // counter1=ME.BMSCList.size();
         // cout<<"THE BMSC cells: "<<counter1<<endl;

         // counter1=ME.MICApopList.size();
         // counter2=ME.MICAbsorbList.size();
         // counter3=ME.MICLivingList.size();

         // cout<<"Living MIC cells: "<<counter3<<endl;
         // cout<<"Being apoptosis MIC cells: "<<counter1<<endl;
         // cout<<"Dead MIC cells: "<<counter2<<endl;
	
         // counter1=ME.PCApopList.size();
         // counter2=ME.PCAbsorbList.size();
         // counter3=ME.PCLivingList.size();
	
         // cout<<"Living PC cells: "<<counter3<<endl;
         // cout<<"Being apoptosis PC cells: "<<counter1<<endl;
         // cout<<"Dead PC cells: "<<counter2<<endl;
	
         // counter1=ME.MMApopList.size();
         // counter2=ME.MMAbsorbList.size();
         // counter3=ME.MMLivingList.size();

         // cout<<"Living MM cells: "<<counter3<<endl;
         // cout<<"Being apoptosis MM cells: "<<counter1<<endl;
         // cout<<"Dead MM cells: "<<counter2<<endl;
	
         // counter1=ME.TMMApopList.size();
         // counter2=ME.TMMAbsorbList.size();
         // counter3=ME.TMMLivingList.size();
	
         // cout<<"Living TMM cells: "<<counter3<<endl;
         // cout<<"Being apoptosis TMM cells: "<<counter1<<endl;
         // cout<<"Dead TMM cells: "<<counter2<<endl;

         //ME.PrintSDF1(TimeStep);
         //ME.OglePic(TimeStep);

         
         while(TimeStep<=par::iMaxTime)
         {
            ME.MicroAddDrug(TimeStep,par::iA,par::iB, par::iC);
            ME.BMSCCellCheck(TimeStep);	  
	  //----------------------------------------
	    if (ME.sw_mic==1) 
	            ME.MICCellCheck(TimeStep,par::iB);
	        
	    if (ME.sw_mm==1)
		    ME.MMCellCheck(TimeStep,par::iB);
	  //----------------------------------------
	    if (ME.sw_cd8==1)  //switch to control CD8+
	            ME.CTLCellCheck(TimeStep,par::iB);

            if (ME.sw_treg==1) //switch to control CD8+ Treg
	            ME.TregCellCheck(TimeStep,par::iB);
			
          //----------------------------------------
	        ME.Diff(TimeStep,lamada_SDF1,ME.tmp,ME.SDF1);
	        ME.DiffTGFb(TimeStep,lamada_TGFb,ME.tmp,ME.TGFb);
			ME.DiffTrd(TimeStep, lamada_SDF1,ME.tmp,ME.SDF1);  //SDF1 degration by the third drug

	 // if (TimeStep%20==0)	
         //     getchar();
          
	        //----------------------------------------
            //ME.OglePic(TimeStep);	
            //ME.OgleBoundary(TimeStep);
            //ME.OglePic1(TimeStep);	
              if (0) 
	        {
	           cout<<"+++++++++++++++++++"<<endl;
	           //cout<<"Time step: "<<TimeStep<<endl;
	           counter1=ME.MICApopList.size();
	           counter2=ME.MICAbsorbList.size();
	           counter3=ME.MICLivingList.size();   
	           cout<<"Living MIC cells: "<<counter3<<endl;
	           cout<<"Being apoptosis MIC cells: "<<counter1<<endl;
	           cout<<"Dead MIC cells: "<<counter2<<endl;
		
	          // counter1=ME.PCApopList.size();
	          // counter2=ME.PCAbsorbList.size();
	          // counter3=ME.PCLivingList.size();		
	          // cout<<"Living PC cells: "<<counter3<<endl;
	          // cout<<"Being apoptosis PC cells: "<<counter1<<endl;
	          // cout<<"Dead PC cells: "<<counter2<<endl;
		
	           counter1=ME.MMApopList.size();
	           counter2=ME.MMAbsorbList.size();
	           counter3=ME.MMLivingList.size();
	           cout<<"Living MM cells: "<<counter3<<endl;
	           cout<<"Being apoptosis MM cells: "<<counter1<<endl;
	           cout<<"Dead MM cells: "<<counter2<<endl;
		
	          // counter1=ME.TMMApopList.size();
	          // counter2=ME.TMMAbsorbList.size();
	          // counter3=ME.TMMLivingList.size();
	          // cout<<"Living TMM cells: "<<counter3<<endl;
	          // cout<<"Being apoptosis TMM cells: "<<counter1<<endl;
	          // cout<<"Dead TMM cells: "<<counter2<<endl;
            }	//	*/
    
			

            printf("\n");
	        printf("\tt: %d \tMIC: \t%ld \tMM: \t%ld \tCD8+: \t%ld \tTreg: \t%ld\n", TimeStep, ME.MICLivingList.size(), ME.MMLivingList.size(), ME.CTLLivingList.size(), ME.TregLivingList.size());
	        printf("\n");
	  	  
	        foutpop1<<TimeStep<<" "<<ME.MICLivingList.size()<<" "<<ME.MMLivingList.size()<<" "<<ME.CTLLivingList.size()<<" "<<ME.TregLivingList.size()<<endl;	
	        foutpop1<<"\n";
      

                // if (TimeStep==100)
               //       getchar();

      // printf("Thread: \t%d, \tt: %d \tlive MIC: \t%ld, \tapop MIC: \t%ld, \tdead MIC: \t%ld; \n",
            // 	     par::iT, TimeStep-1,
            // 	     ME.MICLivingList.size(),ME.MICApopList.size(),ME.MICAbsorbList.size());
            // printf("Thread: \t%d, \tt: %d \tlive PC: \t%L, \tapop PC: \t%L, \tdead PC: \t%L; \n",
            // 	     par::iT, TimeStep-1, ME.PCLivingList.size(),ME.PCApopList.size(),ME.PCAbsorbList.size());
            // printf("Thread: \t%d, \tt: %d \tlive MM: \t%L, \tapop MM: \t%L, \tdead MM: \t%L; \n",
            // 	     par::iT, TimeStep-1, ME.MMLivingList.size(),ME.MMApopList.size(),ME.MMAbsorbList.size());
            // printf("Thread: \t%d, \tt: %d \tlive TMM: \t%L, \tapop TMM: \t%L, \tdead TMM: \t%L; \n",
            // 	     par::iT, TimeStep-1, ME.TMMLivingList.size(),ME.TMMApopList.size(),ME.TMMAbsorbList.size());
            // printf("Thread: \t%d, \tt: %d \tlive Tot: \t%L, \tapop Tot: \t%L, \tdead Tot: \t%L; \n",
            // 	     par::iT, TimeStep-1,
            // 	     ME.MICLivingList.size() + ME.PCLivingList.size() + ME.MMLivingList.size() + ME.TMMLivingList.size(),
            // 	     ME.MICApopList.size() + ME.PCApopList.size() + ME.MMApopList.size() + ME.TMMApopList.size(),
            // 	     ME.MICAbsorbList.size() + ME.PCAbsorbList.size() + ME.MMAbsorbList.size() + ME.TMMAbsorbList.size());
		
            //ME.PrintPC(TimeStep);
            //ME.PrintMIC(TimeStep);
            //ME.PrintBMSC(TimeStep);
            //ME.PrintMM(TimeStep);
            //ME.PrintTMM(TimeStep);
            /*ME.PrintNumber(TimeStep,par::iB,par::iA);
            if (par::bP && TimeStep >= par::iP_Start && TimeStep <= par::iP_End){
	             ME.PrintSTIFF(TimeStep,par::iB,par::iA,par::sM);
	             ME.PrintMIC(TimeStep,par::iB,par::iA,par::sM);
	             getchar();
	             ME.PrintPC(TimeStep,par::iB,par::iA,par::sM);
	             ME.PrintMM(TimeStep,par::iB,par::iA,par::sM);
	             ME.PrintTMM(TimeStep,par::iB,par::iA,par::sM);
            }*/
            //ME.PrintWnt(TimeStep);
	        //ME.PrintDKK(TimeStep);
            //ME.PrintTGFb(TimeStep);
            //ME.PrintSDF1(TimeStep);
            //ME.PrintECM(TimeStep);
         
			//----------------------------------------------------
		    //print a case to visualize the turmor growth		
			if ((TimeStep==50)||(TimeStep==300))
			{
				ME.PrintMIC(TimeStep, par::iA, par::iB,par::sM);
		                ME.PrintMM(TimeStep, par::iA, par::iB,par::sM);
			}
		
		    TimeStep++;
        }
       
      
		//-----------------------------------------------------
       ME.CloseFile();

	   
	}
      
	foutpop1.close();

    printf("Time used for %d steps is: %f Sec. \n", par::iMaxTime,(clock() - t_start)/(double)(CLOCKS_PER_SEC)); 
    //getchar();
    return 0;
}

