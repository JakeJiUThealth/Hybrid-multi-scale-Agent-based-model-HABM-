#ifndef MICROENV_H_
#define MICROENV_H_
#include <iostream>
#include <fstream>
#include <string.h>
#include <math.h> 
#include <stdio.h>
#include <list>
#include <iterator>
#include <stdlib.h>

#include "BMSC.h"
#include "MIC.h"
#include "MM.h"
#include "PC.h"
#include "TMM.h"
#include "CTL.h"
#include "Treg.h"

#include "Location.h"
#include "Utilities.hpp"



#define CYCLEG01BEGIN 1  //The stage of G0/G1 begin
#define CYCLEG01END 6 //The stage of G0/G1 end
#define CYCLESBEGIN 7 //The stage of S begin
#define CYCLESEND 10  //The stage of S end
#define CYCLEG2BEGIN 11   //The stage of G2 begin
#define CYCLEG2END 12   //The stage of G2 end
#define CYCLEMBEGIN 13  //The stage of M begin
#define CYCLEMEND 14  //The stage of M end

#define NOCELL 0  // the value of ECMHigh: NOCELL-> has no cell, CTYPEMIC-> HOLDS A MIC CELL ......
#define CTYPEMIC 1 // THE TYPE OF THE CELL IS MIC
#define CTYPETMM 4 // THE TYPE OF THE CELL IS TMM
#define CTYPEMM 3 // THE TYPE OF THE CELL IS MM
#define CTYPEPC 2 // THE TYPE OF THE CELL IS PC    ====> have been defined in each head file so here do not need to redefine
#define BMSCCELL 5   // For BMSC cube 
#define CTYPECTL 6 // The Type of the cell is CD8 CTL
#define CTYPETREG 7 // The type of the cell is CD8 Treg

//define the Grid scale
#define GridWidth 40//100
#define GridLength 40//100
#define GridHeight 40//100
#define Radius 6 // the initial seedings set radius
#define MaxStep  40

//-----------------------------------------------------------
#define Thal_deg 6
#define lamada_SDF1 0.612 // no-change
#define lamada_TGFb 0.44  // no-change
#define DEG 0.03  //degration rate of SDF1
#define DEGT 0.05 //degration rate of TGFB


#define MICp0 0.01;        // no change
#define MICP_pathway 0.058; // 
#define MICapop0 0.005;    // no change
#define MICapopdrug 0.03;  //0.042

#define myeloma_off_cd8 0.1;   //optimized! 0.55 to turn on; 0.1 to turn off
#define myeloma_on_cd8 0.6;   //0.5 
#define myeloma_btz1 0.4;  //For MIC. should be optimized?
#define myeloma_btz2 0.5;  //For MM . should be optimized?

#define MMp0 0.065;   //0.035      
#define MMP_pathway 0.0;  //no change
#define MMapop0 0.005;     //no change0.01
#define MMapopdrug 0.055;  //no change

#define MICSDF1 0.02;   //
#define	MMTGFb 0.25;    //
#define BMSCTGFb 0.25;  //BMSC secrite TGFb
//-------------------------------------------------
#define	radx 1;// the radius for proliferate 

//-The parameters for immune system-----------------
#define CTLapop0 0.01
#define CTLTGFBapop 0.01  //should be optimized
#define CTLproa 0.01  //should be optimized(induced by Treg and TGFB)
#define CTLprob 0.02  //optimized(induced by Len)0.02
#define Td_LEN  0.47  //optimized
#define Td_Treg 0.5  //optimized(0.55)0.05
#define Td_off_Treg 0.5 //
#define CTLTGFB0 0.5   //should be optimized

#define Tregapop0 0.01
#define TregTGFBapop 0.02  //0.02
#define Tregpro0 0.01
#define Tregprotgfb 0.023  //should be optimized
#define Treg_LEN 0.8;   //LEN can reduce the number of Treg cells
//-------------------------------------------------


// the related 
#define K_A 0.5;
#define K_TGFb 1.0e-5;
#define lamada_IL6 0.52
#define lamada_DKK2 0.44
#define lamada_WNT2 0.54
#define lamada_STIFF 0.43
#define MaxCell 1
#define null 0 
#define PI 3.14159265358979323846
#define c_dkk 1.0;
#define K_WNT 1.0e-3;
#define PCapop0 0.01;
#define TMMapop0 0.01;
#define PCapopdrug 0.06;
#define TMMapopdrug 0.06;
#define PCp0 0.02; 
#define PCP_pathway 0.0;
#define migD0 2;// 2um/timestep
#define migD_max 10; // 10um/timestep 5um/hr=1.388888888888889e-6 cm/sec 1grid/timestep
#define migK_D 0.1;
#define K_prol 0.0;  
#define MICIL6 0.25;// 
#define	MICTGFb 0.55;//
#define	MICSTIFF 0.15;//
#define PCSDF1 0.0025; //  
#define PCIL6 0.25;// 
#define	PCTGFb 0.55;//
#define	PCSTIFF 0.55;//
#define MMSDF1 0.0025; //  
#define MMIL6 0.25;// 
#define	MMSTIFF 0.55;//
#define TMMSDF1 0.0025; //  
#define TMMIL6 0.25;// 
#define	TMMTGFb 0.55;//
#define	TMMSTIFF 0.55;//
#define BMSCSDF1 0.0025; // 



using namespace std;

extern "C" int ODEmodel(double*, double*, double, int, double*, double*);  
extern "C" int ODEmodel2(double*,double*,double, int,double*, double*);

class MicroEnv
{
public:
    //double **TGFbHigh,**WntHigh, **ECMHigh;
	double ***TGFb,***SDF1, ***IL6,***STIFF, ***VF, ***MatrixGrid,***ECMHigh,***ECMCan,***ECM3,***BF, ***U1, ***U, ***WNT2, ***DKK2;
	double ***Protain; // for the interface of Jing's function
	double ***tmp;
	double drug,codrug,trddrug;

	int sw_mic,sw_mm,sw_cd8,sw_treg;

	ofstream fout4,foutsdf,fout41,foutpop;
	ofstream foutsdf1,foutstiff,fouttgfb;

	int i1; // the number that to generate the ID number of each cell.
    list <Location> CpList;
    list <MIC> MICLivingList;  // The list that contains the living MIC cells
    list <MIC> MICApopList;    // The list that contains the being apoptosis MIC cells
	list <MIC> MICAbsorbList;  // The list that contains the dead MIC cells

    list <TMM> TMMLivingList; // The list that contains the living TMM cells
    list <TMM> TMMApopList;   // The list that contains the being apoptosis TMM cells
	list <TMM> TMMAbsorbList; // The list that contains the dead TMM cells
    
	list <PC> PCLivingList;   // The list that contains the living PC cells
    list <PC> PCApopList;     // The list that contains the being apoptosis PC cells
	list <PC> PCAbsorbList;   // The list that contains the dead PC cells
    
	list <MM> MMLivingList;   // The list that contains the living MM cells
    list <MM> MMApopList;     // The list that contains the being apoptosis MM cells
	list <MM> MMAbsorbList;   // The list that contains the dead MM cells

	list <CTL> CTLLivingList; // The list that contains the living CTL cells
	list <CTL> CTLApopList;   // The list that contains the being apoptosis CTL cells
	list <CTL> CTLAbsorbList;  // The list that contains the dead CTL cells

	list <Treg> TregLivingList; // The list that contains the living CTL cells
	list <Treg> TregApopList;   // The list that contains the being apoptosis CTL cells
	list <Treg> TregAbsorbList;  // The list that contains the dead CTL cells

	list <BMSC> BMSCList;// The list that contains the BMSC cells who are ALWAYS there

    list <MIC>:: iterator SaveCellPosition1; //  BACKUP LISTS
    list <TMM>:: iterator SaveCellPosition2;
	list <PC>:: iterator SaveCellPosition3;
    list <MM>:: iterator SaveCellPosition4;
    
	void MicroEnvInit();
	void Microdrug(int drugindex, int coindex, string sM);
	void MicroAddDrug(int Tstep, double drugindex, double coindex,double thirdidx);
	int BoundaryCheck(int,int,int); // THE BOUNDARY IS A CYLINDER
	
	void MICCellCheck(int Tstep, int drugindex);
	void TMMCellCheck(int Tstep,int drugindex);
	void MMCellCheck(int Tstep,int drugindex);
	void PCCellCheck(int Tstep,int drugindex);
	
	void CTLCellCheck(int Tstep,int drugindex);
	void TregCellCheck(int Tstep,int drugindex);

	void BMSCCellCheck(int Tstep);

	double Calculation_Stiff(int, double);  //matlab based ode solver
	double Sur_MIC(int,double,double);	    //matlab based ode solver
	double Adh_MIC(double,double);          //matlab based ode solver

	double Calculation_Stiff2(int, double); //fortan based ode solver
	double Sur_MIC2(int,double,double);	//fortran based ode solver
	double Adh_MIC2(double,double);     //fortran based ode solver

	void PrintMIC(int Tstep, int drugindex, int coindex, string isM);
	void PrintTMM(int Tstep, int drugindex, int coindex, string isM);
	void PrintMM(int Tstep, int drugindex, int coindex, string isM);
	void PrintPC(int Tstep, int drugindex, int coindex, string isM);
	void PrintBMSC(int Tstep, int drugindex, int coindex); 
	void PrintNumber(int Tstep, int drugindex, int coindex);
	void CloseFile();
	void PrintTGFb(int Tstep, int drugindex, int coindex);// ++
	void PrintSDF1(int Tstep, int drugindex, int coindex);// ++
	void PrintIL6(int Tstep, int drugindex, int coindex);// ++
	void PrintSTIFF(int Tstep, int drugindex, int coindex, string isM);// ++
	void PrintECM(int Tstep, int drugindex, int coindex);
	//void Diff(double, double, double, double **); // DIFFUSION FUNCTION
	void Diff(int, double, double ***, double ***);
	void DiffTGFb(int, double, double ***, double ***);
	void DiffTrd(int, double, double ***, double ***);

	int JudgeStatus(int CycleTimer); // JUDGE CELLCYCLE STATUS OF EACH CELL
    void OglePic(int);
	void OglePic1(int);
	void OgleBoundary(int);
	double Distance(int,int,int,int,int,int); // calculate distance in 3-D space
	double getMigD(int,int,int); // get migration radius
	double getRandom(); // random generator
	//double diffusionGauss(double,double,double,double,double);
	double diffusionGauss(int,int,int);

	Location* ProSearch(int,int,int, double); // proliferation search function	
	Location* CompareLocation1(int, int,int); 	

	Location* MSearch(int, int,int,double); // migration search function  MIC
	Location* CompareLocation(int, int,int, double);	 	
	
	Location* MSearch_mm(int, int,int,double); // migration search function  MM
	Location* CompareLocation_migmm(int, int,int,double); 

    Location* ProSearch_ctl(int,int,int, double); // proliferation search function for CTL 
	Location* CompareLocation_proctl(int, int,int);  //for CTL
	Location* MSearch_ctl(int, int,int,double);  //migratrion search function CTL
	Location* CompareLocation_migctl(int, int,int,double);

	Location* MSearch_treg(int, int,int,double);  //migratrion search function Treg
	Location* CompareLocation_migtreg(int, int,int,double);  //for treg
	Location* ProSearch_treg(int,int,int, double); // proliferation search function for Treg
	Location* CompareLocation_protreg(int, int,int);  //for Treg

	MicroEnv();
	virtual ~MicroEnv();

};

#endif /*MICROENV_H_*/
