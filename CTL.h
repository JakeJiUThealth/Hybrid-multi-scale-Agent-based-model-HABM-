#ifndef CTL_H_
#define CTL_H_
#define APOPA 2  // IF Absorb, PCAPOPA
#define APOPY 1  // IF Apoptosis, PCAPOPY
#define APOPN 0  // IF not Apoptosis, PCAPOPN
#define PCMSTATY 1 // M status Yes
#define PCMSTATN 0 // M status No
#define FREESPAY 1 // THERE ARE FREE SPACES. 
#define FREESPAN 0 // THERE ARE NO FREE SPACES.
#define CTYPEMIC 1 // THE TYPE OF THE CELL IS MIC
#define CTYPEPC 2 // THE TYPE OF THE CELL IS PC
#define CTYPEMM 3 // THE TYPE OF THE CELL IS MM
#define CTYPETMM 4 // THE TYPE OF THE CELL IS TMM
#define CTYPECTL 6 // The Type of the cell is CD8 CTL
#define CTYPETREG 7 // The type of the cell is CD8 Treg
#define StaN 0 // Do not enter Cell cycle
#define StaG01 1 // Cell cycle in G0/G1 stage
#define StaS 2 // Cell cycle in S stage
#define StaG2 3 // Cell cycle in G2 stage
#define StaM 4 // Cell cycle in M stage
#define MAXCTL 3 // Cell cycle in M stage

// FOR EACH KIND OF FOUR CELLS

#define T3 3 // PROLIFERATION WAITING TIME OF PC, SET DEFAULT VALUE 3 WHICH CAN BE CHANGED
#include <iostream>
#include <fstream>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <list>
#include <iterator>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <time.h>

#include "Location.h"

// THE CTL CELL CLASS
class CTL
{
 public:
  int CTLApopTimer; // the timer that records the apoptosis time steps of the PC cell
  int CTLWaitProTimer; // the timer that records the waiting time of proliferation
  int CTLCycleTimer; // the timer that records the time of Cell cycle
  int CTLCurrentTimer; // the timer that records the current time of the cell
  int CTLLx,CTLLy,CTLLz;// The location of it

  int CTLApopStatus;// apoptosis status: Absorbed ->PCAPOPA;  Being apoptosis ->PCAPOPY; No->PCAPOPN
  int CTLCellStatus; //The four status of cell cycle: G0/G1 ->StaG01;  S ->StaS;  G2 ->StaG2;  M ->StaM
					// if the value is StaN-> did not enter cell cycle
  int CTLFreeSpace; // Whether there are free spaces: Yes->PCFREESPAY  No->PCFREESPAN
  int CTLGenNo;
  CTL (int,int,int,int);
  virtual ~CTL();


  int ID;  // The ID number of the PC
  int CellType; // The type of the cell.
  double secretSDF1,secretTGFb,secretIL6,secretStiff; // the SDF1, TGF_b and IL6 that the PC cell secret. also of the pathway and Stiffness FOR uE MODULE
};

#endif /*CLT_H_*/
