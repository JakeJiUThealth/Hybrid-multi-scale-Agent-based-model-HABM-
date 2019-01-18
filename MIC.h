#ifndef MIC_H_
#define MIC_H_


#define APOPA 2  // IF Absorb, APOPA
#define APOPY 1  // IF Apoptosis, APOPY
#define APOPN 0  // IF not Apoptosis, APOPN
#define MICMSTATY 1 // M status Yes
#define MICMSTATN 0 // M status No
#define FREESPAY 1 // THERE ARE FREE SPACES. 
#define FREESPAN 0 // THERE ARE NO FREE SPACES.

// TYPES FOR EACH KIND OF FOUR CELLS
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


#define MICPROMIC 1 // MIC -> 2MIC
#define MICPROPC 2 // MIC -> MIC+PC
#define MICCHANGEPC 3 // MIC -> 2PC

#define T1 3 // THE PROLIFERATION WAITING TIME OF MIC, SET DEFAULT VALUE 3
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

// The MIC cell class
class MIC
{
 public:
  int MICApopTimer; // the timer that records the apoptosis time steps of the MIC cell
  int MICWaitProTimer; // the timer that records the time of waiting for proliferation
  int MICCycleTimer; // the timer that records the time of Cell cycle
  int MICCurrentTimer; // the timer that records the current time of the cell
  int MICLx,MICLy,MICLz; // The location of each MIC cell

  int MICApopStatus;// apoptosis status: Absorbed ->MICAPOPA;  Being apoptosis ->MICAPOPY; No apoptosis ->MICAPOPN
  int MICCellStatus; // The four status of cell cycle: G0/G1 ->StaG01;  S ->StaS;  G2 ->StaG2;  M ->StaM
					// if the value is StaN-> did not enter cell cycle
  double MICProSpeed; /* speed of proliferate. */
  int MICFreeSpace; // Whether there are free spaces: Yes->MICFREESPAY  No->MICFREESPAN
  int ProWhat; // determine MIC generates MIC or PC. the value is MICPROPC MICPROMIC, MICCHANGEPC

  MIC (int,int,int,int);
  virtual ~MIC();

  int ID;  // The ID number of the MIC
  int CellType; // The type of the cell.
  double secretSDF1,secretTGFb,secretIL6,secretStiff; // the SDF1 that the MIC cell secret as well as the TGF_b, Stiffness and IL6 decided by the pathway.

};

#endif /*MIC_H_*/
