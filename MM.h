#ifndef MM_H_
#define MM_H_
#define APOPA 2  // IF Absorb, MMAPOPA
#define APOPY 1  // IF Apoptosis, MMAPOPY
#define APOPN 0  // IF not Apoptosis, MMAPOPN
#define MMMSTATY 1 // M status Yes
#define MMMSTATN 0 // M status No
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
#define MAXMM 1 // Cell cycle in M stage


// FOR EACH KIND OF FOUR CELLS

#define T2 3 // PROLIFERATION WAITING TIME OF MM, SET DEFAULT VALUE 3
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
//#include ".\mic.h"
//#include "MIC.h"
////#include ".\pc.h"
//#include "PC.h"
//#include <list>
////#include ".\location.h"
#include "Location.h"

class MM
{
 public:
  int MMApopTimer; // the timer that records the apoptosis time steps of the MM cell
  int MMWaitProTimer; // the timer that records the waiting time of proliferation
  int MMCycleTimer; // the timer that records the time of Cell cycle
  int MMCurrentTimer; // the timer that records the current time of the cell
  int MMLx,MMLy,MMLz;// The location of it

  int MMApopStatus;// apoptosis status: Absorbed ->MMAPOPA;  Being apoptosis ->MMAPOPY; No->MMAPOPN
  int MMCellStatus; //The four status of cell cycle: G0/G1 ->StaG01;  S ->StaS;  G2 ->StaG2;  M ->StaM
					// if the value is StaN-> did not enter cell cycle
  int MMFreeSpace; // Whether there are free spaces: Yes->MMFREESPAY  No->MMFREESPAN
  //int ProWhat; // determine THE CELL TYPES THAT generate
  int MMGenNo;
  MM (int,int,int,int);
  virtual ~MM();


  int ID;  // The ID number of the MM
  int CellType; // The type of the cell.
  double secretSDF1,secretTGFb,secretIL6,secretStiff; // the SDF1, TGF_b and IL6 that the MM cell secret as well as the Stiffness and so on of the pathway.

};

#endif /*MM_H_*/
