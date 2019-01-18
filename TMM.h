#ifndef TMM_H_
#define TMM_H_
#define APOPA 2  // IF Absorb, TMMAPOPA
#define APOPY 1  // IF Apoptosis, TMMAPOPY
#define APOPN 0  // IF not Apoptosis, TMMAPOPN
#define MSTATY 1 // M status Yes
#define MSTATN 0 // M status No
#define FREESPAY 1 // THERE ARE FREE SPACES. 
#define FREESPAN 0 // THERE ARE NO FREE SPACES.
#define CTYPEMIC 1 // THE TYPE OF THE CELL IS MIC
#define CTYPEPC 2 // THE TYPE OF THE CELL IS PC
#define CTYPEMM 3 // THE TYPE OF THE CELL IS MM
#define CTYPETMM 4 // THE TYPE OF THE CELL IS TMM
#define StaN 0 // Do not enter Cell cycle
#define StaG01 1 // Cell cycle in G0/G1 stage
#define StaS 2 // Cell cycle in S stage
#define StaG2 3 // Cell cycle in G2 stage
#define StaM 4 // Cell cycle in M stage


// FOR EACH KIND OF FOUR CELLS

// #define T4 3 // WAITING TIME OF TMM, SET DEFAULT VALUE 3
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

class TMM
{
 public:
  int TMMApopTimer; // the timer that records the apoptosis time steps of the TMM cell
  int TMMWaitProTimer; // the timer that records the time of waiting for proliferation
  int TMMCycleTimer; // the timer that records the time of Cell cycle
  int TMMCurrentTimer; // the timer that records the current time of the cell
  int TMMLx,TMMLy,TMMLz;// The location of it

  int TMMApopStatus;// apoptosis status: Absorbed ->TMMAPOPA;  Being apoptosis ->TMMAPOPY; No->TMMAPOPN
  int TMMCellStatus; //The four status of cell cycle: G0/G1 ->StaG01;  S ->StaS;  G2 ->StaG2;  M ->StaM
					// if the value is StaN-> did not enter cell cycle
  double TMMProSpeed; /* speed of proliferate. */
  int TMMFreeSpace; // Whether there are free spaces: Yes->TMMFREESPAY  No->TMMFREESPAN
  TMM (int,int,int,int);
  virtual ~TMM();


  int ID;  // The ID number of the TMM
  int CellType; // The type of the cell.
  double secretSDF1,secretTGFb,secretIL6; // the SDF1, TGF_b and IL6 that the TMM cell secret.

};

#endif /*TMM_H_*/
