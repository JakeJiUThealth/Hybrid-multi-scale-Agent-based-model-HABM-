#ifndef BMSC_H_
#define BMSC_H_


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
//#include ".\mm.h"
//#include "MM.h"
////#include ".\pc.h"
//#include "PC.h"
////#include ".\tmm.h"
//#include "TMM.h"
//#include <list>
////#include ".\location.h"
#include "Location.h"

// The MIC cell class
class BMSC
{
 public:
  
  int BMSCLx,BMSCLy,BMSCLz; // The location of each BMSC cell
  int BMSCCurrentTimer;
	 // BMSC cells are always alive, so they do not need as many actions as others
  BMSC (int,int,int,int);
  virtual ~BMSC();

  int ID;  // The ID number of the BMSC
  int CellType; // The type of the cell.
};

#endif /*BMSC_H_*/
