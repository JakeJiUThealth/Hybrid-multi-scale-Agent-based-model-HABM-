//#include "stdzafx.h"
#include <stdio.h>
#include <stdlib.h>
#include <list>
#include "BMSC.h"
#define BMSCCELL 5;
//#include ".\bmsc.h"


BMSC::BMSC(int IDBMSC,int x,int y,int z)
{
	ID = IDBMSC;
    BMSCLx = x; BMSCLy = y; BMSCLz = z;
    CellType=BMSCCELL;
	BMSCCurrentTimer=0;
}

		       	

BMSC::~BMSC()
{
}






 
  
