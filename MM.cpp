//#include "stdzafx.h"
#include <stdio.h>
#include <stdlib.h>
#include <list>
#include "MM.h"
//#include ".\mm.h"


MM::MM(int IDMM,int x,int y,int z)
{
	ID = IDMM;
    MMLx = x; MMLy = y;	MMLz = z;
	MMApopStatus = APOPN;
	MMCellStatus=StaN;
	MMApopTimer=0;
	MMCycleTimer=0;
	MMWaitProTimer=0;
	MMCurrentTimer=0;
	MMGenNo=0;
	CellType=CTYPEMM;
}
	
   	

MM::~MM()
{
}






 
  
