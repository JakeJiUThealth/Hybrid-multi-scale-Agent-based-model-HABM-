//#include "stdzafx.h"
#include <stdio.h>
#include <stdlib.h>
#include <list>
#include "MIC.h"
//#include ".\mic.h"


MIC::MIC(int IDMIC,int x,int y,int z)
{
	ID = IDMIC;
    MICLx = x; MICLy = y; MICLz = z;
    MICApopStatus = APOPN;
	MICCellStatus=StaN;
	MICApopTimer=0;
	MICCycleTimer=0;
	MICWaitProTimer=0;
	MICCurrentTimer=0;
	CellType=CTYPEMIC;
}

		       	

MIC::~MIC()
{
}






 
  
