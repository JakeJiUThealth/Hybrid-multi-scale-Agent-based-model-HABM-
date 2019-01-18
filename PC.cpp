//#include "stdzafx.h"
#include <stdio.h>
#include <stdlib.h>
#include <list>
#include "PC.h"
//#include ".\pc.h"


PC::PC(int IDPC,int x,int y,int z)
{
	ID = IDPC;
    PCLx = x; PCLy = y;PCLz = z;
    PCApopStatus = APOPN;
	PCGenNo=0;
	PCCellStatus=StaN;
	PCApopTimer=0;
	PCCycleTimer=0;
	PCWaitProTimer=0;
	PCCurrentTimer=0;
	CellType=CTYPEPC;
}

	    	

PC::~PC()
{
}






 
  
