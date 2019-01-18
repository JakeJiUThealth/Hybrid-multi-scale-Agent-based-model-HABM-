//#include "stdzafx.h"
#include <stdio.h>
#include <stdlib.h>
#include <list>
#include "CTL.h"


CTL::CTL(int IDCTL,int x,int y,int z)
{
	ID = IDCTL;
    CTLLx = x; CTLLy = y;CTLLz = z;
    CTLApopStatus = APOPN;
	CTLGenNo=0;
	CTLCellStatus=StaN;
	CTLApopTimer=0;
	CTLCycleTimer=0;
	CTLWaitProTimer=0;
	CTLCurrentTimer=0;
	CellType=CTYPECTL;
}
	    	

CTL::~CTL()
{
}

