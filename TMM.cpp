//#include "stdzafx.h"
#include <stdio.h>
#include <stdlib.h>
#include <list>
#include "TMM.h"
//#include ".\tmm.h"


TMM::TMM(int IDTMM,int x,int y,int z)
{
	ID = IDTMM;
    TMMLx = x; TMMLy = y;TMMLz = z;
    TMMApopStatus = APOPN;
	TMMCellStatus=StaN;
	TMMApopTimer=0;
	TMMCycleTimer=0;
	TMMWaitProTimer=0;
	TMMCurrentTimer=0;
	CellType=CTYPETMM;
}

		    	
TMM::~TMM()
{
}






 
  
