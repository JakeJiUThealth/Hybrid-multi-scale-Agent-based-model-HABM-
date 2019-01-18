//#include "stdzafx.h"
#include <stdio.h>
#include <stdlib.h>
#include <list>
#include "Treg.h"

Treg::Treg(int IDTreg,int x,int y,int z)
{
	ID = IDTreg;
    TregLx = x; TregLy = y;TregLz = z;
    TregApopStatus = APOPN;
	TregGenNo=0;
	TregCellStatus=StaN;
	TregApopTimer=0;
	TregCycleTimer=0;
	TregWaitProTimer=0;
	TregCurrentTimer=0;
	CellType=CTYPETREG;
}
    	

Treg::~Treg()
{
}