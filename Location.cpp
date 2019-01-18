//#include "stdafx.h"
#include "Location.h"

Location::Location(int x, int y,int z)
{
    PosX = x;
    PosY = y;
    PosZ = z;

}
Location::Location(void)
{
    PosX = 0;
    PosY = 0;
    PosZ = 0;

}
Location::~Location(void)
{
}
