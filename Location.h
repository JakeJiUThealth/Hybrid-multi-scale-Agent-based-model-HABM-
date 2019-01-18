#ifndef LOCATION_H_
#define LOCATION_H_

//The class of location for cells 
class Location
{
public:
 int PosX,PosY,PosZ;
 Location(int,int,int);
 Location (void);
virtual ~Location();
};

#endif /*LOCATION_H_*/
