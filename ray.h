#ifndef RAY_H
#define RAY_H

#include "Gz.h"

//this class is taken directly from raytracing book 1
class ray
{
public:

	GzCoord orig;
	GzCoord dir;


	ray(){}
	ray(GzCoord& origin, GzCoord direction)
	{
		for (int i = 0; i < 3; i++)
		{
			origin[i] = 0;

			direction[i] = 0; 
		}
	}



	//compute what point on ray value t falls on
	void at(double t, GzCoord result) const {
		
		for (int i = 0; i < 3; i++)
		{
			orig[i] + t * dir[i];
		}
	}

	
};

#endif // !RAY_H

