#include "stdafx.h"
#include "cPoint.h"


cPoint::cPoint()
{
}

cPoint::cPoint(double px, double py, double pz, double pvalue)
	: x(px), y(py), z(pz), value(pvalue)
{
}


//the point coordinates are the number between p1 and p2 have the value indicated in value. If pvalue is between p1 and p2 it sets value =pvalue if, sets value=pvalue+10. 
cPoint::cPoint(const cPoint* p1, const cPoint* p2, double pvalue)
{
	double val1 = p1->value;
	double val2 = p2->value;

	//Check that value is within the range of p1 and p2
	if ((val1 <= val2) && ((pvalue <val1) || (pvalue>val2)))
		value = pvalue+10;
	else if ((val1 > val2) && ((pvalue > val1) || (pvalue < val2)))
		value=pvalue+10;
	else
	{
		//Interpolate 
		value = pvalue;
		double lambda = (value - val1) / (val2 - val1);
		x = p1->x + lambda*(p2->x - p1->x);
		y = p1->y + lambda*(p2->y - p1->y);
		z = p1->z + lambda*(p2->z - p1->z);
	}

}


cPoint::~cPoint()
{
}
