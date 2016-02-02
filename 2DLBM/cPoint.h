#pragma once
struct cPoint
{
	cPoint();
	cPoint(double px, double py, double pz, double pvalue);
	//the point coordinates are the number between p1 and p2 have the value indicated in value. If pvalue is between p1 and p2 it sets value =pvalue if, sets value=pvalue+10
	cPoint(const cPoint* p1, const cPoint* p2, double pvalue);
	~cPoint();

	bool SortX(const cPoint& first, const cPoint& second) { return (first.x < second.x); }

	double x=0.;        //x coordinate of the point. Initial value = 0
	double y=0.;        //y coordinate of the point. Initial value = 0
	double z=0.;	     //z coordinate of the point. Initial value = 0
	double value=0.;    //value of the point. Initial value = 0

};

