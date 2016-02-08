//cPoint operates with vectors of lenght 4 representing the coordinates of a point (or any 3D vector) and a value. 
//The = + - and *(for a scalar) operators use all the 4 positions of a cPoint. 
//the * (for dot product) and ^ operators and mag, normalize and distance only use the first 3 positions of cPoint. 
//IF I NEED THEM I CAN MAKE OPERATORS ++ -- ** ^^ THAT WILL THAKE THE 4 POSITIONS OF CPOINT. 

#pragma once

class cPoint
{
public:
	static const cPoint ZERO;
	static const cPoint UNIT_X;
	static const cPoint UNIT_Y;
	static const cPoint UNIT_Z;

	cPoint();
	cPoint(double px, double py, double pz, double pvalue=0.0);
	//the point coordinates are the number between p1 and p2 have the value indicated in value. If pvalue is between p1 and p2 it sets value =pvalue if, sets value=pvalue+10
	cPoint(const cPoint& p1, const cPoint& p2, double pvalue);
	~cPoint();

	double& operator[](int i) { return m_values[i]; }
	double operator[](int i) const { return m_values[i]; }

	cPoint& operator=(const cPoint& other);
	cPoint operator+(const cPoint& other) const;
	cPoint operator-(const cPoint& other) const;
	cPoint operator*(double mult) const;
	double operator*(const cPoint& other) const;	// dot product
	cPoint operator^(const cPoint& other) const;	// cross product

	cPoint& operator+=(const cPoint& other);
	cPoint& operator-=(const cPoint& other);
	cPoint& operator*=(double mult);
	cPoint& operator^=(const cPoint& other);

	double mag() const;
	void normalize();
	cPoint normalize() const;
	double distance(const cPoint& from, const cPoint& to);
	bool SortX(const cPoint& first, const cPoint& second) { return (first[0] < second[0]); }


private:
	double m_values[4]; //x, y and z coordinates of the point and value of the point. 

};

std::ostream& operator<<(std::ostream& out, const cPoint& vector);

