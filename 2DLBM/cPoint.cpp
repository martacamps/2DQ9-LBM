#include "stdafx.h"
#include "cPoint.h"

const cPoint cPoint::ZERO(0.0f, 0.0f, 0.0f, 0.0f);
const cPoint cPoint::UNIT_X(1.0f, 0.0f, 0.0f, 0.0f);
const cPoint cPoint::UNIT_Y(0.0f, 1.0f, 0.0f, 0.0f);
const cPoint cPoint::UNIT_Z(0.0f, 0.0f, 1.0f, 0.0f);


cPoint::cPoint()
{
	for (int i = 0; i < 4; i++)
		m_values[i] = 0;
}

cPoint::cPoint(double px, double py, double pz, double pvalue)
{
	m_values[0] = px;
	m_values[1] = py;
	m_values[2] = pz;
	m_values[3] = pvalue;
}

//the point coordinates are the number between p1 and p2 have the value indicated in value. If pvalue is between p1 and p2 it sets value =pvalue if, sets value=pvalue+10. 
cPoint::cPoint(const cPoint& p1, const cPoint& p2, double pvalue)
{
	double val1 = p1[3];
	double val2 = p2[3];

	//Check that value is within the range of p1 and p2
	if ((val1 <= val2) && ((pvalue <val1) || (pvalue>val2)))
		m_values[3] = pvalue+10;
	else if ((val1 > val2) && ((pvalue > val1) || (pvalue < val2)))
		m_values[3]=pvalue+10;
	else
	{
		//Interpolate 
		m_values[3] = pvalue;
		double lambda = (pvalue - val1) / (val2 - val1);
		for (int i = 0; i < 4; ++i)
			m_values[i] = p1[i] + lambda*(p2[i] - p1[i]);
	}
}

cPoint& cPoint::operator=(const cPoint& other)
{
	for (int i = 0; i < 4; i++)
		m_values[i] = other[i];
	return *this;
}

cPoint cPoint::operator+(const cPoint& other) const
{
	return cPoint(
		m_values[0] + other[0],
		m_values[1] + other[1],
		m_values[2] + other[2],
		m_values[3] + other[3]);
}

cPoint cPoint::operator-(const cPoint& other) const
{
	return cPoint(
		m_values[0] - other[0],
		m_values[1] - other[1],
		m_values[2] - other[2],
		m_values[3] - other[3]);
}

cPoint cPoint::operator*(double mult) const
{
	return cPoint(
		m_values[0] * mult,
		m_values[1] * mult,
		m_values[2] * mult,
		m_values[3] * mult);
}

double cPoint::operator*(const cPoint& other) const
{
	float result = 0.0f;
	for (int i = 0; i < 3; i++)
		result += m_values[i] * other[i];
	return result;
}

cPoint cPoint::operator^(const cPoint& other) const
{
	return cPoint(
		m_values[1] * other[2] - m_values[2] * other[1],
		m_values[2] * other[0] - m_values[0] * other[2],
		m_values[0] * other[1] - m_values[1] * other[0],
		0.0f);
}


cPoint& cPoint::operator+=(const cPoint& other)
{
	for (int i = 0; i < 4; i++)
		m_values[i] += other[i];
	return *this;
}

cPoint& cPoint::operator-=(const cPoint& other)
{
	for (int i = 0; i < 4; i++)
		m_values[i] -= other[i];
	return *this;
}

cPoint& cPoint::operator*=(double mult)
{
	for (int i = 0; i < 4; i++)
		m_values[i] *= mult;
	return *this;
}

cPoint& cPoint::operator^=(const cPoint& other)
{
	cPoint result = *this ^ other;
	for (int i = 0; i < 4; i++)
		m_values[i] = result[i];
	return *this;
}

double cPoint::mag() const
{
	/*float result = 0.0f;
	for (int i = 0; i < 3; i++)
		result += m_values[i] * m_values[i];*/
	return sqrt(*this * *this);
}

void cPoint::normalize()
{
	float m = mag();
	if (m != 0.0f)
		for (int i = 0; i < 3; i++)
			m_values[i] /= m;
}

cPoint cPoint::normalize() const
{
	float m = mag();
	if (m == 0.0f)
		return *this;
	else {
		float values[3];
		for (int i = 0; i < 3; i++)
			values[i] = m_values[i] / m;
		return cPoint(values[0], values[1], values[2], m_values[3]);
	}
}

double cPoint::distance(const cPoint& from, const cPoint& to)
{
	return cPoint(
		from[0] - to[0],
		from[1] - to[1],
		from[2] - to[2],
		from[3] - to[3])
		.mag();
}


cPoint::~cPoint()
{
}

std::ostream& operator<<(std::ostream& out, const cPoint& vector)
{
	return out << "(" << vector[0] << ", " << vector[1] << ", " << vector[2] << ", " << vector[3] << ")";
}
