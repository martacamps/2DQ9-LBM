#include "stdafx.h"
#include "cColourGraph.h"


cColourGraph::cColourGraph()
{
}

cColourGraph::cColourGraph(std::vector<double>* intervals, std::vector<double>* colours, int numColours) : m_intervals(*intervals)
{
	//Check that we have the same number of colours and intervals. 
	if (colours->size() / 3. != intervals->size())
	{
		std::ostringstream strs;
		strs << "The number of colours ( " << colours->size() / 3. << " ) does not match the number of intervals (" << intervals->size() << ")";
		std::string str = strs.str();
		throw(str);
	}


	m_colours.resize((numColours - 1) * 3);
	m_diffColours.resize((numColours - 1) * 3);

	for (int i = 0; i < numColours-1; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			m_colours[3 * i + j] = (*colours)[3 * i + j];
			m_diffColours[3 * i + j] = (*colours)[3 * (i + 1) + j] - (*colours)[3 * i + j];
		}
	}
}

cColourGraph::~cColourGraph()
{
}

bool cColourGraph::pickColour(double value, double *colourx, double* coloury, double* colourz)
{
	//The value is out of range
	if ((value < m_intervals[0]) || (value > m_intervals.back()))
		return false;

	std::vector<double>::iterator low;
	int ilow, ihigh;
	low = std::lower_bound(m_intervals.begin(), m_intervals.end(), value);
	
	ihigh = low - m_intervals.begin();
	if (ihigh == 0)
	{
		ihigh = 1;
		ilow = 0;
	}
	else
		ilow = ihigh - 1;

	//Choose the Colour by linear interpolation between the two colours in the interval of value. 
	double lambda = value / (m_intervals[ihigh] - m_intervals[ilow]) - m_intervals[ilow] / (m_intervals[ihigh] - m_intervals[ilow]);
	*colourx = m_colours[3*ilow]+lambda*m_diffColours[3*ilow];
	*coloury = m_colours[3 * ilow + 1] + lambda*m_diffColours[3 * ilow + 1];
	*colourz = m_colours[3 * ilow + 2] + lambda*m_diffColours[3 * ilow + 2];

	return true;

}

void cColourGraph::drawScale(std::string name)
{
	
}


void cColourGraph::drawScale(std::string name, std::string* intervalNames)
{

}
