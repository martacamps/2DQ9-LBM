//Copyright(c) 2015 Marta Camps Santasmasas
//http://martacamps.weebly.com/
//
//This software is provided 'as-is', without any express or implied
//warranty.In no event will the authors be held liable for any damages
//arising from the use of this software.
//
//Permission is granted to anyone to use this software for any purpose,
//including commercial applications, and to alter it and redistribute it
//freely, subject to the following restrictions :
//
//1. The origin of this software must not be misrepresented; you must not
//claim that you wrote the original software.If you use this software
//in a product, an acknowledgement in the product documentation would be
//appreciated but is not required.
//2. Altered source versions must be plainly marked as such, and must not be
//misrepresented as being the original software.
//3. This notice may not be removed or altered from any source distribution.


#include "stdafx.h"
#include "cColourGraph.h"


cColourGraph::cColourGraph()
{

}

void cColourGraph::SetColours(std::vector<double> *colours, int numColours)
{
	//Check that we have the same number of colours and intervals. 
	if (colours->size() / 3. != m_intervals.size())
	{
		std::ostringstream strs;
		strs << "The number of colours ( " << colours->size() / 3. << " ) does not match the number of intervals (" << m_intervals.size() << ")";
		std::string str = strs.str();
		throw(str);
	}

	m_colours.resize((numColours - 1) * 3);
	m_diffColours.resize((numColours - 1) * 3);

	for (int i = 0; i < numColours - 1; i++)
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

bool cColourGraph::PickColour(double value, double *colourx, double* coloury, double* colourz)
{
	//The value is out of range
	if ((value < m_intervals[0]) || (value > m_intervals.back()))
		return false;

	//Find the interval of the value. 
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

	//Choose the Colour by linear interpolation between the two colours in the interval. 
	double lambda = value / (m_intervals[ihigh] - m_intervals[ilow]) - m_intervals[ilow] / (m_intervals[ihigh] - m_intervals[ilow]);
	*colourx = m_colours[3*ilow]+lambda*m_diffColours[3*ilow];
	*coloury = m_colours[3 * ilow + 1] + lambda*m_diffColours[3 * ilow + 1];
	*colourz = m_colours[3 * ilow + 2] + lambda*m_diffColours[3 * ilow + 2];

	return true;

}

void cColourGraph::DrawScale(std::string name)
{
	float numIntervals = m_intervals.size() - 1;
	//Draw colour bar
	glBegin(GL_TRIANGLE_STRIP);
	for (int i = 0; i < m_intervals.size() - 1; i++)
	{
		glColor3f(m_colours[3 * i], m_colours[3 * i + 1], m_colours[3 * i + 2]);
		glVertex2f(0.0f, (i*0.9) / (numIntervals));
		glColor3f(m_colours[3 * i], m_colours[3 * i + 1], m_colours[3 * i + 2]);
		glVertex2f(0.5f, (i*0.9) / (numIntervals));

		glColor3f(m_colours[3 * i] + m_diffColours[3 * i], m_colours[3 * i + 1] + m_diffColours[3 * i + 1], m_colours[3 * i + 2] + m_diffColours[3 * i + 2]);
		glVertex2f(0.0f, ((i + 1)*0.9) / (numIntervals));
		glColor3f(m_colours[3 * i] + m_diffColours[3 * i], m_colours[3 * i + 1] + m_diffColours[3 * i + 1], m_colours[3 * i + 2] + m_diffColours[3 * i + 2]);
		glVertex2f(0.5f, ((i + 1)*0.9) / (numIntervals));
	}
	glEnd();

	//Draw title and values
	glColor3f(1., 1., 1.);
	DrawText<std::string>(GLUT_BITMAP_HELVETICA_18, 0.0f, 0.93, &name);
	
	for (int i = 0; i < m_intervals.size(); i++)
	{
		DrawText<double>(GLUT_BITMAP_HELVETICA_18, 0.55f, (i*0.9) / (numIntervals), &m_intervals[i]);
	}	
}


template<class T> void cColourGraph::DrawText(GLvoid *fontStyle, float posx, float posy, T *text)
{
	//Convert the type into a character array
	std::ostringstream strs;
	strs << *text;
	std::string str = strs.str();
	const char* cstr = str.c_str();

	//Loop over all the positions of cstr and write them
	glRasterPos2f(posx, posy);
	for (int i = 0; i < strlen(cstr); i++)
	{
		glutBitmapCharacter(fontStyle, cstr[i]);
	}
}