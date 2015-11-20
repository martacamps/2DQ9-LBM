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

//Colour scale and numeric intervals for each colour

#pragma once

class cColourGraph
{
public:
	cColourGraph();
	~cColourGraph();

	//Set the numeric intervals
	void SetIntervals(std::vector<double> *intervals) { m_intervals = *intervals; } 
	//Set the RGB (from 0 to 1) colour code for every interval.
	void SetColours(std::vector<double> *colours, int numColours);                     

	//Draw the colour scale. The string name will be written on top of it.  
	void DrawScale(std::string name);
	//Write any text, numbers...
	template <class T> void DrawText(GLvoid *fontStyle, float posx, float posy, T *text);
	//Colourx, coloury, colourz are the RGB code of the colour corresponding to value. It returns true if value is between the first and the last interval, and false if it isn´t.
	bool PickColour(double value, double *colourx, double* coloury, double* colourz);

private:
	std::vector<double> m_colours;          //First RGB colour for each interval
	std::vector<double> m_diffColours;      //Second RGB colour for each interval minus first RGB colour for each interval
    std::vector<double> m_intervals;        //Value of each interval.

};

