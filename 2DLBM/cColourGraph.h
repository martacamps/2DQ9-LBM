#pragma once


class cColourGraph
{
public:
	cColourGraph();
	~cColourGraph();

	void SetIntervals(std::vector<double> *intervals) { m_intervals = *intervals; }
	void SetColours(std::vector<double> *colours, int numColours);

	void DrawScale(std::string name);
	void DrawScale(std::string name, std::vector<std::string>* intervalNames, float posx, float posy, float sizex, float sizey, float psize);
	template <class T> void DrawText(GLvoid *fontStyle, float posx, float posy, T *text);

	bool PickColour(double value, double *colourx, double* coloury, double* colourz);

private:
	std::vector<double> m_colours;          //First RGB colour for each interval
	std::vector<double> m_diffColours;      //Second RGB colour for each interval minus first RGB colour for each interval
    std::vector<double> m_intervals;        //Value of each interval.

};

