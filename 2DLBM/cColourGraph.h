#pragma once


class cColourGraph
{
public:
	//Just the 5 colours I know for the moment. 
	cColourGraph(std::vector<double>* intervals, std::vector<double>* colours, int numColours);
	cColourGraph();
	~cColourGraph();

	//void setIntervals(double* intervals) { m_intervals = intervals; }
	void drawScale(std::string name);
	void drawScale(std::string name, std::string* intervalNames);

	bool pickColour(double value, double *colourx, double* coloury, double* colourz);

private:
	std::vector<double> m_colours;          //First RGB colour for each interval
	std::vector<double> m_diffColours;      //Second RGB colour for each interval minus first RGB colour for each interval
    std::vector<double> m_intervals;        //Value of each interval.

};

