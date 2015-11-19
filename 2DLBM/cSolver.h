#pragma once

#define GAME_WIDTH	880
#define GAME_HEIGHT 800

#include "cColourGraph.h"


class cSolver
{
public:
	cSolver();
	virtual ~cSolver();

	//Create();
	//Initialize();
	virtual void Input();
	void Init(int i_renderStep);
	bool Loop();
	void Finalize();

	//Input
	void ReadKeyboard(unsigned char key, int x, int y, bool press);
	//void ReadMouse(int button, int state, int x, int y);

	//Process
	virtual void TimeStep(double t);
	//int getNumSteps(){ return numSteps; }

	//Output
	virtual void Render();
	virtual void Save();    //Save the data from the simulation


protected:
	int numSteps, renderStep;
	unsigned char keys[256];
	int t;                            //Number of loops performed by the program
	std::vector<bool> vis;            //Visualization types
	int state;                       //State of the simulation: running, paused or exit. 
	cColourGraph colourScale;        //Current colour scale

	 
	
};

