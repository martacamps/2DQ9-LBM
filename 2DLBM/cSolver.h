#pragma once

#define GAME_WIDTH	1000
#define GAME_HEIGHT 1000

class cSolver
{
public:
	cSolver();
	virtual ~cSolver();

	//Create();
	//Initialize();
	virtual bool Input();
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
	int t;             //Number of loops performed by the program
	 
	
};

