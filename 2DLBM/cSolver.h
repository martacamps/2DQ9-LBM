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

//Base class of a numerical solver.

#pragma once

#define GAME_WIDTH	880
#define GAME_HEIGHT 800

#include "cColourGraph.h"


class cSolver
{
public:
	cSolver();
	virtual ~cSolver();

	//Handles user input
	virtual void Input();
	//Initialize the solver
	void Init(int i_renderStep);
	//Main loop of the solver
	bool Loop();
	void Finalize();

	//Input
	void ReadKeyboard(unsigned char key, int x, int y, bool press);

	//Time step of the numerical solver
	virtual void TimeStep(double t);

	//Output
	virtual void Render();
	//Save the data from the simulation
	virtual void Save();    
	//Draw the pause menu
	virtual void PauseMenu();


protected:
	int numSteps;                     //Number of time steps required to finish the simulation
    int renderStep;                   //Render is called every renderStep time Step
	unsigned char keys[256];          //Keyboard ASCII codes.
	int t;                            //Number of loops performed by the program
	std::vector<bool> vis;            //Visualization modes. 
	int state;                       //State of the simulation: running, paused or exit. 
	cColourGraph colourScale;        //Colour scale used to represent the results of the simulation. 

	 
	
};

