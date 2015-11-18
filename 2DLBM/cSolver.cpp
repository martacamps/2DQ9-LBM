#include "stdafx.h"
#include "cSolver.h"



cSolver::cSolver() : t(0)
{
}


cSolver::~cSolver()
{
}

void cSolver::Init(int i_renderStep)
{
	//Graphics initialization
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0, GAME_WIDTH, 0, GAME_HEIGHT, 0, 1);
	glMatrixMode(GL_MODELVIEW);


	//Solver initialization
	renderStep = i_renderStep;
	//Create();
	//Initialize();
}

bool cSolver::Input()
{
	bool res = true;

	if (keys[27])
		res = false;
	//if(keys[ESPAI]) , PAUSA LA SIMULACIO ....
	//if(keys[1], 2, 3, 4... cambiar el mode de visualització . Per exemple 1 pot ser contorns de velocitat i 2 pot ser vectors de velocitat. 

	return res;
}

bool cSolver::Loop()
{
	bool res = true;

	//Process input
	res = Input();
	
	//Simulation loop
	if (t <= numSteps)
	{
		TimeStep(t);
		if (t%renderStep == 0)
		{
			Render();
			Save();
		}
	}
	else
		Render();     //FOR THE MOMENT I WILL JUST CONTINUOUSLY RENDER THE LAST IMAGE WHEN THE SIMULATION FINISHES. 

	t++;

	return res;
}

void cSolver::Save()
{

}

void cSolver::Finalize()
{

}

//Input
void cSolver::ReadKeyboard(unsigned char key, int x, int y, bool press)
{
	keys[key] = press;
}

//Output
void cSolver::Render()
{
	glClear(GL_COLOR_BUFFER_BIT);
	glLoadIdentity();

	glBegin(GL_TRIANGLES);
		glColor3f(1.0f, 0.0f, 0.0f);	glVertex2i(100, 100);
		glColor3f(0.0f, 1.0f, 0.0f);	glVertex2i(300, 100);
		glColor3f(0.0f, 0.0f, 1.0f);	glVertex2i(200, 200);
	glEnd();

	glutSwapBuffers();
}

void cSolver::TimeStep(double t)
{

}