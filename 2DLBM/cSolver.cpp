#include "stdafx.h"
#include "cSolver.h"



cSolver::cSolver() : t(0), vis(1), state(1)
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

	//3 visualization options
	vis.back() = true;
	for (int i = 0; i < 3; i++)
		vis.push_back(true);

	//Solver initialization
	renderStep = i_renderStep;
	//Create();
	//Initialize();
}

void cSolver::Input()
{
	if (keys[27])      //ESC to exit
		state = 0;
	if (keys[32])      //SPACE to pause simulation
	{
 		state = -state;
		keys[32] = 0;
	}
	if (keys[49])      //Activate / deactivate visualization mode 1
	{
		vis[0] = !vis[0];
		keys[49] = 0;
	}
	if (keys[50])      //Activate / deactivate visualization mode 2
	{
		vis[1] = !vis[1];
		keys[50] = 0;
	}
	if (keys[51])      //Activate / deactivate visualization mode 3 
	{
		vis[2] = !vis[2];
		keys[51] = 0;
	}


}

bool cSolver::Loop()
{
	bool res = true;

	//Process input
	Input();
	switch (state)
	{
	case 0:
		return false;
		break;
	case 1:
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
			Render();     //CONTINUOUSLY RENDER THE LAST IMAGE WHEN THE SIMULATION FINISHES. 
		t++;
		break;
	case -1:
		//Show pause menu
		float win_width = glutGet(GLUT_WINDOW_WIDTH);
		float win_height = glutGet(GLUT_WINDOW_HEIGHT);

	    //glClear(GL_COLOR_BUFFER_BIT);
		glLoadIdentity();
		glColor3f(1., 1., 1.);
		std::string text("[SPACE BAR] Continue simulation");
		colourScale.DrawText<std::string>(GLUT_BITMAP_HELVETICA_18, 0.1*win_width, 0.8*win_height, &text);
		text = "[1] Hide /show velocity lines";
		colourScale.DrawText<std::string>(GLUT_BITMAP_HELVETICA_18, 0.1*win_width, 0.7*win_height, &text);
		text = "[2] Change between colour scale 1 and 2";
		colourScale.DrawText<std::string>(GLUT_BITMAP_HELVETICA_18, 0.1*win_width, 0.6*win_height, &text);
		text = "[3] Change between velocity and cell tag view";
		colourScale.DrawText<std::string>(GLUT_BITMAP_HELVETICA_18, 0.1*win_width, 0.5*win_height, &text);
		text = "[ESC] Exit program";
		colourScale.DrawText<std::string>(GLUT_BITMAP_HELVETICA_18, 0.1*win_width, 0.4*win_height, &text);
		glutSwapBuffers();
		break;
	}

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