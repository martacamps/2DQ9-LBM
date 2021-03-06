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
#include "cSolver.h"



cSolver::cSolver() : m_t(0), m_vis(1), m_state(1)
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
	m_vis.back() = true;
	for (int i = 0; i < 3; i++)
		m_vis.push_back(true);

	//Solver initialization
	m_renderStep = i_renderStep;
}

void cSolver::Input()
{
	if (m_keys[27])      //ESC to exit
		m_state = 0;
	if (m_keys[32])      //SPACE to pause simulation
	{
 		m_state = -m_state;
		m_keys[32] = 0;
		Render();
	}
	if (m_keys[49])      //Activate / deactivate visualization mode 1
	{
		m_vis[0] = !m_vis[0];
		m_keys[49] = 0;
		Render();
	}
	if (m_keys[50])      //Activate / deactivate visualization mode 2
	{
		m_vis[1] = !m_vis[1];
		m_keys[50] = 0;
		Render();
	}
	if (m_keys[51])      //Activate / deactivate visualization mode 3 
	{
		m_vis[2] = !m_vis[2];
		m_keys[51] = 0;
		Render();
	}


}

bool cSolver::Loop()
{
	bool res = true;

	//Process input
	Input();
	switch (m_state)
	{
	case 0:
		return false;
		break;
	case 1:
		//Simulation loop
		if (m_t <= m_numSteps)
		{
			TimeStep(m_t);
			if (m_t%m_renderStep == 0)
			{
				Render();
				Save();
			}
		}
		else
			Render();     //Continuously render the last image when the simulation finishes. 
		m_t++;
		break;
	case -1:
		PauseMenu();
		break;
	}

	return res;
}

void cSolver::PauseMenu()
{
	//Show pause menu
	float win_width = glutGet(GLUT_WINDOW_WIDTH);
	float win_height = glutGet(GLUT_WINDOW_HEIGHT);

	glLoadIdentity();
	glColor3f(1., 1., 1.);
	std::string text("[SPACE BAR] Continue simulation");
	m_colourScale.DrawText<std::string>(GLUT_BITMAP_HELVETICA_18, 0.1*win_width, 0.8*win_height, &text);
	text = "[1] Hide /show velocity lines";
	m_colourScale.DrawText<std::string>(GLUT_BITMAP_HELVETICA_18, 0.1*win_width, 0.7*win_height, &text);
	text = "[2] Change between colour scale 1 and 2";
	m_colourScale.DrawText<std::string>(GLUT_BITMAP_HELVETICA_18, 0.1*win_width, 0.6*win_height, &text);
	text = "[3] Change between velocity and cell tag view";
	m_colourScale.DrawText<std::string>(GLUT_BITMAP_HELVETICA_18, 0.1*win_width, 0.5*win_height, &text);
	text = "[ESC] Exit program";
	m_colourScale.DrawText<std::string>(GLUT_BITMAP_HELVETICA_18, 0.1*win_width, 0.4*win_height, &text);
	glutSwapBuffers();
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
	m_keys[key] = press;
}

//Output
void cSolver::Render()
{

}

void cSolver::TimeStep(float t)
{

}