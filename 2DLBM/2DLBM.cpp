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


// 2DLBM.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include "LBMSolver.h"

LBMSolver cavity;

void AppRender()
{
	cavity.Render();
}
void AppKeyboard(unsigned char key, int x, int y)
{
	cavity.ReadKeyboard(key, x, y, true);
}
void AppIdle()
{
	if (!cavity.Loop()) exit(0);
}

void ChangeSize(int w, int h)
{
	glViewport(0, 0, w, h);      // reset the viewport
	glMatrixMode(GL_PROJECTION); // modify the projection matrix
	glLoadIdentity();            // load an identity matrix into the projection matrix
	glOrtho(0, w, 0, h, 0, 1);   // create new projection matrix

	glMatrixMode(GL_MODELVIEW); // return to the model matrix
	glLoadIdentity();           // load an identity matrix into the model-view matrix
}


int main(int argc, char* argv[])
{
	try
	{
		std::string line;
		int res_x, res_y, pos_x, pos_y, printSteps, testCase;

		//User input and instructions
		if (argc != 3)
		{
			do {
				std::cout << "Drawing frequency: " << std::endl;
				std::cout << "(for example type 10 to draw the results every 10 solver steps)" << std::endl;
				std::cin >> line;
			} while (std::stoi(line) <= 0);
			printSteps = stoi(line);
			std::cout << "Simulation settings: " << std::endl;
			std::cout << "   1: Previous Reynolds number = 250. Normal gravity" << std::endl;
			std::cout << "   2: Previous Reynolds number = 100 low resolution. Gravity in +x" << std::endl;
			std::cout << "   3: Previous Reynolds number = 100 high resolution. Gravity in +x -y" << std::endl;
			std::cout << "   Previous Reynolds number = 400. No gravity" << std::endl;
			std::cout << "   Any other number: Manually input the simulation settings. " << std::endl;
			std::cout << std::endl;
		    std::cout << "  (WARNING: The simulation may crash or show unrealistic results depending on the manually input settings.)" << std::endl;
			std::cout << std::endl;
			std::cin >> line;
			testCase = std::stoi(line);
			std::cout << " NOTE: Press the space bar during the simulation to see the run-time commands" << std::endl;
		}
		else
		{
			printSteps = atoi(argv[1]);
			testCase = atoi(argv[2]);
		}

		//GLUT initialization
		glutInit(&argc, argv);

		//RGBA with double buffer
		glutInitDisplayMode(GLUT_RGBA | GLUT_ALPHA | GLUT_DOUBLE);

		//Create centered window
		res_x = glutGet(GLUT_SCREEN_WIDTH);
		res_y = glutGet(GLUT_SCREEN_HEIGHT);
		pos_x = (res_x >> 1) - (GAME_WIDTH >> 1);
		pos_y = (res_y >> 1) - (GAME_HEIGHT >> 1);

		glutInitWindowPosition(pos_x, pos_y);
		glutInitWindowSize(GAME_WIDTH, GAME_HEIGHT);
		glutCreateWindow("LBM Lid driven cavity flow");

		//Register callback functions
		glutDisplayFunc(AppRender);
		glutKeyboardFunc(AppKeyboard);
		glutReshapeFunc(ChangeSize);
		glutIdleFunc(AppIdle);

		//Solver initializations
		cavity.Init(printSteps);   
		//The solver starts with a default parameter set or a user-defined one. 
		switch (testCase)
		{
		case 1:
			cavity.Create(1e-6, 7e-4, 0, -9.81, 1000, 0.000006, 0.0005, 150);   //Re 250 
			break;
		case 2:
			cavity.Create(1e-6, 7e-4, -9.81,0, 1000, 0.000005, 0.0001, 150);   //Re 100 coarse
			break;
		case 3:
			cavity.Create(1e-6, 7e-4, 4.5, -4.5, 1000, 0.000004, 0.0004, 150);   //Re 100 fine
			break;
		case 4:
			cavity.Create(1e-6, 7e-4, 0, 0, 1000, 0.000006, 0.001, 150);   //Re 400 fine
			break;
		default:
			double nu, fx, fy, dx, size, time;
			std::cout << "Set fluid's kinematic viscosity [m2/s]:" << std::endl;
			std::cin >> line;
			nu = std::stod(line);
			std::cout << "Set the x component of gravity [m/s]:" << std::endl;
			std::cin >> line;
			fx = std::stod(line);
			std::cout << "Set the y component of gravity [m/s]:" << std::endl;
			std::cin >> line;
			fy = std::stod(line);
			std::cout << "Set the cell size [m]:" << std::endl;
			std::cin >> line;
			dx = std::stod(line);
			std::cout << "Set the size of the domain [m]:" << std::endl;
			std::cin >> line;
			size = std::stod(line);
			std::cout << "Set the amount of time to simulate [s]:" << std::endl;
			std::cin >> line;
			time = std::stod(line);
			cavity.Create(nu, 7e-4, fx, fy, 1000, dx, size, time);              //User defined
			break;
		}
	
		cavity.InitialField();

		//Application loop
		glutMainLoop();

	}
	catch (std::exception& e)
	{
		std::cout << e.what() << '\n';
	}
	catch (std::string se)
	{
		std::cout << se << " Interrupting simulation" << '\n';
	}
	

	return 0;
}

