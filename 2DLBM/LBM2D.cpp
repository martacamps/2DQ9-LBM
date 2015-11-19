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
	glViewport(0, 0, w, h); // reset the viewport
	glMatrixMode(GL_PROJECTION); // modify the projection matrix
	glLoadIdentity();            // load an identity matrix into the projection matrix
	glOrtho(0, w, 0, h, 0, 1); // create new projection matrix

	/// Important!!! You need to switch back to the model-view matrix
	/// or else your OpenGL calls are modifying the projection matrix!
	glMatrixMode(GL_MODELVIEW); // return to the model matrix
	glLoadIdentity();           // load an identity matrix into the model-view matrix
}


int main(int argc, char* argv[])
{
	try
	{
		int res_x, res_y, pos_x, pos_y;

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
		//glutKeyboardUpFunc(AppKeyboardUp);
		//glutSpecialFunc(AppSpecialKeys);
		//glutSpecialUpFunc(AppSpecialKeysUp);
		//glutMouseFunc(AppMouse);
		glutIdleFunc(AppIdle);

		//Game initializations
		cavity.Init(10);    //CREATE AND INITIALIZE WILL BE PART OF INIT WHEN THEY READ THEIR INPUTS FROM FILES. THIS WAY ALL OF THEM WILL NEED THE SAME INPUT PARAMETER (THE NAME OF THE FILE). 
		cavity.Create(1e-6, 7e-4, 0.4, 1000, 0.000006, 0.001, 150);
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

