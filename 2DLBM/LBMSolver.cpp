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
#include "LBMSolver.h"


LBMSolver::LBMSolver() :
w({ { 4. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 36., 1. / 36., 1. / 36., 1. / 36. } }),
ex({ { 0, 1, -1, 0, 0, 1, -1, 1, -1 } }),
ey({ { 0, 0, 0, 1, -1, -1, -1, 1, 1 } }),
finv({ { 0, 2, 1, 4, 3, 8, 7, 6, 5 } }),
current(0),
other(1)
{

}

void LBMSolver::Create(double nu, double sig, double v, double rho, double dx, double size, double time)
{
	cellSize = dx;
	lidSpeed = v;
	
	//Compute necessary information from input parameters
	dt = sqrt((1e-4*cellSize) / 9.81);    // So the gravity acceleration in the model is not larger than 1e-4.
	g = -9.81*(dt*dt) / cellSize;
	c = cellSize/dt;

	double nuStar = nu*(dt / (cellSize*cellSize));
	tau = (6 * nuStar + 1) / 2;
	if (tau < 0.5 || tau > 2.5)
	{
		std::ostringstream strs;
		strs << "Relaxation time = " << tau << " out of range, simulation will be too unstable";
		std::string str = strs.str();
		throw(str);
	}
	double dm = rho*cellSize*cellSize*cellSize;
	sigma = sig*(dt*dt) / dm;
	numCells = (int)std::round(size / cellSize);
	numSteps = (int)std::round(time / dt);

	//Create LBM mesh. 
	mesh = new LBMCell*[2];
	mesh[current] = new LBMCell[numCells*numCells];
	mesh[other] = new LBMCell[numCells*numCells];

}

void LBMSolver::InitialField()  
{
	//Set the cell tags for the Boundary Conditions.
	for (int j = 0; j < numCells; j++)
	{
		mesh[current][index(0, j)].tag = NOSLIPBC;
		mesh[current][index(numCells - 1, j)].tag = NOSLIPBC;
		
		mesh[other][index(0, j)].tag = NOSLIPBC;
		mesh[other][index(numCells - 1, j)].tag = NOSLIPBC;
	}
	for (int i = 1; i < (numCells - 1); i++)
	{

		mesh[current][index(i, 0)].tag = NOSLIPBC;
		mesh[current][index(i, numCells - 1)].tag = NOSLIPBC;
		mesh[other][index(i, 0)].tag = NOSLIPBC;
		mesh[other][index(i, numCells - 1)].tag = NOSLIPBC;
	}
	for (int j = 1; j < (numCells - 1); j++)
	{
		mesh[current][index(numCells - 2, j)].BC = FIXEDV;
		mesh[other][index(numCells - 2, j)].BC = FIXEDV;
	}
}

void LBMSolver::TimeStep(double t)
{
	
	for (int i = 1; i < numCells-1; i++)
	{
		for (int j = 1; j < numCells-1; j++)
		{
			int ij = index(i, j);

			//Streaming. 
			for (int l = 0; l < finv.size(); l++)
			{
				int inv = finv[l];
				int previous = index(i + ey[inv], j + ex[inv]);

				//check the type of the neighbour cell and perform streaming
				if (mesh[other][previous].tag == NOSLIPBC)  //The neighbour cell is a non slip wall
				{
					mesh[current][ij].f[l] = mesh[current][ij].f[inv];
				}
				else										//The neighbour cell is fluid
				{
					mesh[current][ij].f[l] = mesh[other][previous].f[l];
				}     
			}

			//Collision
			double rho = 0.0, ux = 0.0, uy = 0.0;
			if (mesh[current][ij].BC == FIXEDV)         //Fix the density and velocity for the cells marked as FIXEDV (the lid cells)
			{
				rho = 1.; ux = lidSpeed; uy = 0.0;
			}
			else                                       //calculate density and velocity 
			{
				double fi;
				for (int l = 0; l < finv.size(); l++)
				{
					fi = mesh[current][ij].f[l];
					rho += fi;
					ux += fi*ex[l];
					uy += fi*ey[l];
				}
				ux = (ux*c)/rho;
				uy = (uy*c)/rho;
			}
			mesh[current][ij].rho = rho;
			mesh[current][ij].u[0] = ux;
			mesh[current][ij].u[1] = uy;

			//Collision: apply gravity
			if (mesh[current][ij].BC != FIXEDV)
			{
				uy += tau*g;
			}
		
			//Collision
			for (int l = 0; l < finv.size(); l++)
			{
				double eDotu = ex[l] * ux + ey[l] * uy;
				double f0 = w[l] *rho*(1 - (3.*(ux*ux + uy*uy)) / (2.*c*c) + (3.*eDotu)/c  + (9.*eDotu*eDotu) / (2.*c*c));  
				mesh[current][ij].f[l] = mesh[current][ij].f[l] - (1 / tau)*(mesh[current][ij].f[l] - f0);
			}

			//Collision: Calculate velocities for display and particle tracing. 
			if (mesh[current][ij].BC != FIXEDV)
			{
				mesh[current][ij].u[1] += 0.5*tau*g;
			}
		}
	}
	//Swap meshes
	other = current;
	current = 1 - other;
}

void LBMSolver::Render()
{
	float win_width = glutGet(GLUT_WINDOW_WIDTH);
	float win_height = glutGet(GLUT_WINDOW_HEIGHT);

	glClear(GL_COLOR_BUFFER_BIT);
	glLoadIdentity();

	//Set colours and intervals for the selected plot. 
	std::vector<double> intervals, colours;
	std::string title;
	if (vis[2])  //Colours by velocity
	{
		for (int i = 0; i < 6; i++)
			intervals.push_back(i*lidSpeed / 5.0);
		colourScale.SetIntervals(&intervals);

		if (vis[1])   //Colour pallette 1
			colours = { 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, };   //blue, cyan, green, yellow, red, magenta
		else          //Colour pallette 2
			colours = { 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1 };   //magenta, red, yellow, green, cyan, blue
		colourScale.SetColours(&colours, 6);

		title = "speed[m/s]";
	}
	else       //Colours by cell tag
	{
		intervals = { FLUID, NOSLIPBC };
		colourScale.SetIntervals(&intervals);
		colours = { 0.2, 0.6, 1.0, 0.4, 0.2, 0. };   //light blue, brown.
		colourScale.SetColours(&colours, 2);
		title = "cell type";
	}

	//Draw and colour the cells of the domain.
	glPushMatrix();
	glScalef((0.9*win_width) / numCells, win_height / numCells, 1.);
	glShadeModel(GL_FLAT);
	std::vector<double> ux, uy, mod;
	double colourx, coloury, colourz;
	bool inside;
	for (int i = 0; i < numCells; i++)
	{
		glBegin(GL_TRIANGLE_STRIP);
		for (int j = 0; j < numCells; j++)
		{
			//Extract velocity. The velocity has not to be divided by c for the lid. 
			if (mesh[other][index(i, j)].BC == FIXEDV)
			{
				ux.push_back(mesh[other][index(i, j)].u[0]);
				uy.push_back(mesh[other][index(i, j)].u[1]);
			}
			else
			{
				ux.push_back(mesh[other][index(i, j)].u[0] / c);
				uy.push_back(mesh[other][index(i, j)].u[1] / c);
			}
			mod.push_back(sqrt(ux.back()*ux.back() + uy.back()*uy.back()));
			if (vis[2])	 //Colours by velocity
			{
				inside = colourScale.PickColour(mod.back(), &colourx, &coloury, &colourz);
			}
			else	//Colours by cell tag
			{
				inside = colourScale.PickColour((double)mesh[other][index(i, j)].tag, &colourx, &coloury, &colourz);
			}
			if (!inside)
			{
				colourx = 0.0; coloury = 0.0; colourz = 0.0;
			}
			glColor3f(colourx, coloury, colourz);
			glVertex2f(j, i + 1);
			glVertex2f(j, i);
			glVertex2f(j + 1, i + 1);
			glVertex2f(j + 1, i);

		}
		glEnd();
	}


	//Draw cell velocity vectors
	if (vis[0])
	{
		glColor3f(0., 0.0, 0.);
		for (int i = 0; i < numCells; i++)
		{
			for (int j = 0; j < numCells; j++)
			{
				glBegin(GL_LINE_STRIP);
				glVertex2f(1 + j, 0.5 + i);
				glVertex2f(1 + j + ux[index(i, j)] / mod[index(i, j)], 0.5 + i + uy[index(i, j)] / mod[index(i, j)]);
				glVertex2f(0.8 + j + ux[index(i, j)] / mod[index(i, j)], 0.8 + i + uy[index(i, j)] / mod[index(i, j)]);   //Arrow head
				glEnd();
			}
		}
	}

	glPopMatrix();

	//Draw colour scale 
	glShadeModel(GL_SMOOTH);
	glPushMatrix();
	glTranslatef(0.9*win_width, 0., 0.);
	glScalef(0.1*win_width, win_height, 1.);
	colourScale.DrawScale(title);
	glPopMatrix();

	//Draw time counter
	std::ostringstream strs;
	strs << std::setprecision(2) << t*dt << "s";
	std::string time = strs.str();
	glColor3f(1., 1., 1.);
	colourScale.DrawText<std::string>(GLUT_BITMAP_HELVETICA_18, 0.9*win_width, 0.96*win_height, &time);

	glutSwapBuffers();
}

LBMSolver::~LBMSolver()
{
	if (mesh)
	{
		delete mesh[0];
		delete mesh[1];
		delete mesh;
	}

}
