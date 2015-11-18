#include "stdafx.h"
#include "LBMSolver.h"
#include "cColourGraph.h"

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
	//dt = cellSize;
	dt = sqrt((1e-4*cellSize) / 9.81);    // So the gravity acceleration in the model is not larger than 1e-4.
	g = -9.81*(dt*dt) / cellSize;

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

	double Re = (lidSpeed*size) / nu;

	//Create LBM mesh. 
	mesh = new LBMCell*[2];
	mesh[current] = new LBMCell[numCells*numCells];
	mesh[other] = new LBMCell[numCells*numCells];
}

void LBMSolver::InitialField()  //This can input the type of simulation or the initial field
{
	//Set the tags for the BC.
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
				else   //The neighbour cell is fluid
				{
					mesh[current][ij].f[l] = mesh[other][previous].f[l];
				}     
			}

			//Collision
			double rho = 0.0, ux = 0.0, uy = 0.0;
			if (mesh[current][ij].BC == FIXEDV)   //Fix the density and velocity for the cells marked as FIXEDV (the lid cells)
			{
				rho = 1.; ux = lidSpeed; uy = 0.0;
			}
			else   //calculate density and velocity 
			{
				double fi;
				for (int l = 0; l < finv.size(); l++)
				{
					fi = mesh[current][ij].f[l];
					rho += fi;
					ux += fi*ex[l];
					uy += fi*ey[l];
				}
				ux = ux / rho;
				uy = uy / rho;
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
				double f0 = w[l] *rho*(1 - (3.*(ux*ux + uy*uy)) / 2. + (3.*eDotu)  + (9.*eDotu*eDotu) / 2.);  
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
	std::cout << t << std::endl;
}

void LBMSolver::Render()
{
	glClear(GL_COLOR_BUFFER_BIT);
	glLoadIdentity();

	//I have to paint mesh[other]
	//Create colour scale (inside the render type switch).
	std::vector<double> intervals, colours;
	for (int i = 0; i < 6; i++)
		intervals.push_back(i*lidSpeed / 5.0);
	//intervals.back() += lidSpeed / 10.0;
	colours = { 0, 0, 1, 1, 1, 0, 1, 0.55, 0, 1, 0, 0, 0.58, 0, 0.83, 1., 1., 1. };     //blue, yellow, orange, red, purple, white
	cColourGraph colourScale(&intervals, &colours, 6);


     double psize = glutGet(GLUT_WINDOW_WIDTH)/numCells;
     glPointSize(psize);   
     glBegin ( GL_POINTS );
	 std::vector<double> ux, uy,mod;
	 bool inside;
     for(int i=0; i<numCells; i++)  // draw point at every grid intersection
     {
		 for (int j = 0; j < numCells; j++)
		 {
			 ux.push_back(mesh[other][index(i, j)].u[0]);
			 uy.push_back(mesh[other][index(i, j)].u[1]);
		     mod.push_back(sqrt(ux.back()*ux.back() + uy.back()*uy.back()));
			 double colourx, coloury, colourz;
			 inside = colourScale.pickColour(mod.back(), &colourx, &coloury, &colourz);
			 if (inside)
				 glColor3f(colourx, coloury, colourz);
			 else
				 glColor3f(0., 0., 0.);           //Black colour means value out of range.

			 glVertex2f(i*psize, j*psize);       
        }
     }
     glEnd();

	 glBegin(GL_LINES);
	 glColor3f(0., 0.0, 0.);
	 for (int i = 0; i < numCells; i++) 
	 {
		 for (int j = 0; j < numCells; j++)
		 {
			 glVertex2f(i * psize,j * psize);
			 glVertex2f(i * psize + (uy[index(i, j)] * psize) / mod[index(i, j)], j * psize + (ux[index(i, j)]*psize)/mod[index(i,j)]);
		 }
	 }
	 glEnd();


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
	

	//std::cout << "deleting mesh " << std::endl;

	//DELETE ALSO THE VECTOR SAVED
}
