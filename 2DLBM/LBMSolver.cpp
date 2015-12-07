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

void LBMSolver::Create(double nu, double sig, double fx, double fy, double rho, double dx, double size, double time)
{
	cellSize = dx;
	//lidSpeed = v;
	
	//Compute necessary information from input parameters
	double gMag = sqrt(fx*fx + fy*fy);
	if (gMag > 0)
		dt = sqrt((0.5e-3*cellSize) / sqrt(fx*fx + fy*fy));    // So the gravity acceleration in the model is not larger than 1e-4.
	else
		dt = cellSize;
	g[0] = fx*(dt*dt) / cellSize;
	g[1] = fy*(dt*dt) / cellSize;
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
	//p = 101300 * (dm / (cellSize*dt*dt));

	numCells = (int)std::round(size / cellSize);
	numSteps = (int)std::round(time / dt);

	//Create LBM mesh. 
	mesh = new LBMCell*[2];
	mesh[current] = new LBMCell[numCells*numCells];
	mesh[other] = new LBMCell[numCells*numCells];

	tags = new cellTag[numCells*numCells];
	for (int i = 0; i < numCells*numCells; i++)
		tags[i] = gas;

}

void LBMSolver::InitialField()  
{
	int fsizex = numCells / 3;
	int fsizey = 3 * numCells / 4;

	//Set the cell tags for the Boundary Conditions.
	for (int j = 0; j < numCells; j++)
	{
		tags[index(0, j)] = noslipbc;
		tags[index(numCells - 1, j)] = noslipbc;
		
		tags[index(0, j)] = noslipbc;
		tags[index(numCells - 1, j)] = noslipbc;
	}
	for (int i = 1; i < (numCells - 1); i++)
	{
		tags[index(i, 0)] = noslipbc;
		tags[index(i, numCells - 1)] = noslipbc;

		tags[index(i, 0)] = noslipbc;
		tags[index(i, numCells - 1)] = noslipbc;
	}

	//Set the initial position of the fluid, gas and interface cells. 
	for (int i = 1; i < fsizey; i++)
	{
		for (int j = 1; j < fsizex; j++)
		{
			tags[index(i, j)] = fluid;
		}
		tags[index(i, fsizex)] = interface;
		mesh[current][index(i, fsizex)].mass = 0.5;
		interfaceCells.push_back(index(i, fsizex));
		tags[index(i, fsizex)] = interface;
		mesh[other][index(i, fsizex)].mass = 0.5;
	}
	for (int j = 1; j <= fsizex; j++)
	{
		tags[index(fsizey, j)] = interface;
		mesh[current][index(fsizey, j)].mass =  0.5;
		interfaceCells.push_back(index(fsizey, j));
		tags[index(fsizey, j)] = interface;
		mesh[other][index(fsizey, j)].mass =  0.5;
	}
}

void LBMSolver::TimeStep(double t)
{
	//TimeStep
	double epsilon, nextEpsilon, f0, f0inv,tempMass;

	int ue = 0;
	for (int i = 1; i < numCells-1; i++)
	{
		for (int j = 1; j < numCells-1; j++)
		{
			int ij = index(i, j);

			if (tags[ij] != gas)
			{

				switch(tags[ij])
				{
				case fluid:
					//Streaming. 
					for (int l = 0; l < finv.size(); l++)
					{
						int inv = finv[l];
						int previous = index(i + ey[inv], j + ex[inv]);

						//check the type of the neighbour cell and perform streaming
						if (tags[previous] == noslipbc)  //The neighbour cell is a non slip wall
						{
							mesh[current][ij].f[l] = mesh[current][ij].f[inv];
						}
						else										//The neighbour cell is fluid
						{
							mesh[current][ij].f[l] = mesh[other][previous].f[l];
						}
					}
					break;
				case interface :
					tempMass = 0.0;
					//Streaming.
					mesh[current][ij].f[0] = mesh[other][ij].f[0];  //Stream
					for (int l = 1; l < finv.size(); l++)
					{
						int inv = finv[l];
						int previous = index(i + ey[inv], j + ex[inv]);
						int next = index(i + ey[l], j + ex[l]);

						//check the type of the neighbour cell and perform streaming. //I COULD DO THIS WITH A FOREACH IF I PUT THE LINES OF CODE INTO A FUNCTION (MAYBE IT IS MORE EFFICIENT)
						switch (tags[previous])
						{
						case fluid:
							mesh[current][ij].f[l] = mesh[other][previous].f[l];  //Stream
							break;
						case interface:   //THE SAME AS FLUID (I THINK) I CAN PUT THEM TOGETHER.
							mesh[current][ij].f[l] = mesh[other][previous].f[l];  //Stream
							break;
						case gas:
							f0inv = Fequi(mesh[other][ij].u[0], mesh[other][ij].u[1], 1.0, inv);
							f0 = Fequi(mesh[other][ij].u[0], mesh[other][ij].u[1], 1.0, l);
							mesh[current][ij].f[l] = f0inv + f0 - mesh[other][ij].f[l];    //Surface reconstruction. FALTA LA SURFACE TENSION I LA PRESSURE. 
							break;
						case noslipbc:    //The neighbour cell is a non slip wall
							mesh[current][ij].f[l] = mesh[current][ij].f[inv];
							break;
						}

						//Mass update
						if (tags[next] == interface ) //|| tags[next] == fluid)
						{
							epsilon = mesh[other][ij].mass / mesh[other][ij].rho;  //Fill level of the current cell
							nextEpsilon = mesh[other][next].mass / mesh[other][next].rho;  // Fill level of the next cell
							//tempMass += nextEpsilon*mesh[other][next].f[inv] - epsilon*mesh[other][ij].f[l];
							tempMass += 0.5*(nextEpsilon+epsilon)*(mesh[other][next].f[inv] - mesh[other][ij].f[l]);
						}
						if (tags[next] == fluid)
						{
							tempMass += mesh[other][next].f[inv] - mesh[other][ij].f[l];
						}
					}
					mesh[current][ij].mass = mesh[other][ij].mass + tempMass;
					break;
				}

				//Collision: calculate cell density and velocity
				double rho = 0.0, ux = 0.0, uy = 0.0;
				double fi;
				for (int l = 0; l < finv.size(); l++)
				{
					fi = mesh[current][ij].f[l];
					rho += fi;
					ux += fi*ex[l];
					uy += fi*ey[l];
				}
				ux = (ux*c) / rho;
				uy = (uy*c) / rho;
				//}
				mesh[current][ij].rho = rho;
				mesh[current][ij].u[0] = ux;
				mesh[current][ij].u[1] = uy;

				//Collision: apply gravity
				ux += tau*g[0];
				uy += tau*g[1];

				//Collision
				for (int l = 0; l < finv.size(); l++)
				{
					//double eDotu = ex[l] * ux + ey[l] * uy;
					//double f0 = w[l] * rho*(1 - (3.*(ux*ux + uy*uy)) / (2.*c*c) + (3.*eDotu) / c + (9.*eDotu*eDotu) / (2.*c*c));
					double f0 = Fequi(ux, uy, rho, l);
					mesh[current][ij].f[l] = mesh[current][ij].f[l] - (1 / tau)*(mesh[current][ij].f[l] - f0);
				}

				//Collision: Calculate velocities for display and particle tracing. 
				mesh[current][ij].u[0] += 0.5*tau*g[0];
				mesh[current][ij].u[1] += 0.5*tau*g[1];
			}
		}
	}

	//Prepare time step: Mark the interface cells that are empty or full
	std::list<int>::iterator it, sideIt;
	for (it = interfaceCells.begin(); it != interfaceCells.end(); ++it)
	{
		if (mesh[current][*it].mass <= 0)
		{
			tags[*it] = iempty;
		}
			
		if (mesh[current][*it].mass >= mesh[current][*it].rho)
		{
			tags[*it] = ifull;
		}
	}

	//Prepare time step: Mass exchange 
	//Done for the 'current' mesh. 
	double excessMass = 0.0;
	for (it = interfaceCells.begin(); it != interfaceCells.end(); ++it)
	{
		//Calculate the position of the cell to be able to search its neighbours. 
		int i = *it / numCells;
		int j = *it - numCells*i;
		int tag = tags[*it];
		std::list<int> sideInterfaces;
		double lackMass = 0.;
		double extraMass = 0.0;

		if (tags[*it] == ifull)
		{
			//Search its neighbours and count and change the surrounding cells to interface cells.
			for (int l = 1; l < 9; l++)
			{
				int inv = finv[l];
				int previous = index(i + ey[inv], j + ex[inv]);
				int previous2 = index(i + 2 * ey[inv], j + 2 * ex[inv]);

				switch (tags[previous])
				{
				case interface :
					if ((l >= 1 && l <= 4) && ((tags[previous2] == fluid) || (tags[previous2] == noslipbc)))  // If the cell next to the interface is fluid or a wall, I have to fill the interface cell with fluid.
					{
						tags[previous] = ifull;
						lackMass += mesh[current][previous].rho - mesh[current][previous].mass;
					}
					else
					{
						//sideInterfaces.push_back(previous);
					}
					break;
				case gas:  
					tags[previous] = interface;
					mesh[current][previous].f = mesh[current][*it].f;  //Copy the distribution functions from the filled interface cell. 
					mesh[current][previous].mass = 0.0;
					interfaceCells.push_back(previous);
					sideInterfaces.push_back(previous);
				}
			}
		}
		if (tags[*it] == iempty)
		{ 
			//Search its neighbours and count and change the surrounding cells to interface cells.
			for (int l = 1; l < 9; l++)
			{
				int inv = finv[l];
				int previous = index(i + ey[inv], j + ex[inv]);
				int previous2 = index(i + 2 * ey[inv], j + 2 * ex[inv]);

				switch (tags[previous])
				{
				case interface :
					if ((l>=1 && l<=4) && ((tags[previous2] == gas) || (tags[previous2]==noslipbc)))  // If the cell next to the interface is fluid or a wall, I have to fill the interface cell with fluid.
					{
						tags[previous] = iempty;
						extraMass += mesh[current][previous].mass;
					}
					else
					{
						//sideInterfaces.push_back(previous);
					}
					break;
				case fluid:  
					tags[previous] = interface;
					mesh[current][previous].mass = mesh[current][previous].rho;
					interfaceCells.push_back(previous);
					sideInterfaces.push_back(previous);
					break;
				}
			}
		}

		//Calculate the extra mass and share it among the surrounding interface cells.
		double mass = mesh[current][*it].mass;
		if (mass<0 || mass > mesh[current][*it].rho)
		{
			//Calculate the extra mass for this cell
			if (mass < 0)
				extraMass += mass-lackMass;
			else
				extraMass += (mass - mesh[current][*it].rho)-lackMass;

			//Share it between the surrounding interface cells or keep it to put it to all interface cells if there are no surrounding cells. 
			int sideSize = sideInterfaces.size();
			if (sideSize == 0)
			{
				excessMass += extraMass;
			}
			else
			{
				for (sideIt = sideInterfaces.begin(); sideIt != sideInterfaces.end(); ++sideIt)
					mesh[current][*sideIt].mass += extraMass / (double)sideSize;
			}
		}

	}

	//Prepare timestep: Cell type update and remove cells from interface list.
	std::list<int> notInterfaceCells;
	for (it = interfaceCells.begin(); it != interfaceCells.end(); ++it)
	{
		switch (tags[*it])
		{
		case ifull:
			tags[*it] = fluid;
			mesh[current][*it].mass = 1.0;
			notInterfaceCells.push_back(*it);
			break;
		case iempty:
			tags[*it] = gas;
			notInterfaceCells.push_back(*it);
			mesh[current][*it].u[0] = 0.0;
			mesh[current][*it].u[1] = 0.0;
			mesh[current][*it].mass = 0.0;
			mesh[current][*it].rho = 1.0;
			mesh[current][*it].f = { { 4. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 36., 1. / 36., 1. / 36., 1. / 36. } };
			mesh[other][*it].u[0] = 0.0;
			mesh[other][*it].u[1] = 0.0;
			mesh[other][*it].mass = 0.0;
			mesh[other][*it].rho = 1.0;
			mesh[other][*it].f = { { 4. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 36., 1. / 36., 1. / 36., 1. / 36. } };
			break;
		}
	}

	//Prepare next timestep: remove the cells that are no longer interface from the interface list. 
	for (it = notInterfaceCells.begin(); it != notInterfaceCells.end(); ++it)
	{
		interfaceCells.remove(*it);
	}

	//Prepare timestep: If there is excess mass, distribuite it into all the interface cells. 
	if (excessMass > 0)
	{
		int intSize = interfaceCells.size();
		for (it = interfaceCells.begin(); it != interfaceCells.end(); ++it)
		{
			mesh[current][*it].mass += excessMass / (double)intSize;
		}
	}

	//Swap meshes
	other = current;
	current = 1 - other;
	std::cout << t << std::endl;
}

double LBMSolver::Fequi(double ux, double uy, double rho, int l)
{
	double eDotu = ex[l] * ux + ey[l] * uy;
	return w[l] * rho*(1 - (3.*(ux*ux + uy*uy)) / (2.*c*c) + (3.*eDotu) / c + (9.*eDotu*eDotu) / (2.*c*c));
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
		double maxSpeed = 0.1;
		for (int i = 0; i < 6; i++)
			intervals.push_back(i*maxSpeed / 5.0);
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
		intervals = { fluid, gas, interface, ifull, iempty, slipbc, noslipbc };
		colourScale.SetIntervals(&intervals);
		colours = { 0.2, 0.6, 1.0,   1., 1., 1.,    0.,0.,0.8,  1.0,0.,0.,    0.,1.,0.,   0.4,0.5,0.6,     0.4, 0.2, 0. };   //light blue, white, dark blue, red, green, grey, brown.
		colourScale.SetColours(&colours, 7);
		title = "cell type";
	}

	//Draw and colour the cells of the domain.
	double psize = 1;
	glPointSize(psize);
	glPushMatrix();
	glScalef((0.9*win_width) / (psize*numCells), win_height / (psize*numCells), 1.);
	glShadeModel(GL_FLAT);
	std::vector<double> ux, uy, mod;
	double colourx, coloury, colourz;
	bool inside;
	for (int i = 0; i < numCells; i++)
	{
		glBegin(GL_TRIANGLE_STRIP);
		for (int j = 0; j < numCells; j++)
		{
			ux.push_back(mesh[current][index(i, j)].u[0] / c);
			uy.push_back(mesh[current][index(i, j)].u[1] / c);
			mod.push_back(sqrt(ux.back()*ux.back() + uy.back()*uy.back()));
			if (vis[2])	 //Colours by velocity
			{
				inside = colourScale.PickColour(mod.back(), &colourx, &coloury, &colourz);
			}
			else	//Colours by cell tag
			{
				inside = colourScale.PickColour((double)tags[index(i, j)], &colourx, &coloury, &colourz);
			}
			if (!inside)
			{
				colourx = 0.0; coloury = 0.0; colourz = 0.0;
			}
			glColor3f(colourx, coloury, colourz);
			glVertex2f(j, i + psize);
			glVertex2f(j, i);
			glVertex2f(j + 1, i + psize);
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
