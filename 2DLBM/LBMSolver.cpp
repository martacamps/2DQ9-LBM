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
	
	//Compute necessary information from input parameters
	dt = sqrt((1e-4*cellSize) / sqrt(fx*fx+fy*fy));    // So g is not larger than 1e-4. 
	//dt = 1e-4 / sqrt(fx*fx + fy*fy);
	//cellSize = dt;
	//fx = 0;
	//fy = 0;
	//dt = dx;
	g[0] = fx*(dt*dt) / cellSize;
	g[1] = fy*(dt*dt) / cellSize;
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

	double Re = (0.4*size) / nu;

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
	//for (int t = 0; t < numSteps; t++)
	//{
		for (int i = 1; i < numCells-1; i++)
		{
			for (int j = 1; j < numCells-1; j++)
			{
				int ij = index(i, j);
				//MASS EXCHANGE AND CELL TYPE UPDATE. 

				//THERE IS AN IF (OR A CASE) MISSING,SINCE GAS, LIQUID AND INTERFACE CELLS HAVE TO BE TREATED DIFFERENTLY. 

				//Streaming. (I HAVE TO CHECK IF I CAN DO ALL THE L LOOPS IN ONE OR NOT. I KNOW I CAN'T IN THE COLLISION STEP, BUT MAYBE I CAN IN ALL THE OTHER STEPS. 
				for (int l = 0; l < finv.size(); l++)
				{
					int inv = finv[l];
					int previous = index(i + ey[inv], j + ex[inv]);
					//check the type of the neighbour cell
					if (mesh[other][previous].tag == NOSLIPBC)
					{
						mesh[current][ij].f[l] = mesh[current][ij].f[inv];
					}
					else if (mesh[other][previous].tag == SLIPBC)
					{
						//IMPLEMENT SLIP BC
						//IT HAS TO HAVE CHECKS FOR INTERFACE CELLS ALSO, AND MAYBE FOR SOME MORE.
					}
					else   //The neighbour cell is fluid
					{
						mesh[current][ij].f[l] = mesh[other][previous].f[l];
					}     
					//DISTRIBUTION FUNCTIONS RECONSTRUCTION AND SUPERFICIAL FORCES.
					//MASS UPDATE AND MARK EMPTY/FULL CELL
					//BUBBLE VOLUME CALCULATION
				}

				//Collision: density and velocity 
				double rho = 0.0, ux = 0.0, uy = 0.0;
				if (mesh[current][ij].BC == FIXEDV)   //THIS PART IS JUST TO TEST. IT WILL NOT BE HERE FOR THE DEFINITIVE MODEL, SINCE IT DOES NOT CONSERVE MASS.
				{
					rho = 1.; ux = 0.4; uy = 0.0;
				}
				else
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

				//Collision: apply body forces
				if (mesh[current][ij].BC != FIXEDV)  
				{
					ux += tau*g[0];
					uy += tau*g[1];
				}
				
				//Collision: equilibrium function (Eq 3.57 of A single-phase free surface LBM by Nils Thurey)
				double c = cellSize / dt;
				for (int l = 0; l < finv.size(); l++)
				{
					double eDotu = ex[l] * ux + ey[l] * uy;
					//double f0 = w[l] *( rho - 3. / 2.*(ux*ux + uy*uy)/(c*c) + 3.*eDotu/(c*c) + 9. / 2.*(eDotu*eDotu)/(c*c*c*c));

					//double f0 = w[l] *rho*(1 - (3.*(ux*ux + uy*uy)) / (2*c*c) + (3.*eDotu) / (c*c) + (9.*eDotu*eDotu) / (2.*c*c*c*c));  //Eq 3.45 of A single-phase free surface LBM by Nils Thurey
					//double f0 = w[l] * rho*(1 - (3.*(ux*ux + uy*uy)) / (2 * c*c) + (3.*eDotu) / (c) + (9.*eDotu*eDotu) / (2.*c*c));		
					
					
					//double f0 = w[l] * (rho - 3. / 2.*(ux*ux + uy*uy) + 3.*eDotu  + 9. / 2.*(eDotu*eDotu));
					
					double f0 = w[l] * rho*(1 - (3.*(ux*ux + uy*uy)) / 2. + (3.*eDotu) + (9.*eDotu*eDotu) / 2.);
					mesh[current][ij].f[l] = mesh[current][ij].f[l] - (1 / tau)*(mesh[current][ij].f[l] - f0);
				}

				//Collision: Calculate velocities for display and particle tracing. 
				if (mesh[current][ij].BC != FIXEDV)
				{
					mesh[current][ij].u[0] += 0.5*tau*g[0];
					mesh[current][ij].u[1] += 0.5*tau*g[1];
				}
			}
		}
		//Swap meshes
		other = current;
		current = 1 - other;
		std::cout << t << std::endl;
	//}
}

void LBMSolver::Render()
{
	glClear(GL_COLOR_BUFFER_BIT);
	glLoadIdentity();

	//I have to paint mesh[other]

     // b. points
	 float psize = 6;
     glPointSize(psize);   // 3 pixel point. why it only works outside of glBegin-glEnd.
     glBegin ( GL_POINTS );
	 std::vector<double> ux, uy,mod;
     for(int i=0; i<numCells; i++)  // draw point at every grid intersection
     {
		 for (int j = 0; j < numCells; j++)
		 {
			 ux.push_back(mesh[other][index(i, j)].u[0]);
			 uy.push_back(mesh[other][index(i, j)].u[1]);
		     mod.push_back(sqrt(ux.back()*ux.back() + uy.back()*uy.back()));
			 double lambda;// = mod.back() / 0.2;
			 double colourx, coloury, colourz;
			 if (mod.back() <= 0.04)
			 {
				 lambda = mod.back() / 0.04;
				 colourx = lambda;
				 coloury = lambda;
				 colourz = 1-lambda;
			 }
			 else if ((mod.back()> 0.04) &&  (mod.back()<= 0.08))
			 {
				 lambda = mod.back() / 0.04 - 1;
				 colourx = 1;
				 coloury = 1-0.45*lambda;
				 colourz = 0;
			 }
			 else if ((mod.back()> 0.08) && (mod.back() <= 0.12))
			 {
				 lambda = mod.back() / 0.04 - 2;
				 colourx = 1;
				 coloury = 0.55 - 0.55*lambda;
				 colourz = 0;
			 }
			 else if ((mod.back()> 0.12) && (mod.back() <= 0.16))
			 {
				 lambda = mod.back() / 0.04 - 3;
				 colourx = 1-0.42*lambda;
				 coloury = 0;
				 colourz = 0.83*lambda;
			 }
			 else if ((mod.back()> 0.16) && (mod.back() <= 0.4))
			 {
				 lambda = mod.back() / 0.24 - 0.67;
				 colourx = 0.58 + 0.42*lambda;
				 coloury = lambda;
				 colourz = 0.83+0.17*lambda;
			 }
			 else
			 {
				 colourx = 0.29;
				 coloury = 0;
				 colourz = 0.51;
				 double wtf = mod.back();
			 }
			 
			 if (mesh[other][index(i, j)].tag == NOSLIPBC)
				 glColor3f(0.0, 1.0, 0.0);
			 else
				 glColor3f(colourx, coloury, colourz);

			 glVertex2f(i*psize, j*psize);       
			 /*if (mesh[other][index(i, j)].tag == NOSLIPBC)
			 {
				 glColor3f(1.0, 0.0, 0.0);
			 } 
			 else
			 {
				 if (mesh[other][index(i, j)].BC == FIXEDV)
				 {
					 glColor3f(0.0, 0.0, 1.0);
				 }
				 else
				 {
					glColor3f(0.0, 1.0, 0.0);
				 }
				 
			 }	
			 glVertex2f(j * 20, i * 20);*/
        }
     }
     glEnd();

	 glBegin(GL_LINES);
	 glColor3f(0., 0.0, 0.);
	 for (int i = 0; i < numCells; i++)  // draw point at every grid intersection
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
