#include "stdafx.h"
#include "LBMSolver.h"


LBMSolver::LBMSolver(double nu, double sig, double g, double rho, double dx, double size, double time) :
w({ { 4. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 36., 1. / 36., 1. / 36., 1. / 36. } }),
ex({ { 0, 1, -1, 0, 0, 1, -1, 1, -1 } }),
ey({ { 0, 0, 0, 1, -1, -1, -1, 1, 1 } }),
finv({ { 0, 2, 1, 4, 3, 8, 7, 6, 5 } }),
current(0),
other(1),
cellSize(dx),
f(g)
{
	//Compute necessary information from input parameters
	dt = sqrt((1e-4*cellSize) / f);
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

void LBMSolver::init()  //This can input the type of simulation or the initial field
{
	//Set the tags for the BC.
	for (int j = 0; j < numCells; j++)
	{
		mesh[current][index(0, j)].tag = NOSLIPBC;
		mesh[current][index(numCells - 1, j)].tag = NOSLIPBC;
		mesh[current][index(numCells - 2, j)].BC = FIXEDV;
		int ue = index(numCells - 2, j);
		mesh[other][index(0, j)].tag = NOSLIPBC;
		mesh[other][index(numCells - 1, j)].tag = NOSLIPBC;
		mesh[other][index(numCells - 2, j)].BC = FIXEDV;
	}
	for (int i = 1; i < (numCells - 1); i++)
	{
		mesh[current][index(i, 0)].tag = NOSLIPBC;
		mesh[current][index(i, numCells - 1)].tag = NOSLIPBC;
		mesh[other][index(i, 0)].tag = NOSLIPBC;
		mesh[other][index(i, numCells - 1)].tag = NOSLIPBC;
	}
}

void LBMSolver::mainLoop()
{
	for (int t = 0; t < numSteps; t++)
	{
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
					mesh[current][ij].rho = 1;
					mesh[current][ij].u[0] = 0.01;
					mesh[current][ij].u[1] = 0;
				}
				else
				{
					for (int l = 0; l < finv.size(); l++)
					{
						double fi = mesh[current][ij].f[l];
						rho += fi;
						ux += fi*ex[l];
						uy += fi*ey[l];
					}
					mesh[current][ij].rho = rho;
					mesh[current][ij].u[0] = ux;
					mesh[current][ij].u[1] = uy;
				}
				

				//COLLISION: APPLY BODY FORCES.

				//Collision: equilibrium function (Eq 3.57 of A single-phase free surface LBM by Nils Thurey)
				double c2 = pow(2, cellSize / dt);
				for (int l = 0; l < finv.size(); l++)
				{
					double eDotu = ex[l] * ux + ey[l] * uy;
					double f0 = w[l] *( rho - 3. / 2.*(ux*ux + uy*uy)/c2 + 3.*eDotu/c2 + 9. / 2.*(eDotu*eDotu)/(c2*c2));
					mesh[current][ij].f[l] = mesh[current][ij].f[l] - 1 / tau*(mesh[current][ij].f[l] - f0);
				}
			}
		}
		//Swap meshes
		other = current;
		current = 1 - other;
		std::cout << t << std::endl;
	}
}

LBMSolver::~LBMSolver()
{
	delete mesh[0];
	delete mesh[1];
	delete mesh;
}
