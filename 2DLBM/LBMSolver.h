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

//2DQ9 LBM Solver for a lid driven cavity flow


#pragma once

#include "LBMCell.h"
#include "cSolver.h"
#include "cPoint.h"

#define BUFFER 1e-2

enum cellTag
{
	fluid, gas, interface, ifull, iempty, inew, ifluid, igas, slipbc, noslipbc
};



class LBMSolver : public cSolver
{
public:
	LBMSolver();

	//Input the solver parameters and prepare the solver to start
	void Create(double nu, double sig, double fx, double fy, double rho, double dx, double size, double time);
	//Set the initial conditions for a lid driven cavity flow
	void InitialField();       
	//2DQ9 LBM time step
	void TimeStep(double t);
	//
	//Draws the state of the simulation.
	void Render();            
	~LBMSolver();

private:
	double tau;							//relaxation time
	double c;							//lattice speed
	double sigma;						//surface tension parameter
	double p;							//pressure parameter
	double cellSize;	                //cell Size in meters
	int numCells;						//number of cells in each direction
	double dt;                          //lenght of each time step in seconds
	int current, other;					//To differenciate between the current mesh and the auxiliar mesh
	LBMCell **mesh;                     //LBM mesh
	cellTag *tags;						//Cell type for each cell in the grid.
	int index(int i, int j) { return i*numCells + j; }        //To access the mesh array as it was a 2D matrix.
	void separateIndex(int ij, int*i, int*j) {*i = ij / numCells; *j = ij - numCells*(*i);}
	double Fequi(double ux, double uy, double rho, int l);    //equilibrium function
	double MassDistribute(int type, const std::vector<int>& cells);		//Performs the mass exchange for the ifull or iempty interface cells. Returns the excess mass that can't be asigned to the neighbouring interface cells. 
	void MassExchange(int ij);          //Mass exchange process for one cell
	void SurfaceNormal(int ij);			//Update the normal vector of cell ij with the normal to the surface using the gradient of the filling, approximated by central diferences. 
	std::array<double, 2> g;			//acceleration of gravity
	const std::array<double, 9> w;		//Weights
	const std::array<int, 9> ex;        //x components of the velocity vectors
	const std::array<int, 9> ey;        //y components of the velocity vectors
	const std::array<double, 9> modE;   //Modulus of each direction vector e. 
	const std::array<int, 9> finv;      //Index of the velocity vector pointing in the opposite direction for each of the 9 velocity vectors. 
	std::list<int> changeTag;           //List of the cells that are going to change tag at the begining of the next step
	std::list<int> interfaceCells;      //List of the cells in the interface. 
	std::list<cPoint> freeSurface;		//List with points on the free surface

};

