#pragma once

#include "LBMCell.h"
#include "cSolver.h"


class LBMSolver : public cSolver
{
public:
	//LBMSolver(double nu, double sig, double g, double rho, double dx, double size, double time);
	LBMSolver();
	void Create(double nu, double sig, double dum, double rho, double dx, double size, double time);
	void InitialField();       //This initialization is only valid for the cavity flow problem for the moment.I have to allow the user to set the starting field of the fluid.
	void TimeStep(double t);
	void Render();            //Draws the state of the simulation using OpenGL with glut. 
	~LBMSolver();

private:
	double tau;							//relaxation time
	double sigma;						//surface tension parameter
	double cellSize;	                //cell Size in meters
	int numCells;						//number of cells in each direction
	double dt;                          //lenght of each time step in seconds
	double lidSpeed;                    //Horizontal Speed of the lid
	double g;							//acceleration of gravity
	int current, other;
	LBMCell **mesh;
	size_t index(int i, int j) { return i*numCells + j; }  //To access the data array as it was a 2DMatrix. Having this here instead of in another class will decrease the calls to other classes every time I want to access a position of the array. 
	const std::array<double, 9> w;
	const std::array<int, 9> ex;
	const std::array<int, 9> ey;
	const std::array<int, 9> finv;
};

