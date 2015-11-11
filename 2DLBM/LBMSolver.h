#pragma once

#include "LBMCell.h"
//#include "Matrix3D.h"

class LBMSolver
{
public:
	LBMSolver(double nu, double sig, double g, double rho, double dx, double size, double time);
	void init();       //This initialization is only valid for the cavity flow problem for the moment.I have to allow the user to set the starting field of the fluid.
	void mainLoop();
	~LBMSolver();

private:
	double tau;							//relaxation time
	double sigma;						//surface tension parameter
	double f;							//body force acceleration
	double cellSize;	
	int numCells;						//number of cells in each direction
	double dt;
	double numSteps;  
	int current, other;
	//Matrix3D<LBMCell> *mesh;
	LBMCell **mesh;
	size_t index(int i, int j) { return i*numCells + j; }  //To access the data array as it was a 3DMatrix. Having this here instead of in another class will decrease the calls to other classes every time I want to access a position of the array. CHECK IF THIS IS OK (THE WAY THE CALCULUS IS DONE)!!!! 
	const std::array<double, 9> w;
	const std::array<int, 9> ex;
	const std::array<int, 9> ey;
	const std::array<int, 9> finv;
};

