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


class LBMSolver : public cSolver
{
public:
	LBMSolver();

	//Input the solver parameters and prepare the solver to start
	void Create(float nu, float sig, float dum, float rho, float dx, float size, float time);
	//Set the initial conditions for a lid driven cavity flow
	void InitialField();       
	//2DQ9 LBM time step
	void TimeStep(float t);
	//Draws the state of the simulation.
	void Render();            
	~LBMSolver();

private:
	float m_tau;							//relaxation time
	float m_c;							//lattice speed
	float m_sigma;						//surface tension parameter
	float m_cellSize;	                //cell Size in meters
	int m_numCells;						//number of cells in each direction
	float m_dt;                          //lenght of each time step in seconds
	float m_lidSpeed;                    //Horizontal Speed of the lid
	float m_g;							//acceleration of gravity
	int m_current, m_other;					//To differenciate between the current mesh and the auxiliar mesh
	LBMCell **m_mesh;                     //LBM mesh
	int index(int i, int j) { return i*m_numCells + j; }  //To access the mesh array as it was a 2D matrix.
	const std::array<float, 9> m_w;		//Weights
	const std::array<int, 9> m_ex;        //x direction of the velocity vectors
	const std::array<int, 9> m_ey;        //y direction of the velocity vectors
	std::valarray<float> m_exMod;        //x components of the velocity vectors
	std::valarray<float> m_eyMod;        //y components of the velocity vectors
	const std::array<int, 9> m_finv;      //Index of the velocity vector pointing in the opposite direction for each of the 9 velocity vectors. 
};

