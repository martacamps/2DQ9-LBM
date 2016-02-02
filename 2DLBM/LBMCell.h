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

//Data stored in each D2Q9 Lattice Boltzmenn Method (LBM) cell


#pragma once
#include "cPoint.h"

//enum cellTag
//{
//	fluid, gas, interface, ifull, iempty, slipbc, noslipbc
//};

struct LBMCell
{
	double mass = 0;           //cell mass
	double rho = 1;            //cell density
	std::array<double, 2> u;   //cell velocity
	std::array<double, 9> f;   //particle distribution functions
	std::array<double,2>  n;   //normal vector to the free surface
	double curvature;		   //curvature of the free surface.
	//double epsilon = 0.0;	   //How filled with water is the cell
	cPoint coord;              //Coordinates of the cell in the mesh. It also stores its fill level. 
	//cellTag tag = gas;         //cell type
	// bool newInterface = false; //true if the cell will change to a interface cell in the next step

	LBMCell::LBMCell() : 
		f({ { 4. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 36., 1. / 36., 1. / 36., 1. / 36. } }),
		u({ { 0., 0. } }),
		n({ { 0., 0. } }){}

	void InitFluid();         //Initialize the cell as fluid
	void InitGas();			  //Initialize the cell as gas

}; 
