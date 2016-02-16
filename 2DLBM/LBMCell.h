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
	LBMCell();

	double mass = 0;           //cell mass
	double rho = 1;            //cell density
	cPoint u;				   //cell velocity
	std::array<double, 9> f;   //particle distribution functions
	cPoint n;				   //normal vector to the free surface
	double curvature;		   //curvature of the free surface.
	cPoint coord;              //Coordinates of the cell in the mesh. It also stores its fill level. 

	void InitFluid();         //Initialize the cell as fluid
	void InitGas();			  //Initialize the cell as gas
	void UpdateFill();        //Calculates the fill level for the cell 
}; 
