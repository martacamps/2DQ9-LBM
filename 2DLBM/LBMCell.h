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

#define FLUID 0
#define GAS 1
#define INTERFACE 2
#define NOSLIPBC 6
#define FIXEDV 7
#define INTERNAL 8

struct LBMCell
{
	float mass = 1;           //cell mass
	float rho = 1;            //cell density
	std::array<float, 2> u;   //cell velocity
	std::array<float, 9> f;   //particle distribution functions
	int tag = FLUID;           //cell type
	int BC = INTERNAL;         //FIXEDV if the cell velocity is fixed
	bool newInterface = false; //true if the cell will change to a interface cell in the next step

	LBMCell() : f({ { 4.f / 9.f, 1.f / 9.f, 1.f / 9.f, 1.f / 9.f, 1.f / 9.f, 1.f / 36.f, 1.f / 36.f, 1.f / 36.f, 1.f / 36.f } }), u({ { 0.f, 0.f } }){}

}; 
