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
#include "LBMCell.h"

LBMCell::LBMCell() :
f({ { 4. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 36., 1. / 36., 1. / 36., 1. / 36. } })
{
	feq = f;
}

void LBMCell::InitFluid()
{
	mass = rho;
	n[0] = 0.0;
	n[1] = 0.0;
	curvature = -1;
	coord[3] = 1;
}

void LBMCell::InitGas()
{
	mass = 0.0;
	u[0] = 0.0;
	u[1] = 0.0;
	n[0] = 0.0;
	n[1] = 0.0;
	curvature = -1;
	coord[3] = 0.0;
}

void LBMCell::UpdateFill()
{
	coord[3] = mass / rho;
}

