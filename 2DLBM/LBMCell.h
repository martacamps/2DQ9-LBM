#pragma once

#define FLUID 0
#define GAS 1
#define INTERFACE 2
#define IFULL 3
#define IEMPTY 4
#define SLIPBC 5
#define NOSLIPBC 6
#define FIXEDV 7
#define INTERNAL 8

struct LBMCell
{
	double mass = 1;
	double rho = 1;
	std::array<double, 2> u;   //cell velocity
	std::array<double, 9> f;
	int tag = FLUID;
	int BC = INTERNAL;
	bool newInterface = false;

	LBMCell::LBMCell() : f({ { 4. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 36., 1. / 36., 1. / 36., 1. / 36. } }), u({ { 0., 0. } }){}

}; 
