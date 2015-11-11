// 2DLBM.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include "LBMSolver.h"


int main(int argc, _TCHAR* argv[])
{
	try
	{
		LBMSolver cavity(1e-6, 7e-4, 9.81, 1000, 0.0001, 0.01, 0.5);
		cavity.init();
		cavity.mainLoop();
		
	}
	catch (std::exception& e)
	{
		std::cout << e.what() << '\n';
	}
	catch (std::string se)
	{
		std::cout << se << " Interrupting simulation" << '\n';
	}
	

	return 0;
}

