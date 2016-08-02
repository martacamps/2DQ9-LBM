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
#include "LBMSolver.h"


LBMSolver::LBMSolver() :
m_w({ { 4.f / 9.f, 1.f / 9.f, 1.f / 9.f, 1.f / 9.f, 1.f / 9.f, 1.f / 36.f, 1.f / 36.f, 1.f / 36.f, 1.f / 36.f } }),
m_ex({ { 0, 1, -1, 0, 0, 1, -1, 1, -1 } }),
m_ey({ { 0, 0, 0, 1, -1, -1, -1, 1, 1 } }),
m_exMod({ { 0, 1, -1, 0, 0, 1, -1, 1, -1 } }),
m_eyMod({ { 0, 0, 0, 1, -1, -1, -1, 1, 1 } }),
m_finv({ { 0, 2, 1, 4, 3, 8, 7, 6, 5 } }),
m_current(0),
m_other(1)
{

}

void LBMSolver::Create(float nu, float sig, float v, float rho, float dx, float size, float time)
{
	m_cellSize = dx;
	m_lidSpeed = v;
	
	//Compute necessary information from input parameters
	m_dt = sqrt((1e-4*m_cellSize) / 9.81);    // So the gravity acceleration in the model is not larger than 1e-4.
	m_g = -9.81*(m_dt*m_dt) / m_cellSize;
	m_c = m_cellSize/m_dt;
	m_exMod = m_c*m_exMod;
	m_eyMod = m_c*m_eyMod;

	float nuStar = nu*(m_dt / (m_cellSize*m_cellSize));
	m_tau = (6 * nuStar + 1) / 2;
	if (m_tau < 0.5 || m_tau > 2.5)
	{
		std::ostringstream strs;
		strs << "Relaxation time = " << m_tau << " out of range, simulation will be too unstable";
		std::string str = strs.str();
		throw(str);
	}
	float dm = rho*m_cellSize*m_cellSize*m_cellSize;
	m_sigma = sig*(m_dt*m_dt) / dm;
	m_numCells = (int)std::round(size / m_cellSize);
	m_numSteps = (int)std::round(time / m_dt);

	//Create LBM mesh. 
	m_mesh = new LBMCell*[2];
	m_mesh[m_current] = new LBMCell[m_numCells*m_numCells];
	m_mesh[m_other] = new LBMCell[m_numCells*m_numCells];

}

void LBMSolver::InitialField()  
{
	//Set the cell tags for the Boundary Conditions.
	for (int j = 0; j < m_numCells; j++)
	{
		m_mesh[m_current][index(0, j)].tag = NOSLIPBC;
		m_mesh[m_current][index(m_numCells - 1, j)].tag = NOSLIPBC;
		
		m_mesh[m_other][index(0, j)].tag = NOSLIPBC;
		m_mesh[m_other][index(m_numCells - 1, j)].tag = NOSLIPBC;
	}
	for (int i = 1; i < (m_numCells - 1); i++)
	{

		m_mesh[m_current][index(i, 0)].tag = NOSLIPBC;
		m_mesh[m_current][index(i, m_numCells - 1)].tag = NOSLIPBC;
		m_mesh[m_other][index(i, 0)].tag = NOSLIPBC;
		m_mesh[m_other][index(i, m_numCells - 1)].tag = NOSLIPBC;
	}
	for (int j = 1; j < (m_numCells - 1); j++)
	{
		m_mesh[m_current][index(m_numCells - 2, j)].BC = FIXEDV;
		m_mesh[m_other][index(m_numCells - 2, j)].BC = FIXEDV;
	}
}

void LBMSolver::TimeStep(float t)
{
	
	for (int i = 1; i < m_numCells-1; i++)
	{
		for (int j = 1; j < m_numCells-1; j++)
		{
			int ij = index(i, j);

			//Streaming. 
			for (int l = 0; l < m_finv.size(); l++)
			{
				int inv = m_finv[l];
				int previous = index(i + m_ey[inv], j + m_ex[inv]);

				//check the type of the neighbour cell and perform streaming
				if (m_mesh[m_other][previous].tag == NOSLIPBC)  //The neighbour cell is a non slip wall
				{
					m_mesh[m_current][ij].f[l] = m_mesh[m_current][ij].f[inv];
				}
				else										//The neighbour cell is fluid
				{
					m_mesh[m_current][ij].f[l] = m_mesh[m_other][previous].f[l];
				}     
			}

			//Collision
			float rho = 0.0, ux = 0.0, uy = 0.0;
			if (m_mesh[m_current][ij].BC == FIXEDV)         //Fix the density and velocity for the cells marked as FIXEDV (the lid cells)
			{
				rho = 1.; ux = m_lidSpeed; uy = 0.0;
			}
			else                                       //calculate density and velocity 
			{
				float fi;
				for (int l = 0; l < m_finv.size(); l++)
				{
					fi = m_mesh[m_current][ij].f[l];
					rho += fi;
					ux += fi*m_exMod[l];
					uy += fi*m_eyMod[l];
				}
				ux = ux/rho;
				uy = uy/rho;
			}
			m_mesh[m_current][ij].rho = rho;
			m_mesh[m_current][ij].u[0] = ux;
			m_mesh[m_current][ij].u[1] = uy;

			//Collision: apply gravity
			if (m_mesh[m_current][ij].BC != FIXEDV)
			{
				uy += m_tau*m_g;
			}
		
			//Collision
			for (int l = 0; l < m_finv.size(); l++)
			{
				float eDotu = m_exMod[l] * ux + m_eyMod[l] * uy;
				float f0 = m_w[l] *rho*(1 - (3.*(ux*ux + uy*uy)) / (2.*m_c*m_c) + (3.*eDotu)/(m_c*m_c)  + (9.*eDotu*eDotu) / (2.*m_c*m_c*m_c*m_c));  
				//float f0 = m_w[l] * (rho - (3.*(ux*ux + uy*uy)) / (2.*m_c*m_c) + (3.*eDotu) / (m_c*m_c) + (9.*eDotu*eDotu) / (2.*m_c*m_c*m_c*m_c));
				m_mesh[m_current][ij].f[l] = m_mesh[m_current][ij].f[l] - (1 / m_tau)*(m_mesh[m_current][ij].f[l] - f0);
			}

			//Collision: Calculate velocities for display and particle tracing. 
			if (m_mesh[m_current][ij].BC != FIXEDV)
			{
				m_mesh[m_current][ij].u[1] += 0.5*m_tau*m_g;
			}
		}
	}
	//Swap m_meshes
	m_other = m_current;
	m_current = 1 - m_other;
}

void LBMSolver::Render()
{
	float win_width = glutGet(GLUT_WINDOW_WIDTH);
	float win_height = glutGet(GLUT_WINDOW_HEIGHT);

	glClear(GL_COLOR_BUFFER_BIT);
	glLoadIdentity();

	//Set colours and intervals for the selected plot. 
	std::vector<float> intervals, colours;
	std::string title;
	if (m_vis[2])  //Colours by velocity
	{
		for (int i = 0; i < 6; i++)
			intervals.push_back(i*m_lidSpeed / 5.0);
		m_colourScale.SetIntervals(&intervals);

		if (m_vis[1])   //Colour pallette 1
			colours = { 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, };   //blue, cyan, green, yellow, red, magenta
		else          //Colour pallette 2
			colours = { 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1 };   //magenta, red, yellow, green, cyan, blue
		m_colourScale.SetColours(&colours, 6);

		title = "speed[m/s]";
	}
	else       //Colours by cell tag
	{
		intervals = { FLUID, NOSLIPBC };
		m_colourScale.SetIntervals(&intervals);
		colours = { 0.2f, 0.6f, 1.0f, 0.4f, 0.2f, 0.f };   //light blue, brown.
		m_colourScale.SetColours(&colours, 2);
		title = "cell type";
	}

	//Draw and colour the cells of the domain.
	glPushMatrix();
	glScalef((0.9*win_width) / m_numCells, win_height / m_numCells, 1.);
	glShadeModel(GL_FLAT);
	std::vector<float> ux, uy, mod;
	float colourx, coloury, colourz;
	bool inside;
	for (int i = 0; i < m_numCells; i++)
	{
		glBegin(GL_TRIANGLE_STRIP);
		for (int j = 0; j < m_numCells; j++)
		{
			//Extract velocity. The velocity has not to be divided by c for the lid. 
			if (m_mesh[m_other][index(i, j)].BC == FIXEDV)
			{
				ux.push_back(m_mesh[m_other][index(i, j)].u[0]);
				uy.push_back(m_mesh[m_other][index(i, j)].u[1]);
			}
			else
			{
				ux.push_back(m_mesh[m_other][index(i, j)].u[0] / m_c);
				uy.push_back(m_mesh[m_other][index(i, j)].u[1] / m_c);
			}
			mod.push_back(sqrt(ux.back()*ux.back() + uy.back()*uy.back()));
			if (m_vis[2])	 //Colours by velocity
			{
				inside = m_colourScale.PickColour(mod.back(), &colourx, &coloury, &colourz);
			}
			else	//Colours by cell tag
			{
				inside = m_colourScale.PickColour((float)m_mesh[m_other][index(i, j)].tag, &colourx, &coloury, &colourz);
			}
			if (!inside)
			{
				colourx = 0.0; coloury = 0.0; colourz = 0.0;
			}
			glColor3f(colourx, coloury, colourz);
			glVertex2f(j, i + 1);
			glVertex2f(j, i);
			glVertex2f(j + 1, i + 1);
			glVertex2f(j + 1, i);

		}
		glEnd();
	}


	//Draw cell velocity vectors
	if (m_vis[0])
	{
		glColor3f(0., 0.0, 0.);
		for (int i = 0; i < m_numCells; i++)
		{
			for (int j = 0; j < m_numCells; j++)
			{
				glBegin(GL_LINE_STRIP);
				glVertex2f(1 + j, 0.5 + i);
				glVertex2f(1 + j + ux[index(i, j)] / mod[index(i, j)], 0.5 + i + uy[index(i, j)] / mod[index(i, j)]);
				glVertex2f(0.8 + j + ux[index(i, j)] / mod[index(i, j)], 0.8 + i + uy[index(i, j)] / mod[index(i, j)]);   //Arrow head
				glEnd();
			}
		}
	}

	glPopMatrix();

	//Draw colour scale 
	glShadeModel(GL_SMOOTH);
	glPushMatrix();
	glTranslatef(0.9*win_width, 0., 0.);
	glScalef(0.1*win_width, win_height, 1.);
	m_colourScale.DrawScale(title);
	glPopMatrix();

	//Draw time counter
	std::ostringstream strs;
	strs << std::setprecision(2) << m_t*m_dt << "s";
	std::string time = strs.str();
	glColor3f(1., 1., 1.);
	m_colourScale.DrawText<std::string>(GLUT_BITMAP_HELVETICA_18, 0.9*win_width, 0.96*win_height, &time);

	glutSwapBuffers();


}

LBMSolver::~LBMSolver()
{
	if (m_mesh)
	{
		delete m_mesh[0];
		delete m_mesh[1];
		delete m_mesh;
	}

}
