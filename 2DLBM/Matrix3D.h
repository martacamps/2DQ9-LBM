#pragma once
template <class Type> class Matrix3D
{
public:
	Matrix3D(int width, int height, int depth) : m_width(width), m_height(height), m_depth(depth){ data = new Type[m_width*m_height*m_depth]; }
	~Matrix3D() { delete data; }
	//T* get(int x, int y, int z) {return;
	//void set(int x, int y, int z, T& d);

	size_t index(int x, int y, int z) { return x + m_width*y + m_height*z; }  //CHECK IF THIS IS OK (THE WAY THE CALCULUS IS DONE)!!!!
	int m_width, m_height, m_depth;
	Type *data;
};

