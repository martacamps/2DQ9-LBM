#pragma once
template<typename T> class cMatrix2D
{
public:
	cMatrix2D(){}
	cMatrix2D(int rows, int cols)
	{
		MatrixResize(rows, cols);
	}
	~cMatrix2D(){}

	size_t GetCols(){ return m_cols; }
	size_t GetRows(){ return m_rows; }
	std::valarray<T>& operator[](int i) { return m_values[i]; }
	std::valarray<T> operator[](int i) const { return m_values[i]; }

	void MatrixResize(int rows, int cols)          //Resizes the matrix. It will erase all its content. 
	{
		m_values.resize(rows);
		for (auto it = begin(m_values); it != end(m_values); ++it)
			it->resize(cols,0.);

		m_cols = cols;
		m_rows = rows;
	}

	cMatrix2D<T>& operator=(const cMatrix2D<T>& other)
	{
		//WHAT HAPPENS IF THE SIZES OF THE ARRAYS ARE NTO THE SAME? WILL THE EQUAL FROM VALARRAY THROW AN EXCEPTION? OR DO THIS FUNCTION HAS TO THROW IT?
		for (int i = 0; i < rows; ++i)
			m_values[i] = other[i];
		return *this;
	}
	T operator%(const cMatrix2D<T>& other) const      //Inner product
	{
		std::valarray<T> aux(m_rows);
		std::valarray<T> auxCol(m_cols);
		for (int i = 0; i < m_rows; ++i)
		{
			auxCol = m_values[i] * other[i];
			aux[i] = auxCol.sum();
		}
		return aux.sum();
	}

private:
	size_t m_cols;
	size_t m_rows;
	std::valarray<std::valarray<T>> m_values;
};


