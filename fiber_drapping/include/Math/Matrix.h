#pragma once
#include <iostream>

typedef unsigned int u_int;

template <class T>
class Matrix;


template <class T>
class Row
{
	T* pValues;
	u_int size;
public:
	friend Matrix<T>;
	Row() : pValues(nullptr), size(0) {}
	Row(u_int _size, const T _fillItem) : size(_size) 
	{
		pValues = new T[size];
		for (u_int i = 0; i < size; ++i)
		{
			pValues[i] = _fillItem;
		}
	}
	~Row()
	{
		delete[] pValues;
		pValues = nullptr;
	}
	Row(const Row& _obj)
	{
		size = _obj.size;
		if (size)
		{
			pValues = new T[size];

			for (u_int i = 0; i < size; ++i)
			{
				pValues[i] = _obj.pValues[i];
			}
		}
		else pValues = nullptr;
	}
	Row operator=(const Row<T>& _obj)
	{
		if (pValues) delete[] pValues;

		size = _obj.size;
		if (size)
		{
			pValues = new T[size];

			for (u_int i = 0; i < size; ++i)
			{
				pValues[i] = _obj.pValues[i];
			}
		}

		return *this;
	}
	Row operator+(const Row<T>& _obj1)
	{
		if (this->size != _obj1.size) throw "Row::binary_plus unequal size!\n";
		Row res(size, 0);
		for (u_int i = 0; i < size; ++i)
		{
			res.pValues[i] = pValues[i] + _obj1.pValues[i];
		}
		return res;
	}
	Row operator-(const Row<T>& _obj1)
	{
		if (size != _obj1.size) throw "Row::binary_minus unequal size!\n";
		Row res(size, 0);
		for (u_int i = 0; i < size; ++i)
		{
			res.pValues[i] = pValues[i] - _obj1.pValues[i];
		}
		return res;
	}
	Row operator*(double _multiply)
	{
		Row res(size, 0);
		for (u_int i = 0; i < size; ++i)
		{
			res.pValues[i] = pValues[i] * _multiply;
		}
		return res;
	}
	Row operator/(double _multiply)
	{
		if (!_multiply) throw "Matrix::Row::operator/: Dividing by zero!\n";
		Row res(size, 0);
		for (u_int i = 0; i < size; ++i)
		{
			res.pValues[i] = pValues[i] / _multiply;
		}
		return res;
	}
	T& operator[](u_int _colomn)
	{
		if (size <= _colomn) throw "Row::operator[] : wrong size!\n";
		return pValues[_colomn];
	}

	u_int getSize() { return size; }
};

template <class T>
class Matrix
{
	Row<T>* pRows;
	u_int size;//coloums count
public:
	Matrix() : pRows(nullptr), size(0){}
	Matrix(u_int _size, u_int _RowSize, T _fillItem) : size(_size)
	{
		if ((_size && _RowSize) == 0) throw "Matrix: size == 0";
		
		pRows = new Row<T>[size];
		
		Row<T>* pStandart = new Row<T>(_RowSize, _fillItem);
		for (u_int i = 0; i < size; ++i)
		{
			pRows[i] = *pStandart;
		}

		delete pStandart;
	}
	~Matrix()
	{
		delete[] pRows;
		pRows = nullptr;
	}
	Matrix(const Matrix& _obj)
	{
		if (_obj.size)
		{
			size = _obj.size;
			pRows = new Row<T>[size];

			for (u_int i = 0; i < size; ++i)
			{
				pRows[i] = _obj.pRows[i];
			}
		}
		else
		{
			size = 0;
			pRows = nullptr;
		}
	}
	Matrix<T> operator=(const Matrix<T> _obj)
	{
		if (_obj.size)
		{
			if (pRows) delete[] pRows;

			size = _obj.size;
			pRows = new Row<T>[size];

			for (u_int i = 0; i < size; ++i)
			{
				pRows[i] = _obj.pRows[i];
			}
		}
		return *this;
	}
	
	Matrix<T> operator+(const Matrix<T>& _obj)
	{
		if (_obj.size != size) throw "Matrix::binary plus : unequal size!\n";
		Matrix<T> res(size, pRows[0].size, 0);

		for (u_int i = 0; i < size; ++i)
		{
			for (u_int j = 0; j < pRows[0].size; ++j)
			{
				res.pRows[i].pValues[j] = pRows[i].pValues[j] + _obj.pRows[i].pValues[j];
			}
		}

		return res;
	}
	Matrix<T> operator-(const Matrix<T>& _obj)
	{
		if (_obj.size != size) throw "Matrix::binary minus : unequal size!\n";
		Matrix<T> res(size, pRows[0].size, 0);

		for (u_int i = 0; i < size; ++i)
		{
			for (u_int j = 0; j < pRows[0].size; ++j)
			{
				res.pRows[i].pValues[j] = pRows[i].pValues[j] - _obj.pRows[i].pValues[j];
			}
		}

		return res;
	}
	Matrix<T> operator*(const Matrix<T>& _obj)
	{
		if (pRows[0].size != _obj.size) throw "Matrix::binary multiply : unequal size!\n";
		u_int rowCount = size;
		u_int colomsCount = _obj.pRows[0].size;
		Matrix<T> res(rowCount, colomsCount, 0);

		for (u_int i = 0; i < rowCount; ++i)
		{
			for (u_int j = 0; j < colomsCount; ++j)
			{
				res.pRows[i].pValues[j] = 0;
				for (u_int r = 0; r < _obj.size; ++r)
				{
					res.pRows[i].pValues[j] += pRows[i].pValues[r] * _obj.pRows[r].pValues[j];
				}
			}
		}

		return res;
	}
	Matrix<T> operator*(double _s)
	{
		Matrix<T> res(*this);
		for (u_int i = 0; i < size; ++i)
		{
			for (u_int j = 0; j < pRows[0].size; ++j)
			{
				res.pRows[i].pValues[j]  = pRows[i].pValues[j] * _s;
			}
		}

		return res;
	}
	
	Row<T>& operator[](u_int _row)
	{
		if (size <= _row) throw "Matrix::operator[] : wrong size!\n";
		return pRows[_row];
	}

	void transpose()
	{
		Matrix<T> res(pRows[0].size, size, 0);

		for (u_int i = 0; i < size; ++i)
		{
			for (u_int j = 0; j < pRows[0].size; ++j)
			{
				res.pRows[j].pValues[i] = pRows[i].pValues[j];
			}
		}
		*this = res;
	}

	void print()
	{
		for (u_int i = 0; i < size; ++i)
		{
			for (u_int j = 0; j < pRows[0].size; ++j)
			{
				std::cout << pRows[i].pValues[j] << "\t";
			}
			std::cout << std::endl;
		}
	}
	void fill(const T& _item, u_int _row, u_int _colomn)
	{
		if ((_row >= size) && (pRows[0].size > _colomn)) throw "Matrix::fill() : wrong size";
		pRows[_row].pValues[_colomn] = _item;
	}

	u_int getSize() { return size; }
	double norm()
	{
		double res = 0;

		for (u_int i = 0; i < size; ++i) res = pow(pRows[i][0], 2);

		return sqrt(res);
	}

	Matrix<T> makeUpperDiagionalMx()
	{
		Matrix<T> result(*this);
		if (size != pRows->getSize()) return *this;
		else
		{
			for (int i = 1; i < size; ++i)
			{
				for (int j = 0; j < i; ++j)
				{
					result[i][j] = 0;
				}
			}
		}
		return result;
	}

	Matrix<T> makeLowerDiagionalMx()
	{
		Matrix<T> result(*this);
		if (size != pRows->getSize()) return *this;
		else
		{
			for (int i = 0; i < size - 1; ++i)
			{
				for (int j = i + 1; j < size; ++j)
				{
					result[i][j] = 0;
				}
			}
		}
		return result;
	}
};


template <class T>
Matrix<T> calculateObr( Matrix<T> _obj)
{
	Matrix<T> res(_obj.getSize(), _obj.getSize(), 0);
	
	for (int i = 0; i < res.getSize(); ++i)
	{
		res[i][i] = 1;
	}

	//int i = 2;
	for (int i = 0; i < res.getSize() ; ++i)
	{
		 //Non zero index
		res[i] = res[i] / _obj[i][i];
		_obj[i] = _obj[i] / _obj[i][i];
		for (int j = i + 1; j < res.getSize(); ++j)
		{
			res[j] = res[j] - res[j - 1] * _obj[j][i];
			_obj[j] = _obj[j] - _obj[j - 1] * _obj[j][i];
		}
	}

	for (int i = res.getSize() - 1; i > 0 ; --i)
	{
		res[i] = res[i] / _obj[i][i];
		_obj[i] = _obj[i] / _obj[i][i];
		//Non zero index
		for (int j = i - 1; j >= 0; --j)
		{
			res[j] = res[j] - res[j + 1] * _obj[j][i];
			_obj[j] = _obj[j] - _obj[j + 1] * _obj[j][i];
		}
	}
	return res;
}



//https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
template <class T>
Matrix<T> tridigionalMtx(Matrix<T>& _left, Matrix<T>& _right)
{
	Matrix<T> res(_left.getSize(), 1, 0);

	u_int dim = _right.getSize() - 1;
	double* C = new double[_left.getSize() - 1];
	double* D = new double[_left.getSize() - 1];
	memset(C, 0, sizeof(double) * _left.getSize());
	memset(D, 0, sizeof(double) * _left.getSize());
	
	//Forward
	C[0] = _left[0][1] / _left[0][0];
	D[0] = _right[0][0] / _left[0][0];
	for (int i = 1; i < dim; ++i)
	{
		C[i] = _left[i][i + 1] / (_left[i][i] - _left[i][i - 1] * C[i-1]);
		D[i] = (_right[i][0] - _left[i][i - 1] * D[i - 1]) / (_left[i][i] - _left[i][i - 1] * C[i - 1]);
	}
	//Backward
	res[dim][0] = (_right[dim][0] - _left[dim][dim - 1] * D[dim]) / (_left[dim][dim] - _left[dim][dim - 1] * C[dim - 1]);

	for (int i = dim - 1; i >= 0; i--)
	{
		res[i][0] = D[i] - C[i] * res[i + 1][0];
	}
	return res;
}