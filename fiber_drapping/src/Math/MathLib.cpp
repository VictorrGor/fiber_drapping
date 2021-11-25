// This is a personal academic project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: http://www.viva64.com
#include "Math/MathLib.h"


void LUDecomposition(double** A, size_t q, double*** L, double*** U)
{

#ifdef LOG_ON
	std::cout.width(5);
	std::cout << "A is: \n";
	for (size_t i = 0; i < q; ++i)
	{
		for (size_t j = 0; j < q; ++j)
			std::cout << std::setprecision(4) << A[i][j] << "\t";
		std::cout << "\n";
	}
#endif

	for (size_t i = 0; i < q; ++i)
	{
		(*U)[i] = (double*)calloc(q, sizeof(double));
		(*L)[i] = (double*)calloc(q, sizeof(double));
		(*U)[i][i] = 1;
	}

	for (int i = 0; i < q; ++i)
	{
		for (int j = 0; j < q; ++j)
		{
			double sum = 0;
			if (i < j)
			{
				for (int k = 0; k < i; ++k)
					sum += (*L)[i][k] * (*U)[k][j];

				if (!(*L)[i][i]) throw "LUDecomposition: division by zero!";
				(*U)[i][j] = (A[i][j] - sum) / (*L)[i][i];
			}
			else
			{
				for (int k = 0; k < j; ++k)
					sum += (*L)[i][k] * (*U)[k][j];

				(*L)[i][j] = (A[i][j] - sum);
			}
		}
	}
}

double* LUForwardBackward(double** L, double** U, double* b, size_t q)
{
	double* y = (double*)calloc(q, sizeof(double));

#ifdef LOG_ON
	std::cout.width(5);
	std::cout << "U is: \n";
	for (size_t i = 0; i < q; ++i)
	{
		for (size_t j = 0; j < q; ++j)
			std::cout << std::setprecision(4) << U[i][j] << "\t";
		std::cout << "\n";
	}
	std::cout << "L is: \n";
	for (size_t i = 0; i < q; ++i)
	{
		for (size_t j = 0; j < q; ++j)
			std::cout << std::setprecision(4) << L[i][j] << "\t";
		std::cout << "\n";
	}
#endif

	//L Forward
	for (size_t i = 0; i < q; ++i)
	{
		for (size_t j = 0; j < i; ++j)
		{
			b[i] -= L[i][j] * y[j];
		}
		if (!L[i][i])
		{
			delete[] y;
			throw "LUForwardBackward: L[i][i] is null!";
		}
		b[i] /= L[i][i];
		y[i] = b[i];
	}
	//U Backward
	double* res = (double*)calloc(q, sizeof(double));
	for (int i = q - 1; i >= 0; --i)
	{
		for (size_t j = q - 1; j > i; --j)
		{
			y[i] -= U[i][j] * res[j];
		}
		if (!U[i][i]) throw "LUForwardBackward: U[i][i] is null!";
		y[i] /= U[i][i];
		res[i] = y[i];
	}
	delete[] y;
	return res;
}
