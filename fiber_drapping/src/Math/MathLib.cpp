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



void LUDecomposition(double* A, size_t q, double*** L, double*** U)
{

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
				(*U)[i][j] = (A[i*q + j] - sum) / (*L)[i][i];
			}
			else
			{
				for (int k = 0; k < j; ++k)
					sum += (*L)[i][k] * (*U)[k][j];

				(*L)[i][j] = (A[i*q + j] - sum);
			}
		}
	}
}


void LUDecomposition(double* A, size_t q, double* L, double* U)
{
	memset(L, 0, sizeof(double) * q * q);
	memset(U, 0, sizeof(double) * q * q);
	for (size_t i = 0; i < q; ++i)
	{
		U[i*q+i] = 1;
	}
	for (int i = 0; i < q; ++i)
	{
		for (int j = 0; j < q; ++j)
		{
			double sum = 0;
			if (i < j)
			{
				for (int k = 0; k < i; ++k)
					sum += L[i*q + k] * U[k * q + j];

				if (!L[i * q + i]) throw "LUDecomposition: division by zero!";
				U[i * q + j] = (A[i * q + j] - sum) / L[i * q + i];
			}
			else
			{
				for (int k = 0; k < j; ++k)
					sum += L[i * q + k] * U[k * q + j];

				L[i * q + j] = (A[i * q + j] - sum);
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

void LUForwardBackward(double* L, double* U, double* b, size_t q, double* res, double* y)
{
	memset(y, 0, sizeof(double) * q);
	memset(res, 0, sizeof(double) * q);
	//L Forward
	for (size_t i = 0; i < q; ++i)
	{
		for (size_t j = 0; j < i; ++j)
		{
			b[i] -= L[i*q + j] * y[j];
		}
		if (!L[i * q + i])
		{
			throw "LUForwardBackward: L[i][i] is null!";
		}
		b[i] /= L[i * q + i];
		y[i] = b[i];
	}
	//U Backward
	for (int i = q - 1; i >= 0; --i)
	{
		for (size_t j = q - 1; j > i; --j)
		{
			y[i] -= U[i * q + j] * res[j];
		}
		if (!U[i * q + i]) throw "LUForwardBackward: U[i][i] is null!";
		y[i] /= U[i * q + i];
		res[i] = y[i];
	}
}


double bisectionU(const surfInfo* sfI, const bSplinePt& Pij, const bSplinePt& Pim1j, double b, double eps, double A)
{
	bSplinePt ya(Pij), yb(Pim1j);
	d_vertex vx_a, vx_b, vx_midle;
	double dist_a, dist_b, dist_middle, a = Pim1j.u, d = b-a, middle;

	ya.u = a;
	yb.u = b;
	ya.pt = &vx_a;
	yb.pt = &vx_b;
	vx_a = SurfacePoint(sfI, ya.u, ya.v);
	vx_b = SurfacePoint(sfI, yb.u, yb.v);

	dist_a = getDistance(*ya.pt, *Pim1j.pt) - A;
	dist_b = getDistance(*yb.pt, *Pim1j.pt) - A;
	if (signbit(dist_a) == signbit(dist_b)) 
	{
		if (dist_b < eps) return 1.;
		throw("bisectionU singbit(dist1) == signbit(dist2)"); ///@todo If this exception was thrown, than use different method(such Newton)
	}
	do
	{
		d /= 2;
		middle = ya.u + d;
		vx_midle = SurfacePoint(sfI, middle, Pim1j.v);
		dist_middle = getDistance(vx_midle, *Pim1j.pt) - A;
		
		if (signbit(dist_middle) == signbit(dist_a)) ya.u = middle;
		else yb.u = middle;

	} while (fabs(dist_middle) > eps);

	return middle;
}

double bisectionV(const surfInfo* sfI, const bSplinePt& Pij, const bSplinePt& Pijm1, double b, double eps, double B)
{
	bSplinePt ya(Pij), yb(Pijm1);
	d_vertex vx_a, vx_b, vx_midle;
	double dist_a, dist_b, dist_middle, a = Pijm1.v, d = b - a, middle;

	ya.v = a;
	yb.v = b;
	ya.pt = &vx_a;
	yb.pt = &vx_b;
	vx_a = SurfacePoint(sfI, ya.u, ya.v);
	vx_b = SurfacePoint(sfI, yb.u, yb.v);

	dist_a = getDistance(*ya.pt, *Pijm1.pt) - B;
	dist_b = getDistance(*yb.pt, *Pijm1.pt) - B;
	if (signbit(dist_a) == signbit(dist_b))
	{
		if (dist_b < eps) return 1.;
		throw("bisectionV singbit(dist1) == signbit(dist2)"); ///@todo If this exception was thrown, than use different method(such Newton)
	}
	do
	{
		d /= 2;
		middle = ya.v + d;
		vx_midle = SurfacePoint(sfI, Pijm1.u, middle);
		dist_middle = getDistance(vx_midle, *Pijm1.pt) - B;

		if (signbit(dist_middle) == signbit(dist_a)) ya.v = middle;
		else yb.v = middle;

	} while (fabs(dist_middle) > eps);

	return middle;
}