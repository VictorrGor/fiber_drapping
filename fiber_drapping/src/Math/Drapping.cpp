﻿// This is a personal academic project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: http://www.viva64.com
#include "../Math/Drapping.h"

void getF(double* f, const d_vertex* ptIJ, const d_vertex* ptIm1J, const d_vertex* ptIJm1)
{
	f[0] = -(pow((ptIJ->x - ptIm1J->x), 2) + pow((ptIJ->y - ptIm1J->y), 2) + pow((ptIJ->z - ptIm1J->z), 2) - pow(A, 2));
	f[1] = -(pow((ptIJ->x - ptIJm1->x), 2) + pow((ptIJ->y - ptIJm1->y), 2) + pow((ptIJ->z - ptIJm1->z), 2) - pow(B, 2));
}

bool getBSplineDrapPoint(double** W, double** invW, bSplinePt* ptIJ, bSplinePt* ptIm1J, bSplinePt* ptIJm1, surfInfo* sfI)
{
	const size_t maxIterCt = 10000;
	bool flag = true;
	size_t i = 0;

	size_t dim = 2;
	double* f = new double[dim];
	double* dx;

	double epsilon = 0.00001; // Варьировать сравнение ошибки в зависимости от накопленной ошибки. Т.е. если накопилась большая ошибка, то уменьшать макисмально взоможную, или прижимать в другую сторону
	double len = 100000;
	
	bool corrupted = false;

	int counter = 0;

	double** L = new double* [2];
	double** U = new double* [2];
	while (flag && (maxIterCt > i))
	{
		counter++;

		corrupted = false;
		//(*ptIJ->pt) = SurfacePoint(sfI, ptIJ->u, ptIJ->v);

		d_vertex** IJder = SurfaceDerivsAlg1(sfI, ptIJ->u, ptIJ->v, 1);
		getJakobain(W, ptIJ, ptIm1J, ptIJm1, IJder);
		getF(f, ptIJ->pt, ptIm1J->pt, ptIJm1->pt);///@todo Optimize memory using. At each iteration memory allocation happens. It too much heavy
		
		//std::cout << "\n f is: " << f[0] << "\n" << f[1] << "\n";

		LUDecomposition(W, dim, &L, &U);
		dx = LUForwardBackward(L, U, f, 2);
		dx[0] *= 0.1;
		dx[1] *= 0.1;
		//if ((ptIJ->u + dx[0] < 0) || (ptIJ->v + dx[1] < 0) || (ptIJ->u + dx[0] > 1) || (ptIJ->v + dx[1] > 1)) corrupted = true;
		d_vertex oldPoint = *(ptIJ->pt);
		
		if (ptIJ->v + dx[1] < 0)
		{
			dx[1] = 0;
			ptIJ->v = 0;
		}
		else
			if (ptIJ->v + dx[1] > 1)
			{
				dx[1] = 0;
				ptIJ->v = 1;
			}
		if (ptIJ->u + dx[0] < 0)
		{
			dx[0] = 0;
			ptIJ->u = 0;
		}
		else
			if (ptIJ->u + dx[0] > 1)
			{
				dx[0] = 0;
				ptIJ->u = 1;
			}

			
		ptIJ->u += dx[0];
		ptIJ->v += dx[1];
		(*ptIJ->pt) = SurfacePoint(sfI, ptIJ->u, ptIJ->v);
		

		for (size_t ct1 = 0; ct1 < dim; ++ct1)
		{
			delete[] L[ct1];
			delete[] U[ct1];
			delete[] IJder[ct1];
		}
		delete[] dx;
		delete[] IJder;

		if ((f[0] <= A * epsilon) && (f[1] <= B * epsilon) && !corrupted)
		//if ((fabs(f[0]) < 0.01) && (fabs(f[1]) < 0.01))
		{
			flag = false;
			break;
		}
		++i;

	}
	delete[] L;
	delete[] U;

	delete[] f;

	return !flag;
}

void getJakobain(double** W, const bSplinePt* ptIJ, const bSplinePt* ptIm1J, const bSplinePt* ptIJm1, const d_vertex* const* IJder)
{	W[0][0] = 2 * (ptIJ->pt->x - ptIm1J->pt->x) * IJder[1][0].x + 2 * (ptIJ->pt->y - ptIm1J->pt->y) * IJder[1][0].y +
		+2 * (ptIJ->pt->z - ptIm1J->pt->z) * IJder[1][0].z;
	W[0][1] = 2 * (ptIJ->pt->x - ptIm1J->pt->x) * IJder[0][1].x + 2 * (ptIJ->pt->y - ptIm1J->pt->y) * IJder[0][1].y +
		+2 * (ptIJ->pt->z - ptIm1J->pt->z) * IJder[0][1].z;

	W[1][0] = 2 * (ptIJ->pt->x - ptIJm1->pt->x) * IJder[1][0].x + 2 * (ptIJ->pt->y - ptIJm1->pt->y) * IJder[1][0].y +
		+2 * (ptIJ->pt->z - ptIJm1->pt->z) * IJder[1][0].z;
	W[1][1] = 2 * (ptIJ->pt->x - ptIJm1->pt->x) * IJder[0][1].x + 2 * (ptIJ->pt->y - ptIJm1->pt->y) * IJder[0][1].y +
		+2 * (ptIJ->pt->z - ptIJm1->pt->z) * IJder[0][1].z;
}

double getAngle(const d_vertex* const* gird, size_t i, size_t j, size_t p, size_t q, size_t s, size_t h)
{
	struct { float fi; float teta;} p0, p1, p2;
	XMFLOAT3 v1, v2;
	v1.x = gird[p][q].x - gird[i][j].x;
	v1.y = gird[p][q].y - gird[i][j].y;
	v1.z = gird[p][q].z - gird[i][j].z;


	v2.x = gird[s][h].x - gird[i][j].x;
	v2.y = gird[s][h].y - gird[i][j].y;
	v2.z = gird[s][h].z - gird[i][j].z;
	
	XMVECTOR angle = XMVector3AngleBetweenVectors(XMLoadFloat3(&v1), XMLoadFloat3(&v2));
	return angle.m128_f32[0] * 57.2958;
}

void drapping_part(RenderSys* _rs, const drappingInit& _is)
{
	UINT size = GIRD_SIZE;
	bSplinePt** P = new bSplinePt * [size]; //points warper

	std::cout << "A is: " << A << "; B is: " << B << ";\n";
	d_vertex** Q = new d_vertex * [size];
	for (UINT i = 0; i < size; ++i)
	{
		P[i] = new bSplinePt[size];

		Q[i] = new d_vertex[size];
		for (UINT j = 0; j < size; ++j)
		{
			P[i][j].pt = &Q[i][j];
		}
	}

	UINT dim = 2;
	double** W = new double* [dim]; //Jakobian
	double** invW = new double* [dim]; //inverse
	for (UINT i = 0; i < dim; ++i)
	{
		W[i] = new double[dim];
		invW[i] = new double[dim];
	}
	UINT err_ct = 0;

	double cycle_step = 1. / (size - 1);

	//Generating initial lines
	for (UINT i = 0; i < size; ++i)
	{
		if (_is.isU1)
		{
			P[0][i].u = _is.u1;
			P[0][i].v = i * cycle_step;
		}
		else
		{
			P[0][i].u = i * cycle_step;
			P[0][i].v = _is.v1;
		}
		if (_is.isU2)
		{
			P[i][0].u = _is.u2;
			P[i][0].v = i * cycle_step;
		}
		else
		{
			P[i][0].u = i * cycle_step;
			P[i][0].v = _is.v2;
		}
		Q[0][i] = SurfacePoint(_is.sfI, P[0][i].u, P[0][i].v);
		Q[i][0] = SurfacePoint(_is.sfI, P[i][0].u, P[i][0].v);
	}

	double delta_u = 0.01;

	for (UINT i = 0; i < size - 1; ++i)
	{
		for (UINT j = 0; j < size - 1; ++j)
		{
			if ((P[i][j].u < 0) || (P[i][j + 1].u < 0) || (P[i + 1][j].u < 0))
				continue;

			bSplinePt* ptIJ = &P[i + 1][j + 1];
			bSplinePt* ptIm1J = &P[i][j + 1];
			bSplinePt* ptIJm1 = &P[i + 1][j];

			ptIJ->u = (ptIJm1->u + ptIm1J->u) / 2;//ptIJm1->u;
			ptIJ->v = max(ptIJm1->v, ptIm1J->v);// ptIJm1->v - 2* delta_u;
			if ((ptIJ->u == ptIJm1->u) && (ptIJ->v == ptIJm1->v) || (ptIJ->u == ptIm1J->u) && (ptIJ->v == ptIm1J->v))
			{
				ptIJ->u = (ptIJm1->u + ptIm1J->u) / 2;
			}
			if (ptIJ->u > 1) ptIJ->u -= 2 * delta_u;
			if (ptIJ->v > 1) ptIJ->v -= 2 * delta_u;
			if (ptIJ->u < 0) ptIJ->u = (ptIJm1->u + ptIm1J->u) / 2;//0;
			if (ptIJ->v < 0) ptIJ->v = min(ptIJm1->v, ptIm1J->v);// 0;// 1 + ptIJ->v;

			(*ptIJ->pt) = SurfacePoint(_is.sfI, ptIJ->u, ptIJ->v);

#ifdef _DEBUG
			std::cout << "\ti:" << i << "; j:" << j << "\n";
#endif

			
			if (!getBSplineDrapPoint(W, invW, ptIJ, ptIm1J, ptIJm1, _is.sfI))
			{
#ifdef _DEBUG
				std::cout << "\tu:" << ptIJ->u << "; v:" << ptIJ->v << "\n";
				std::cout << "\tx:" << ptIJ->pt->x << "; y:" << ptIJ->pt->y << "; z:" << ptIJ->pt->z << "\n";
				std::cout << "Distance XJ, X_1J: " << getDistance(*(ptIJ->pt), *(ptIm1J->pt)) << "\n";
				std::cout << "Distance XJ, XJ_1: " << getDistance(*(ptIJ->pt), *(ptIJm1->pt)) << "\n";
#endif
				ptIJ->u = -1;
				ptIJ->v = -1;
				++err_ct;
				std::cout << err_ct << "\n";
				std::cout << "\ti:" << i << "; j:" << j << "\n";
			}
			else
			{ 
#ifdef _DEBUG
				std::cout << "\tu:" << ptIJ->u << "; v:" << ptIJ->v << "\n";
				std::cout << "\tx:" << ptIJ->pt->x << "; y:" << ptIJ->pt->y << "; z:" << ptIJ->pt->z << "\n";
				std::cout << "Distance XJ, X_1J: " << getDistance(*(ptIJ->pt), *(ptIm1J->pt)) << "\n";
				std::cout << "Distance XJ, XJ_1: " << getDistance(*(ptIJ->pt), *(ptIJm1->pt)) << "\n";
#endif
			}
		}
	}
	std::cout << "\nDrapping errors: " << err_ct << "\n";


	vertex* triangle = new vertex[3];
	for (UINT i = 0; i < size - 1; ++i)
	{
		for (UINT j = 0; j < size - 1; ++j)
		{
			if ((P[i][j].u < 0) || (P[i][j + 1].u < 0) || (P[i + 1][j].u < 0))
				continue;
			float angle = getAngle(Q, i, j, i + 1, j, i, j + 1);
			float red, green, blue, coeff;
			coeff = 1;
			if (angle < 90)
				red = 1 - angle * coeff / 90;
			else
				red = 0;
			if (angle > 90)
			{
				green = 1 - (angle - 90) * coeff / 90;
				blue = (angle - 90) * coeff / 90;
			}
			else
			{
				green = angle * coeff / 90;
				blue = 0;
			}

			triangle[0] = Q[i][j];
			triangle[1] = Q[i][j + 1];
			triangle[2] = Q[i + 1][j];
			triangle[0].normal = calculateTriangleNormal(Q[i][j], Q[i][j + 1], Q[i + 1][j]);
			triangle[1].normal = triangle[0].normal;
			triangle[2].normal = triangle[0].normal;

			triangle[0].Color = vec4(red, green, blue, 1);
			triangle[1].Color = vec4(red, green, blue, 1);
			triangle[2].Color = vec4(red, green, blue, 1);
			_rs->drawTriangle(triangle);
		}
	}
	for (UINT i = 0; i < size - 1; ++i)
	{
		for (UINT j = 0; j < size - 1; ++j)
		{
			float angle = getAngle(Q, i + 1, j + 1, i + 1, j, i, j + 1);
			float red, green, blue, coeff;
			coeff = 1;
			if (angle < 90)
				red = 1 - angle * coeff / 90;
			else
				red = 0;
			if (angle > 90)
			{
				green = 1 - (angle - 90) * coeff / 90;
				blue = (angle - 90) * coeff / 90;
			}
			else
			{
				green = angle * coeff / 90;
				blue = 0;
			}


			vertex* triangle = new vertex[3];
			if ((P[i + 1][j + 1].u < 0) || (P[i][j + 1].u < 0) || (P[i + 1][j].u < 0))
				continue;
			triangle[0] = Q[i + 1][j];
			triangle[1] = Q[i][j + 1];
			triangle[2] = Q[i + 1][j + 1];
			triangle[0].normal = calculateTriangleNormal(Q[i + 1][j], Q[i][j + 1], Q[i + 1][j + 1]);
			triangle[1].normal = triangle[0].normal;
			triangle[2].normal = triangle[0].normal;

			triangle[0].Color = vec4(red, green, blue, 1);
			triangle[1].Color = vec4(red, green, blue, 1);
			triangle[2].Color = vec4(red, green, blue, 1);
			_rs->drawTriangle(triangle);
		}
	}


	for (UINT i = 0; i < size; ++i)
	{
		delete[] Q[i];
		delete[] P[i];
	}
	for (UINT i = 0; i < dim; ++i)
	{
		delete W[i];
		delete invW[i];
	}
	delete[] Q;
	delete[] P;
	delete[] W;
	delete[] invW;
}
