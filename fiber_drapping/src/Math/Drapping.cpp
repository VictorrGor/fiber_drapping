// This is a personal academic project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: http://www.viva64.com
#include "../Math/Drapping.h"

void getF(double* f, const d_vertex* ptIJ, const d_vertex* ptIm1J, const d_vertex* ptIJm1, double A, double B)
{
	f[0] = -(pow((ptIJ->x - ptIm1J->x), 2) + pow((ptIJ->y - ptIm1J->y), 2) + pow((ptIJ->z - ptIm1J->z), 2) - pow(B, 2));
	f[1] = -(pow((ptIJ->x - ptIJm1->x), 2) + pow((ptIJ->y - ptIJm1->y), 2) + pow((ptIJ->z - ptIJm1->z), 2) - pow(A, 2));
}

bool getBSplineDrapPoint(double** W, double** invW, drappingCell& cell)
{
	const size_t maxIterCt = 10000;
	bool flag = true;
	size_t i = 0;

	size_t dim = 2;
	double* f = new double[dim];
	double* dx;

	double epsilon = 0.0001; // Варьировать сравнение ошибки в зависимости от накопленной ошибки. Т.е. если накопилась большая ошибка, то уменьшать макисмально взоможную, или прижимать в другую сторону
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

		d_vertex** IJder = SurfaceDerivsAlg1(cell.si->sfI, cell.ptIJ->u, cell.ptIJ->v, 1);
		getJakobain(W, cell.ptIJ, cell.ptIm1J, cell.ptIJm1, IJder);
		getF(f, cell.ptIJ->pt, cell.ptIm1J->pt, cell.ptIJm1->pt, cell.si->A, cell.si->B);///@todo Optimize memory using. At each iteration memory allocation happens. It too much heavy
		
		//std::cout << "\n f is: " << f[0] << "\n" << f[1] << "\n";

		LUDecomposition(W, dim, &L, &U);
		dx = LUForwardBackward(L, U, f, 2);
		dx[0] *= 0.01;
		dx[1] *= 0.01;
		//if ((ptIJ->u + dx[0] < 0) || (ptIJ->v + dx[1] < 0) || (ptIJ->u + dx[0] > 1) || (ptIJ->v + dx[1] > 1)) corrupted = true;
		d_vertex oldPoint = *(cell.ptIJ->pt);
		
		if (cell.ptIJ->v + dx[1] < 0)
		{
			dx[1] = 0;
			cell.ptIJ->v = 0;
		}
		else
			if (cell.ptIJ->v + dx[1] > 1)
			{
				dx[1] = 0;
				cell.ptIJ->v = 1;
			}
		if (cell.ptIJ->u + dx[0] < 0)
		{
			dx[0] = 0;
			cell.ptIJ->u = 0;
		}
		else
			if (cell.ptIJ->u + dx[0] > 1)
			{
				dx[0] = 0;
				cell.ptIJ->u = 1;
			}

			
		cell.ptIJ->u += dx[0];
		cell.ptIJ->v += dx[1];
		(*cell.ptIJ->pt) = SurfacePoint(cell.si->sfI, cell.ptIJ->u, cell.ptIJ->v);
		

		for (size_t ct1 = 0; ct1 < dim; ++ct1)
		{
			delete[] L[ct1];
			delete[] U[ct1];
			delete[] IJder[ct1];
		}
		delete[] dx;
		delete[] IJder;

		double bufA = (cell.si->A) * epsilon - cell.accumULen;
		double bufB = (cell.si->B) * epsilon - cell.accumVLen;

		//if ((f[0] <= (cell.si->A) * epsilon - cell.accumULen) && (f[1] <= (cell.si->B) * epsilon - cell.accumVLen) && !corrupted)
		//if ((fabs(f[0]) < 0.01) && (fabs(f[1]) < 0.01))
		if ((fabs(getDistance(*cell.ptIJ->pt, *cell.ptIJm1->pt) - cell.si->A) < cell.si->A * epsilon) 
			&& (fabs(getDistance(*cell.ptIJ->pt, *cell.ptIm1J->pt) - cell.si->B) < cell.si->B* epsilon))
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


void generateInitialLines(bSplinePt** P, d_vertex** Q, const drappingInit& _is)
{

	double cycle_step = 1. / (_is.gird_size - 1);
	std::cout << "Generate initial lines:\n";

	P[0][0].u = 0;
	P[0][0].v = 0;
	
	if (_is.isU1) P[0][0].u = _is.u1;
	else P[0][0].v = _is.v1;
	if (_is.isU2) P[0][0].u = _is.u2;
	else P[0][0].v = _is.v2;

	Q[0][0] = SurfacePoint(_is.sfI, P[0][0].u, P[0][0].v);
	double eps = 0.0001;
	for (UINT i = 1; i < _is.gird_size; ++i)
	{
		if (_is.isU1)
		{
			P[0][i].u = _is.u1;
			P[0][i].v = i * cycle_step;
			P[0][i].v = bisectionV(_is.sfI, P[0][i], P[0][i - 1], 1., eps*_is.B, _is.B);
		}
		else
		{
			P[0][i].u = i * cycle_step;
			P[0][i].v = _is.v1;
			P[0][i].u = bisectionU(_is.sfI, P[0][i], P[0][i - 1], 1., eps * _is.A, _is.A);
		}
		if (_is.isU2)
		{
			P[i][0].u = _is.u2;
			P[i][0].v = i * cycle_step;
			P[i][0].v = bisectionV(_is.sfI, P[i][0], P[i-1][0], 1., eps * _is.B, _is.B);
		}
		else
		{
			P[i][0].u = i * cycle_step;
			P[i][0].v = _is.v2;
			P[i][0].u = bisectionU(_is.sfI, P[i][0], P[i - 1][0], 1., eps * _is.A, _is.A);
		}
		Q[i][0] = SurfacePoint(_is.sfI, P[i][0].u, P[i][0].v);
		Q[0][i] = SurfacePoint(_is.sfI, P[0][i].u, P[0][i].v);
#ifdef _DEBUG
		std::cout << "i: " << i << ";\n\t(Q[0][i];Q[0][i - 1]) dist is: " << getDistance(Q[0][i], Q[0][i - 1]) << "\n";
		std::cout << "\t(Q[i][0];Q[i-1][0]) dist is: " << getDistance(Q[i][0], Q[i-1][0]) << "\n";
#endif
	}
}

void drawDrappedCell(RenderSys* _rs, bSplinePt** P, d_vertex** Q, size_t size)
{
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
}

void makeDrappedGird(RenderSys* _rs, const drappingInit& _is)
{
	std::cout << "A is: " << _is.A << "; B is: " << _is.B << ";\n";
	
	UINT size = _is.gird_size;
	bSplinePt** P = new bSplinePt * [size]; //points warper
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
	double delta_u = 0.01;
	drappingCell cell = { nullptr, nullptr, nullptr, &_is };
	//Length accumalation array by u and v coordinate;
	size_t accum_size = (size) * (size);
	double* uAcc = new double[accum_size];
	double* vAcc = new double[accum_size];
	memset(uAcc, 0, sizeof(double) * accum_size);
	memset(vAcc, 0, sizeof(double) * accum_size);

	{
		TimeBench tb;
		generateInitialLines(P, Q, _is);
	}
	{
		TimeBench tb;
		for (UINT i = 0; i < size - 1; ++i)
		{
			for (UINT j = 0; j < size - 1; ++j)
			{
				if ((P[i][j].u < 0) || (P[i][j + 1].u < 0) || (P[i + 1][j].u < 0))
					continue;

				bSplinePt* ptIJ = &P[i + 1][j + 1];
				bSplinePt* ptIm1J = &P[i][j + 1];
				bSplinePt* ptIJm1 = &P[i + 1][j];

				cell.ptIJ = ptIJ;
				cell.ptIJm1 = ptIJm1;
				cell.ptIm1J = ptIm1J;
				cell.accumULen = uAcc[(i + 1) * size + j];
				cell.accumVLen = vAcc[i * size + j + 1];

				//ptIJ->u = max(ptIJm1->u, ptIm1J->u);
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
				if (!getBSplineDrapPoint(W, invW, cell))
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
					uAcc[(i + 1) * size + (j + 1)] = uAcc[(i + 1) * size + j] + getDistance(*(ptIJ->pt), *(ptIJm1->pt)) - _is.A;
					vAcc[(i + 1) * size + (j + 1)] = vAcc[i * size + j + 1] + getDistance(*(ptIJ->pt), *(ptIm1J->pt)) - _is.B;
#ifdef _DEBUG
					std::cout << "\tu:" << ptIJ->u << "; v:" << ptIJ->v << "\n";
					std::cout << "\tx:" << ptIJ->pt->x << "; y:" << ptIJ->pt->y << "; z:" << ptIJ->pt->z << "\n";
					std::cout << "Distance XJ, X_1J: " << getDistance(*(ptIJ->pt), *(ptIm1J->pt)) << "\n";
					std::cout << "Distance XJ, XJ_1: " << getDistance(*(ptIJ->pt), *(ptIJm1->pt)) << "\n";
					std::cout << "uAccum: " << uAcc[(i + 1) * size + (j + 1)] << "; vAccum: " << vAcc[(i + 1) * size + (j + 1)] << ";\n";
#endif									   
				}
			}
		}
	}
	std::cout << "\nDrapping errors: " << err_ct << "\n";

	std::ofstream pFileU, pFileV;
	pFileU.open("matrixUOut.xls");
	pFileV.open("matrixVOut.xls");
	for (int i = 0; i < size; ++i)
	{
		for (int j = 0; j < size; ++j)
		{
			pFileU << uAcc[i * size + j] << "\t";
			pFileV << vAcc[i * size + j] << "\t";
		}
		pFileU << "\n";
		pFileV << "\n";
	}
	pFileU.close();
	pFileV.close();

	drawDrappedCell(_rs, P, Q, size);

	for (UINT i = 0; i < size; ++i) delete[] Q[i], P[i];
	for (UINT i = 0; i < dim; ++i)  delete W[i], invW[i];

	delete[] Q, P, W, invW, uAcc, vAcc;
}

/// 
/// Optimized Version of functions
///

void getJakobain(double* W, const bSplinePt* ptIJ, const bSplinePt* ptIm1J, const bSplinePt* ptIJm1, const d_vertex* const* IJder)
{
	W[0] = 2 * (ptIJ->pt->x - ptIm1J->pt->x) * IJder[1][0].x + 2 * (ptIJ->pt->y - ptIm1J->pt->y) * IJder[1][0].y +
		+2 * (ptIJ->pt->z - ptIm1J->pt->z) * IJder[1][0].z;
	W[1] = 2 * (ptIJ->pt->x - ptIm1J->pt->x) * IJder[0][1].x + 2 * (ptIJ->pt->y - ptIm1J->pt->y) * IJder[0][1].y +
		+2 * (ptIJ->pt->z - ptIm1J->pt->z) * IJder[0][1].z;

	W[2] = 2 * (ptIJ->pt->x - ptIJm1->pt->x) * IJder[1][0].x + 2 * (ptIJ->pt->y - ptIJm1->pt->y) * IJder[1][0].y +
		+2 * (ptIJ->pt->z - ptIJm1->pt->z) * IJder[1][0].z;
	W[3] = 2 * (ptIJ->pt->x - ptIJm1->pt->x) * IJder[0][1].x + 2 * (ptIJ->pt->y - ptIJm1->pt->y) * IJder[0][1].y +
		+2 * (ptIJ->pt->z - ptIJm1->pt->z) * IJder[0][1].z;
}

double getAngle(const d_vertex& vx_ij, const d_vertex& vx_pq, const d_vertex& vx_sh)
{
	struct { float fi; float teta; } p0, p1, p2;
	XMFLOAT3 v1, v2;
	v1.x = vx_pq.x - vx_ij.x;
	v1.y = vx_pq.y - vx_ij.y;
	v1.z = vx_pq.z - vx_ij.z;


	v2.x = vx_sh.x - vx_ij.x;
	v2.y = vx_sh.y - vx_ij.y;
	v2.z = vx_sh.z - vx_ij.z;

	XMVECTOR angle = XMVector3AngleBetweenVectors(XMLoadFloat3(&v1), XMLoadFloat3(&v2));
	return angle.m128_f32[0] * 57.2958;
}


void drawDrappedCell(RenderSys* _rs, bSplinePt* P, d_vertex* Q, size_t size)
{
	vertex* triangle = new vertex[3];
	for (UINT i = 0; i < size - 1; ++i)
	{
		for (UINT j = 0; j < size - 1; ++j)
		{
			if ((P[i*size + j].u < 0) || (P[i * size + j + 1].u < 0) || (P[(i + 1) * size + j].u < 0))
				continue;
			float angle = getAngle(Q[i * size + j], Q[(i + 1) * size + j], Q[i * size + j + 1]);
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

			triangle[0] = Q[i * size + j];
			triangle[1] = Q[i * size + j + 1];
			triangle[2] = Q[(i + 1) * size +j];
			triangle[0].normal = calculateTriangleNormal(Q[i * size + j], Q[i * size + j + 1], Q[(i + 1) * size + j]);
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
			float angle = getAngle(Q[(i + 1) * size + j + 1], Q[(i + 1) * size + j], Q[i * size + j + 1]);
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
			if ((P[(i + 1) * size + j + 1].u < 0) || (P[i * size + j + 1].u < 0) || (P[(i + 1) * size + j].u < 0))
				continue;
			triangle[0] = Q[(i + 1) * size + j];
			triangle[1] = Q[i * size + j + 1];
			triangle[2] = Q[(i + 1) * size + j + 1];
			triangle[0].normal = calculateTriangleNormal(Q[(i + 1) * size + j], Q[i * size + j + 1], Q[(i + 1) * size + j + 1]);
			triangle[1].normal = triangle[0].normal;
			triangle[2].normal = triangle[0].normal;

			triangle[0].Color = vec4(red, green, blue, 1);
			triangle[1].Color = vec4(red, green, blue, 1);
			triangle[2].Color = vec4(red, green, blue, 1);
			_rs->drawTriangle(triangle);
		}
	}
}

bool getBSplineDrapPoint(double* W, double* invW, drappingCell& cell)
{
	const size_t maxIterCt = 10000;
	bool flag = true;
	size_t i = 0;

	size_t dim = 2;
	double* f = new double[dim];
	double* dx;

	double epsilon = 0.0001; // Варьировать сравнение ошибки в зависимости от накопленной ошибки. Т.е. если накопилась большая ошибка, то уменьшать макисмально взоможную, или прижимать в другую сторону
	double len = 100000;

	bool corrupted = false;

	int counter = 0;
	double** L = new double* [2];
	double** U = new double* [2];
	//Derivation Init
	size_t derivation_degree = 1;
	d_vertex** IJder = new d_vertex * [derivation_degree + 1];
	for (size_t i = 0; i < derivation_degree + 1; ++i) IJder[i] = new d_vertex[derivation_degree + 1];

	while (flag && (maxIterCt > i))
	{
		counter++;

		corrupted = false;
		//(*ptIJ->pt) = SurfacePoint(sfI, ptIJ->u, ptIJ->v);

		SurfaceDerivsAlg1(cell.si->sfI, cell.ptIJ->u, cell.ptIJ->v, 1, IJder);
		getJakobain(W, cell.ptIJ, cell.ptIm1J, cell.ptIJm1, IJder);
		getF(f, cell.ptIJ->pt, cell.ptIm1J->pt, cell.ptIJm1->pt, cell.si->A, cell.si->B);///@todo Optimize memory using. At each iteration memory allocation happens. It too much heavy

		//std::cout << "\n f is: " << f[0] << "\n" << f[1] << "\n";

		LUDecomposition(W, dim, &L, &U);
		dx = LUForwardBackward(L, U, f, 2);
		dx[0] *= 0.01;
		dx[1] *= 0.01;
		//if ((ptIJ->u + dx[0] < 0) || (ptIJ->v + dx[1] < 0) || (ptIJ->u + dx[0] > 1) || (ptIJ->v + dx[1] > 1)) corrupted = true;
		d_vertex oldPoint = *(cell.ptIJ->pt);

		if (cell.ptIJ->v + dx[1] < 0)
		{
			dx[1] = 0;
			cell.ptIJ->v = 0;
		}
		else
			if (cell.ptIJ->v + dx[1] > 1)
			{
				dx[1] = 0;
				cell.ptIJ->v = 1;
			}
		if (cell.ptIJ->u + dx[0] < 0)
		{
			dx[0] = 0;
			cell.ptIJ->u = 0;
		}
		else
			if (cell.ptIJ->u + dx[0] > 1)
			{
				dx[0] = 0;
				cell.ptIJ->u = 1;
			}


		cell.ptIJ->u += dx[0];
		cell.ptIJ->v += dx[1];
		(*cell.ptIJ->pt) = SurfacePoint(cell.si->sfI, cell.ptIJ->u, cell.ptIJ->v);


		for (size_t ct1 = 0; ct1 < dim; ++ct1)
		{
			delete[] L[ct1];
			delete[] U[ct1];
		}
		delete[] dx;

		double bufA = (cell.si->A) * epsilon - cell.accumULen;
		double bufB = (cell.si->B) * epsilon - cell.accumVLen;

		//if ((f[0] <= (cell.si->A) * epsilon - cell.accumULen) && (f[1] <= (cell.si->B) * epsilon - cell.accumVLen) && !corrupted)
		//if ((fabs(f[0]) < 0.01) && (fabs(f[1]) < 0.01))
		if ((fabs(getDistance(*cell.ptIJ->pt, *cell.ptIJm1->pt) - cell.si->A) < cell.si->A * epsilon)
			&& (fabs(getDistance(*cell.ptIJ->pt, *cell.ptIm1J->pt) - cell.si->B) < cell.si->B * epsilon))
		{
			flag = false;
			break;
		}
		++i;

	}
	delete[] L;
	delete[] U;
	delete[] f;

	for (int i = 0; i < derivation_degree; ++i) delete[] IJder[i];
	delete[] IJder;

	return !flag;
}

bool getBSplineDrapPoint_optmized(double* W, double* invW, drappingCell& cell)
{
	const size_t maxIterCt = 10000;
	bool flag = true;
	size_t i = 0;

	size_t dim = 2;
	double* f = new double[dim];

	double epsilon = 0.0001; // Варьировать сравнение ошибки в зависимости от накопленной ошибки. Т.е. если накопилась большая ошибка, то уменьшать макисмально взоможную, или прижимать в другую сторону
	double len = 100000;

	bool corrupted = false;

	int counter = 0;
	double* L =  new double[dim*dim];
	double* U =  new double[dim*dim];
	double* dx = new double[dim];
	double* LU_staff_y = new double[dim];
	//Derivation Init
	DerivationInit* derInit = initDerivationInitStruct(cell.si->sfI, 1);

	while (flag && (maxIterCt > i))
	{
		counter++;

		corrupted = false;
		//(*ptIJ->pt) = SurfacePoint(sfI, ptIJ->u, ptIJ->v);

		SurfaceDerivsAlg1(cell.si->sfI, cell.ptIJ->u, cell.ptIJ->v, derInit);
		getJakobain(W, cell.ptIJ, cell.ptIm1J, cell.ptIJm1, derInit->SKL);
		getF(f, cell.ptIJ->pt, cell.ptIm1J->pt, cell.ptIJm1->pt, cell.si->A, cell.si->B);///@todo Optimize memory using. At each iteration memory allocation happens. It too much heavy

		//std::cout << "\n f is: " << f[0] << "\n" << f[1] << "\n";

		LUDecomposition(W, dim, L, U);
		LUForwardBackward(L, U, f, dim, dx, LU_staff_y);
		dx[0] *= 0.01;
		dx[1] *= 0.01;
		//if ((ptIJ->u + dx[0] < 0) || (ptIJ->v + dx[1] < 0) || (ptIJ->u + dx[0] > 1) || (ptIJ->v + dx[1] > 1)) corrupted = true;
		d_vertex oldPoint = *(cell.ptIJ->pt);

		if (cell.ptIJ->v + dx[1] < 0)
		{
			dx[1] = 0;
			cell.ptIJ->v = 0;
		}
		else
			if (cell.ptIJ->v + dx[1] > 1)
			{
				dx[1] = 0;
				cell.ptIJ->v = 1;
			}
		if (cell.ptIJ->u + dx[0] < 0)
		{
			dx[0] = 0;
			cell.ptIJ->u = 0;
		}
		else
			if (cell.ptIJ->u + dx[0] > 1)
			{
				dx[0] = 0;
				cell.ptIJ->u = 1;
			}


		cell.ptIJ->u += dx[0];
		cell.ptIJ->v += dx[1];
		(*cell.ptIJ->pt) = SurfacePoint(cell.si->sfI, cell.ptIJ->u, cell.ptIJ->v);

		if ((fabs(getDistance(*cell.ptIJ->pt, *cell.ptIJm1->pt) - cell.si->A) < cell.si->A * epsilon)
			&& (fabs(getDistance(*cell.ptIJ->pt, *cell.ptIm1J->pt) - cell.si->B) < cell.si->B * epsilon))
		{
			flag = false;
			break;
		}
		++i;

	}
	delete[] L, U, f, dx, LU_staff_y;

	releaseDerivationInitStruct(cell.si->sfI, derInit);

	return !flag;
}


void generateInitialLines(bSplinePt* P, d_vertex* Q, const drappingInit& _is)
{

	double cycle_step = 1. / (_is.gird_size - 1);
	std::cout << "Generate initial lines:\n";

	P[0].u = 0;
	P[0].v = 0;

	if (_is.isU1) P[0].u = _is.u1;
	else P[0].v = _is.v1;
	if (_is.isU2) P[0].u = _is.u2;
	else P[0].v = _is.v2;

	Q[0] = SurfacePoint(_is.sfI, P[0].u, P[0].v);
	double eps = 0.0001;
	for (UINT i = 1; i < _is.gird_size; ++i)
	{
		if (_is.isU1)
		{
			P[i].u = _is.u1;
			P[i].v = i * cycle_step;
			P[i].v = bisectionV(_is.sfI, P[i], P[i - 1], 1., eps * _is.B, _is.B);
		}
		else
		{
			P[i].u = i * cycle_step;
			P[i].v = _is.v1;
			P[i].u = bisectionU(_is.sfI, P[i], P[i - 1], 1., eps * _is.A, _is.A);
		}
		if (_is.isU2)
		{
			P[i*_is.gird_size].u = _is.u2;
			P[i*_is.gird_size].v = i * cycle_step;
			P[i*_is.gird_size].v = bisectionV(_is.sfI, P[i * _is.gird_size], P[(i - 1)*_is.gird_size], 1., eps * _is.B, _is.B);
		}
		else
		{
			P[i*_is.gird_size].u = i * cycle_step;
			P[i*_is.gird_size].v = _is.v2;
			P[i*_is.gird_size].u = bisectionU(_is.sfI, P[i * _is.gird_size], P[(i - 1) * _is.gird_size], 1., eps * _is.A, _is.A);
		}
		Q[i * _is.gird_size] = SurfacePoint(_is.sfI, P[i * _is.gird_size].u, P[i * _is.gird_size].v);
		Q[i] = SurfacePoint(_is.sfI, P[i].u, P[i].v);
#ifdef _DEBUG
		std::cout << "i: " << i << ";\n\t(Q[0][i];Q[0][i - 1]) dist is: " << getDistance(Q[i], Q[i - 1]) << "\n";
		std::cout << "\t(Q[i][0];Q[i-1][0]) dist is: " << getDistance(Q[i * _is.gird_size], Q[(i - 1) * _is.gird_size]) << "\n";
#endif
	}
}
void makeDrappedGird_optimized(RenderSys* _rs, const drappingInit& _is)
{
	std::cout << "A is: " << _is.A << "; B is: " << _is.B << ";\n";

	UINT size = _is.gird_size;
	bSplinePt* P = new bSplinePt[size * size]; //points warper
	d_vertex* Q = new d_vertex[size * size];
	for (UINT i = 0; i < size; ++i)
	{
		for (UINT j = 0; j < size; ++j)
		{
			P[i*size + j].pt = &Q[i*size + j];
		}
	}
	UINT dim = 2;
	double* W = new double[dim*dim]; //Jakobian
	double* invW = new double[dim*dim]; //inverse

	UINT err_ct = 0;
	double delta_u = 0.01;
	drappingCell cell = { nullptr, nullptr, nullptr, &_is };
	//Length accumalation array by u and v coordinate;
	size_t accum_size = (size) * (size);
	double* uAcc = new double[accum_size];
	double* vAcc = new double[accum_size];
	memset(uAcc, 0, sizeof(double) * accum_size);
	memset(vAcc, 0, sizeof(double) * accum_size);

	{
		TimeBench tb;
		generateInitialLines(P, Q, _is);
	}
	{
		TimeBench tb;
		for (UINT i = 0; i < size - 1; ++i)
		{
			for (UINT j = 0; j < size - 1; ++j)
			{
				if ((P[i*size + j].u < 0) || (P[i * size + j + 1].u < 0) || (P[(i + 1) * size + j].u < 0))
					continue;

				bSplinePt* ptIJ = &P[(i + 1) * size + j + 1];
				bSplinePt* ptIm1J = &P[i * size + j + 1];
				bSplinePt* ptIJm1 = &P[(i + 1) * size + j];

				cell.ptIJ = ptIJ;
				cell.ptIJm1 = ptIJm1;
				cell.ptIm1J = ptIm1J;
				cell.accumULen = uAcc[(i + 1) * size + j];
				cell.accumVLen = vAcc[i * size + j + 1];

				//ptIJ->u = max(ptIJm1->u, ptIm1J->u);
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
				if (!getBSplineDrapPoint(W, invW, cell))
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
					uAcc[(i + 1) * size + (j + 1)] = uAcc[(i + 1) * size + j] + getDistance(*(ptIJ->pt), *(ptIJm1->pt)) - _is.A;
					vAcc[(i + 1) * size + (j + 1)] = vAcc[i * size + j + 1] + getDistance(*(ptIJ->pt), *(ptIm1J->pt)) - _is.B;
#ifdef _DEBUG
					std::cout << "\tu:" << ptIJ->u << "; v:" << ptIJ->v << "\n";
					std::cout << "\tx:" << ptIJ->pt->x << "; y:" << ptIJ->pt->y << "; z:" << ptIJ->pt->z << "\n";
					std::cout << "Distance XJ, X_1J: " << getDistance(*(ptIJ->pt), *(ptIm1J->pt)) << "\n";
					std::cout << "Distance XJ, XJ_1: " << getDistance(*(ptIJ->pt), *(ptIJm1->pt)) << "\n";
					std::cout << "uAccum: " << uAcc[(i + 1) * size + (j + 1)] << "; vAccum: " << vAcc[(i + 1) * size + (j + 1)] << ";\n";
#endif									   
				}
			}
		}
	}
	std::cout << "\nDrapping errors: " << err_ct << "\n";


	//drawDrappedCell(_rs, P, Q, size);

	delete[] Q, P, W, invW, uAcc, vAcc;
}

void makeDrappedGird_optimized_v2(RenderSys* _rs, const drappingInit& _is)
{
	std::cout << "A is: " << _is.A << "; B is: " << _is.B << ";\n";

	UINT size = _is.gird_size;
	bSplinePt* P = new bSplinePt[size * size]; //points warper
	d_vertex* Q = new d_vertex[size * size];
	for (UINT i = 0; i < size; ++i)
	{
		for (UINT j = 0; j < size; ++j)
		{
			P[i * size + j].pt = &Q[i * size + j];
		}
	}
	UINT dim = 2;
	double* W = new double[dim * dim]; //Jakobian
	double* invW = new double[dim * dim]; //inverse

	UINT err_ct = 0;
	double delta_u = 0.01;
	drappingCell cell = { nullptr, nullptr, nullptr, &_is };
	//Length accumalation array by u and v coordinate;
	size_t accum_size = (size) * (size);
	double* uAcc = new double[accum_size];
	double* vAcc = new double[accum_size];
	memset(uAcc, 0, sizeof(double) * accum_size);
	memset(vAcc, 0, sizeof(double) * accum_size);

	{
		TimeBench tb;
		generateInitialLines(P, Q, _is);
	}
	{
		TimeBench tb;
		for (UINT i = 0; i < size - 1; ++i)
		{
			for (UINT j = 0; j < size - 1; ++j)
			{
				if ((P[i * size + j].u < 0) || (P[i * size + j + 1].u < 0) || (P[(i + 1) * size + j].u < 0))
					continue;

				bSplinePt* ptIJ = &P[(i + 1) * size + j + 1];
				bSplinePt* ptIm1J = &P[i * size + j + 1];
				bSplinePt* ptIJm1 = &P[(i + 1) * size + j];

				cell.ptIJ = ptIJ;
				cell.ptIJm1 = ptIJm1;
				cell.ptIm1J = ptIm1J;
				cell.accumULen = uAcc[(i + 1) * size + j];
				cell.accumVLen = vAcc[i * size + j + 1];

				//ptIJ->u = max(ptIJm1->u, ptIm1J->u);
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
				if (!getBSplineDrapPoint_optmized(W, invW, cell))
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
					uAcc[(i + 1) * size + (j + 1)] = uAcc[(i + 1) * size + j] + getDistance(*(ptIJ->pt), *(ptIJm1->pt)) - _is.A;
					vAcc[(i + 1) * size + (j + 1)] = vAcc[i * size + j + 1] + getDistance(*(ptIJ->pt), *(ptIm1J->pt)) - _is.B;
#ifdef _DEBUG
					std::cout << "\tu:" << ptIJ->u << "; v:" << ptIJ->v << "\n";
					std::cout << "\tx:" << ptIJ->pt->x << "; y:" << ptIJ->pt->y << "; z:" << ptIJ->pt->z << "\n";
					std::cout << "Distance XJ, X_1J: " << getDistance(*(ptIJ->pt), *(ptIm1J->pt)) << "\n";
					std::cout << "Distance XJ, XJ_1: " << getDistance(*(ptIJ->pt), *(ptIJm1->pt)) << "\n";
					std::cout << "uAccum: " << uAcc[(i + 1) * size + (j + 1)] << "; vAccum: " << vAcc[(i + 1) * size + (j + 1)] << ";\n";
#endif									   
				}
			}
		}
	}
	std::cout << "\nDrapping errors: " << err_ct << "\n";


	drawDrappedCell(_rs, P, Q, size);

	delete[] Q, P, W, invW, uAcc, vAcc;
}