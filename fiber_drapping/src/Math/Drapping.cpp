// This is a personal academic project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: http://www.viva64.com
#include "../Math/Drapping.h"

void getF(double* f, const d_vertex* ptIJ, const d_vertex* ptIm1J, const d_vertex* ptIJm1, double A, double B)
{
	f[0] = -(pow((ptIJ->x - ptIm1J->x), 2) + pow((ptIJ->y - ptIm1J->y), 2) + pow((ptIJ->z - ptIm1J->z), 2) - pow(B, 2));
	f[1] = -(pow((ptIJ->x - ptIJm1->x), 2) + pow((ptIJ->y - ptIJm1->y), 2) + pow((ptIJ->z - ptIJm1->z), 2) - pow(A, 2));
}

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

bool getBSplineDrapPoint(double* W, double* invW, drappingCell& cell, drapPointInit* _is)
{
	const size_t maxIterCt = 10000;
	bool flag = true;
	size_t i = 0;

	size_t dim = 2;

	double epsilon = 0.0001; // Варьировать сравнение ошибки в зависимости от накопленной ошибки. Т.е. если накопилась большая ошибка, то уменьшать макисмально взоможную, или прижимать в другую сторону
	double len = 100000;

	while (flag && (maxIterCt > i))
	{
		SurfaceDerivsAlg1(cell.si->sfI, cell.ptIJ->u, cell.ptIJ->v, _is->derInit);
		getJakobain(W, cell.ptIJ, cell.ptIm1J, cell.ptIJm1, _is->derInit->SKL);
		getF(_is->f, cell.ptIJ->pt, cell.ptIm1J->pt, cell.ptIJm1->pt, cell.si->A, cell.si->B);


		LUDecomposition(W, dim, _is->L, _is->U);
		LUForwardBackward(_is->L, _is->U, _is->f, dim, _is->dx, _is->LU_staff_y);
		//dx[0] *= 0.01;
		//dx[1] *= 0.01;

		if (cell.ptIJ->v + _is->dx[1] < 0)
		{
			_is->dx[1] = 0;
			cell.ptIJ->v = 0;
		}
		else
			if (cell.ptIJ->v + _is->dx[1] > 1)
			{
				_is->dx[1] = 0;
				cell.ptIJ->v = 1;
			}
		if (cell.ptIJ->u + _is->dx[0] < 0)
		{
			_is->dx[0] = 0;
			cell.ptIJ->u = 0;
		}
		else
			if (cell.ptIJ->u + _is->dx[0] > 1)
			{
				_is->dx[0] = 0;
				cell.ptIJ->u = 1;
			}


		cell.ptIJ->u += _is->dx[0];
		cell.ptIJ->v += _is->dx[1];
		(*cell.ptIJ->pt) = SurfacePoint(cell.si->sfI, cell.ptIJ->u, cell.ptIJ->v);

		if ((fabs(getDistance(*cell.ptIJ->pt, *cell.ptIJm1->pt) - cell.si->A) < cell.si->A * epsilon)
			&& (fabs(getDistance(*cell.ptIJ->pt, *cell.ptIm1J->pt) - cell.si->B) < cell.si->B * epsilon))
		{
			flag = false;
			break;
		}
		++i;

	}
	return !flag;
}

void makeDrappedGird(RenderSys* _rs, const drappingInit& _is)
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
#ifdef BENCHMARK_DRAPPING
		TimeBench tb;
#endif
		generateInitialLines(P, Q, _is);
		drapPointInit* drapMem = initDrapPointInitStruct(_is.sfI, dim, 1);
		for (UINT i = 0; i < size - 1; ++i)
		{
			for (UINT j = 0; j < size - 1; ++j)
			{
				if ((P[i * size + j + 1].u < 0) || (P[(i + 1) * size + j].u < 0))
				{
					continue;
				}
				bSplinePt* ptIJ = &P[(i + 1) * size + j + 1];
				bSplinePt* ptIm1J = &P[i * size + j + 1];
				bSplinePt* ptIJm1 = &P[(i + 1) * size + j];

				cell.ptIJ = ptIJ;
				cell.ptIJm1 = ptIJm1;
				cell.ptIm1J = ptIm1J;

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
				if (!getBSplineDrapPoint(W, invW, cell, drapMem))
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
					//std::cout << err_ct << "\n";
					//std::cout << "\ti:" << i << "; j:" << j << "\n";
				}
#ifdef _DEBUG
				else
				{
					std::cout << "\tu:" << ptIJ->u << "; v:" << ptIJ->v << "\n";
					std::cout << "\tx:" << ptIJ->pt->x << "; y:" << ptIJ->pt->y << "; z:" << ptIJ->pt->z << "\n";
					std::cout << "Distance XJ, X_1J: " << getDistance(*(ptIJ->pt), *(ptIm1J->pt)) << "\n";
					std::cout << "Distance XJ, XJ_1: " << getDistance(*(ptIJ->pt), *(ptIJm1->pt)) << "\n";   
				}
#endif		
			}
		}
		releaseDrapPointInitStruct(_is.sfI, drapMem, 1);
	}
	std::cout << "\nDrapping errors: " << err_ct << "\n";

	drawDrappedCell(_rs, P, Q, size);

	delete[] Q;
	delete[] P;
	delete[] W;
	delete[] invW;
	delete[] uAcc;
	delete[] vAcc;
}

drapPointInit* initDrapPointInitStruct(const surfInfo* sfI, size_t dim, size_t der_degree)
{
	drapPointInit* res = new drapPointInit();
	res->L = new double[dim * dim];
	res->U = new double[dim * dim];
	res->dx = new double[dim];
	res->LU_staff_y = new double[dim];
	res->f = new double[dim];
	res->W = new double[dim * dim];
	res->invW = new double[dim * dim];
	res->derInit = initDerivationInitStruct(sfI, der_degree);
	return res;
}

void releaseDrapPointInitStruct(const surfInfo* sfI, drapPointInit* obj, size_t der_degree)
{
	releaseDerivationInitStruct(sfI, obj->derInit);
	delete[] obj->dx;
	delete[] obj->f;
	delete[] obj->L;
	delete[] obj->LU_staff_y;
	delete[] obj->U;
	delete[] obj->W;
	delete[] obj->invW;
	delete obj;
}


//
//Paralleld section
//

int calcOnePoint(const drappingInit& _is, drappingCell& cell, drapPointInit* drapMem, size_t i, size_t j, bSplinePt* P)
{
	static double delta_u = 0.01;

	cell.ptIJ->u = (cell.ptIJm1->u + cell.ptIm1J->u) / 2;//ptIJm1->u;
	cell.ptIJ->v = max(cell.ptIJm1->v, cell.ptIm1J->v);// ptIJm1->v - 2* delta_u;
	if ((cell.ptIJ->u == cell.ptIJm1->u) && (cell.ptIJ->v == cell.ptIJm1->v) || (cell.ptIJ->u == cell.ptIm1J->u) && (cell.ptIJ->v == cell.ptIm1J->v))
	{
		cell.ptIJ->u = (cell.ptIJm1->u + cell.ptIm1J->u) / 2;
	}
	if (cell.ptIJ->u > 1) cell.ptIJ->u -= 2 * delta_u;
	if (cell.ptIJ->v > 1) cell.ptIJ->v -= 2 * delta_u;
	if (cell.ptIJ->u < 0) cell.ptIJ->u = (cell.ptIJm1->u + cell.ptIm1J->u) / 2;//0;
	if (cell.ptIJ->v < 0) cell.ptIJ->v = min(cell.ptIJm1->v, cell.ptIm1J->v);// 0;// 1 + ptIJ->v;

	(*cell.ptIJ->pt) = SurfacePoint(_is.sfI, cell.ptIJ->u, cell.ptIJ->v);

	if (!getBSplineDrapPoint(drapMem->W, drapMem->invW, cell, drapMem))
	{
		cell.ptIJ->u = -1;
		cell.ptIJ->v = -1;
		return 1;
	}
	return 0;
}

int calcGirdLine_fixedI(const drappingInit& _is, drapPointInit* drapMem, size_t i, bSplinePt* P)
{
	static double delta_u = 0.01;
	drappingCell cell;
	cell.si = &_is;

	for (size_t j = i + 1; j < _is.gird_size - 1; ++j)
	{
		if ((P[i * _is.gird_size + j + 1].u < 0) || (P[(i + 1) * _is.gird_size + j].u < 0)) return 1;

		cell.ptIJ = &P[(i + 1) * _is.gird_size + j + 1];
		cell.ptIm1J = &P[i * _is.gird_size + j + 1];
		cell.ptIJm1 = &P[(i + 1) * _is.gird_size + j];

		cell.ptIJ->u = (cell.ptIJm1->u + cell.ptIm1J->u) / 2;//ptIJm1->u;
		cell.ptIJ->v = max(cell.ptIJm1->v, cell.ptIm1J->v);// ptIJm1->v - 2* delta_u;
		if ((cell.ptIJ->u == cell.ptIJm1->u) && (cell.ptIJ->v == cell.ptIJm1->v) || (cell.ptIJ->u == cell.ptIm1J->u) && (cell.ptIJ->v == cell.ptIm1J->v))
		{
			cell.ptIJ->u = (cell.ptIJm1->u + cell.ptIm1J->u) / 2;
		}
		if (cell.ptIJ->u > 1) cell.ptIJ->u -= 2 * delta_u;
		if (cell.ptIJ->v > 1) cell.ptIJ->v -= 2 * delta_u;
		if (cell.ptIJ->u < 0) cell.ptIJ->u = (cell.ptIJm1->u + cell.ptIm1J->u) / 2;//0;
		if (cell.ptIJ->v < 0) cell.ptIJ->v = min(cell.ptIJm1->v, cell.ptIm1J->v);// 0;// 1 + ptIJ->v;

		(*cell.ptIJ->pt) = SurfacePoint(_is.sfI, cell.ptIJ->u, cell.ptIJ->v);

		if (!getBSplineDrapPoint(drapMem->W, drapMem->invW, cell, drapMem))
		{
			cell.ptIJ->u = -1;
			cell.ptIJ->v = -1;
			return 1;
		}
	}
	return 0;
}

int calcGirdLine_fixedJ(const drappingInit& _is, drapPointInit* drapMem, size_t j, bSplinePt* P)
{
	static double delta_u = 0.01;
	drappingCell cell;
	cell.si = &_is;

	for (size_t i = j + 1; i < _is.gird_size - 1; ++i)
	{
		if ((P[i * _is.gird_size + j + 1].u < 0) || (P[(i + 1) * _is.gird_size + j].u < 0)) return 1;

		cell.ptIJ = &P[(i + 1) * _is.gird_size + j + 1];
		cell.ptIm1J = &P[i * _is.gird_size + j + 1];
		cell.ptIJm1 = &P[(i + 1) * _is.gird_size + j];

		cell.ptIJ->u = (cell.ptIJm1->u + cell.ptIm1J->u) / 2;//ptIJm1->u;
		cell.ptIJ->v = max(cell.ptIJm1->v, cell.ptIm1J->v);// ptIJm1->v - 2* delta_u;
		if ((cell.ptIJ->u == cell.ptIJm1->u) && (cell.ptIJ->v == cell.ptIJm1->v) || (cell.ptIJ->u == cell.ptIm1J->u) && (cell.ptIJ->v == cell.ptIm1J->v))
		{
			cell.ptIJ->u = (cell.ptIJm1->u + cell.ptIm1J->u) / 2;
		}
		if (cell.ptIJ->u > 1) cell.ptIJ->u -= 2 * delta_u;
		if (cell.ptIJ->v > 1) cell.ptIJ->v -= 2 * delta_u;
		if (cell.ptIJ->u < 0) cell.ptIJ->u = (cell.ptIJm1->u + cell.ptIm1J->u) / 2;//0;
		if (cell.ptIJ->v < 0) cell.ptIJ->v = min(cell.ptIJm1->v, cell.ptIm1J->v);// 0;// 1 + ptIJ->v;

		(*cell.ptIJ->pt) = SurfacePoint(_is.sfI, cell.ptIJ->u, cell.ptIJ->v);

		if (!getBSplineDrapPoint(drapMem->W, drapMem->invW, cell, drapMem))
		{
			cell.ptIJ->u = -1;
			cell.ptIJ->v = -1;
			return 1;
		}
	}
	return 0;
}

void makeDrappedGird_paralled(RenderSys* _rs, const drappingInit& _is)
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
#ifdef BENCHMARK_DRAPPING
		TimeBench tb;
#endif

		generateInitialLines(P, Q, _is);
		drapPointInit* drapMem_fixedI = initDrapPointInitStruct(_is.sfI, dim, 1);
		drapPointInit* drapMem_fixedJ = initDrapPointInitStruct(_is.sfI, dim, 1);
		for (UINT k = 0; k < size - 1; ++k)
		{
			size_t i = k;
			size_t j = k;
			if ((P[i * size + j + 1].u < 0) || (P[(i + 1) * size + j].u < 0))
				continue;

			bSplinePt* ptIJ = &P[(i + 1) * size + j + 1];
			bSplinePt* ptIm1J = &P[i * size + j + 1];
			bSplinePt* ptIJm1 = &P[(i + 1) * size + j];

			cell.ptIJ = ptIJ;
			cell.ptIJm1 = ptIJm1;
			cell.ptIm1J = ptIm1J;
			if (!calcOnePoint(_is, cell, drapMem_fixedI, k, k, P))
			{
				auto fixedI = std::async(&calcGirdLine_fixedI, _is, drapMem_fixedI, k, P);
				auto fixedJ = std::async(&calcGirdLine_fixedJ, _is, drapMem_fixedJ, k, P);
				err_ct += fixedI.get() + fixedJ.get();
			}
			else
			{
				++err_ct;
			}
		}
		releaseDrapPointInitStruct(_is.sfI, drapMem_fixedI, 1);
		releaseDrapPointInitStruct(_is.sfI, drapMem_fixedJ, 1);
	}
	std::cout << "\nDrapping errors: " << err_ct << "\n";

	drawDrappedCell(_rs, P, Q, size);

	delete[] Q;
	delete[] P;
	delete[] W;
	delete[] invW;
	delete[] uAcc;
	delete[] vAcc;
}
