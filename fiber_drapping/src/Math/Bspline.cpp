// This is a personal academic project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: http://www.viva64.com
#include "../Math/Bspline.h"

size_t splineDegree = 3;

vertex * getLinePoints(size_t _ptCount, double _step, double _r, double _z)
{
	vertex* res = new vertex[_ptCount];
	double r = _r;

	_step = DirectX::XM_PI / (_ptCount - 1);

	for (int i = 0; i < _ptCount; ++i)
	{
		res[i].pos.x = r * cos(_step * i);
		res[i].pos.y = r * sin(_step * i);
		res[i].pos.z = _z;
		res[i].Color = vec4(1.0, 1.0, 1.0, 1.0);
	}

	return res;
}
vertex * getVerticalLinePoints(size_t _ptCount, double _step, double _r, double _x)
{
	vertex* res = new vertex[_ptCount];
	double r = _r;

	_step = DirectX::XM_PI / (_ptCount - 1);

	for (int i = 0; i < _ptCount; ++i)
	{
		res[i].pos.x = _x;
		res[i].pos.y = r * sin(_step * i);
		res[i].pos.z = r * cos(_step * i);
		res[i].Color = vec4(1.0, 1.0, 1.0, 1.0);
	}

	return res;
}

vertex* getDerivatePoints(size_t _ptCount, vertex* vtx)
{
	vertex* res = new vertex[2];

	res[0] = vtx[0];
	res[0].pos.x = 0;// 2 * res[0].pos.x;

	res[0].pos.y = 0;// 2 * res[0].pos.y;
	res[0].pos.z = 0;//2 * res[0].pos.z;

	res[1] = vtx[_ptCount - 1];

	res[1].pos.x = 0;//2 * res[1].pos.x;
	res[1].pos.y = 0;//2 * res[1].pos.y;
	res[1].pos.z = 0;//2 * res[1].pos.z;
	return res;
}

void getInterpolationKnotVector(size_t pointCount, double* forwardU, double* backwardU, size_t knotVectorSize, vertex* vtx)
{

	//Calculate forwardU
	//double d = 0;
	/*for (size_t i = 1; i < pointCount; ++i)
	{
	d += sqrt((vtx[i] - vtx[i - 1]).getLength());
	}*/

	for (size_t i = 0; i < pointCount; ++i)
	{
		forwardU[i] = (double)i / (double)(pointCount - 1);
	}
	/*forwardU[0] = 0;
	for (size_t i = 1; i < pointCount - 1; ++i)
	{
	forwardU[i] = forwardU[i - 1] + sqrt((vtx[i] - vtx[i - 1]).getLength()) / d;
	}
	forwardU[pointCount - 1] = 1;*/

	for (size_t i = 1; i < pointCount - 1; ++i)
	{
		backwardU[i + 3] = forwardU[i];
	}
	for (size_t i = 0; i <= splineDegree; ++i)
	{
		backwardU[i] = 0;
		backwardU[knotVectorSize - 1 - i] = 1;
	}

}

splineInfo addInterpolationSpline(vertex* vtx, size_t pointCount)
{
	char* buf = new char[128];
	splineDegree = 3;

	HRESULT hRes = S_OK;

	double step = 0.05;
	vertex* derivatePoints = getDerivatePoints(pointCount, vtx);

	size_t knotVectorSize = (pointCount - 2) + 2 * (splineDegree + 1);
	double* forwardU = new double[pointCount];				//Индексы для уравнения Qk = Sum(0..n)[Ni,p(forwardU)*Pi]
	double* backwardU = new double[knotVectorSize];//Индексы для построения полученного спалйна
	getInterpolationKnotVector(pointCount, forwardU, backwardU, knotVectorSize, vtx);


	size_t n = pointCount - 1;
	vertex* R = new vertex[pointCount];
	for (size_t i = 3; i < n; ++i) R[i] = vtx[i - 1];

	double* dd = new double[pointCount];
	double pRes[4];
	double den = 0;

	vertex* P = new vertex[pointCount + 2];
	P[0] = vtx[0];
	P[1] = (backwardU[4] / 3) * derivatePoints[0] + P[0];
	P[pointCount + 1] = vtx[pointCount - 1];
	P[pointCount] = P[pointCount + 1] - (1 - backwardU[pointCount + 1]) / 3 * derivatePoints[1];



	BasisFuncs(4, backwardU[4], 3, backwardU, pRes);
	den = pRes[1];
	P[2] = (vtx[1] - pRes[0] * P[1]) / den;
	for (size_t i = 3; i < n; ++i)
	{
		dd[i] = pRes[2] / den;
		BasisFuncs(i + 2, backwardU[i + 2], 3, backwardU, pRes);
		den = pRes[1] - pRes[0] * dd[i];
		P[i] = (R[i] - pRes[0] * P[i - 1]) / den;
	}
	dd[n] = pRes[2] / den;
	BasisFuncs(n + 2, backwardU[n + 2], 3, backwardU, pRes);
	den = pRes[1] - pRes[0] * dd[n];
	P[n] = (vtx[n - 1] - pRes[2] * P[n + 1] - pRes[0] * P[n - 1]) / den;
	for (size_t i = n - 1; i >= 2; --i) P[i] = P[i] - dd[i + 1] * P[i + 1];


	delete[] R;
	delete[] dd;


	delete[] buf;
	//delete[] forwardU;
	delete[] derivatePoints;

	splineInfo res;
	res.controlPoints = P;
	res.cpCount = pointCount + 2;
	res.knotLength = knotVectorSize;
	res.knotVector = backwardU;
	res.forwardU = forwardU;

	return res;
}

void BasisFuncs(size_t _i, double _u, size_t _p, double * _knots, double* pRes)
{
	pRes[0] = 1.0;

	double* left = new double[_p + 1];
	double* right = new double[_p + 1];
	double saved = 0;
	double temp = 0;

	memset(left, 0, sizeof(double) * (_p + 1));
	memset(right, 0, sizeof(double) * (_p + 1));

	for (size_t j = 1; j <= _p; ++j)
	{
		left[j] = _u - _knots[_i + 1 - j];
		right[j] = _knots[_i + j] - _u;
		saved = 0;
		for (size_t r = 0; r < j; ++r)
		{
			temp = pRes[r] / (right[r + 1] + left[j - r]);
			pRes[r] = saved + right[r + 1] * temp;
			saved = left[j - r] * temp;
		}
		pRes[j] = saved;
	}

	delete[] left;
	delete[] right;
}

size_t FindSpan(size_t _n, size_t _p, double _u, double * _knot, size_t _knotSize)
{
	_n -= 1;

	if (_u == _knot[_n + 1]) return _n;
	size_t low = _p;
	size_t high = _n + 1;
	size_t mid = (low + high) / 2;

	while ((_u < _knot[mid]) || (_u >= _knot[mid + 1]))
	{
		if (_u < _knot[mid])
		{
			high = mid;
		}
		else
		{
			low = mid;
		}
		mid = (low + high) / 2;
	}

	return mid;
}

double* makeKnotVector(size_t _ptCount, size_t _q)
{
	size_t size = _q + _ptCount;
	double* knots = new double[size];

	for (size_t i = 0; i < _q; ++i)
	{
		knots[i] = 0;
	}

	for (size_t i = _q; i < _ptCount; ++i)
	{
		knots[i] = i - (int)_q + 1;
	}
	for (size_t i = _ptCount; i < size; ++i)
	{
		knots[i] = (int)_ptCount - (int)_q + 2;
	}


	return knots;
}

vertex CurvePoint(splineInfo _spI, size_t _p, double _u)
{
	vertex res;
	res.pos = vec3(0, 0, 0);
	res.Color = vec4(1, 1, 1, 1);
	double* Nq = new double[_p + 1];

	size_t span = FindSpan(_spI.cpCount, _p, _u, _spI.knotVector, _spI.knotLength);
	BasisFuncs(span, _u, _p, _spI.knotVector, Nq);
	for (size_t i = 0; i <= _p; ++i)
	{
		res = res + Nq[i] * _spI.controlPoints[span - _p + i];
	}

	delete[] Nq;
	return res;
}

vertex * makeBSpline(size_t _vxCount, size_t _p, splineInfo _spI)
{
	vertex* result = new vertex[_vxCount];

	double* _knots = _spI.knotVector;
	double dt = _spI.knotVector[_p - 1];
	double step = (_knots[_spI.cpCount] - _knots[_p - 1]) / (_vxCount - 1);
	double Nq = 0;

	for (size_t j = 0; j < _vxCount; ++j, dt += step)
	{
		result[j] = CurvePoint(_spI, _p, dt);
	}
	return result;
}

//U - knot vector; ders - result; n (n<= p) - derivate degree
double** DersBasisFuns(size_t i, double u, int p, int n, double* U)
{
	double** ders = new double*[n + 1];
	for (int q = 0; q <= n; ++q) ders[q] = new double[p + 1];
	
	double** ndu = new double*[p + 1];
	for (int q = 0; q < p + 1; ++q) ndu[q] = new double[p + 1];

	double** a = new double*[2];
	for (int q = 0; q < 2; ++q) a[q] = new double[p + 1];

	double* left = new double[p + 1];
	double* right = new double[p + 1];
	double saved = 0;
	double temp = 0;

	ndu[0][0] = 1.0; 
	for (int j = 1; j <= p; ++j)
	{
		left[j] = u - U[i + 1 - j];
		right[j] = U[i + j] - u;
		saved = 0.0;
		for (int r = 0; r < j; ++r)
		{
			ndu[j][r] = right[r + 1] + left[j - r];
			temp = ndu[r][j - 1] / ndu[j][r];

			ndu[r][j] = saved + right[r + 1] * temp;
			saved = left[j - r] * temp;
		}
		ndu[j][j] = saved;
	}

	for (int j = 0; j <= p; ++j) ders[0][j] = ndu[j][p];

	for (int r = 0; r <= p; ++r)
	{
		int s1 = 0;
		int s2 = 1;
		a[0][0] = 1.0;
		int j1 = 0;
		int j2 = 0;
		for (int k = 1; k <= n; ++k)
		{
			double d = 0;
			int rk = r - k;
			int pk = p - k;
			if (r >= k)
			{
				a[s2][0] = a[s1][0] / ndu[pk + 1][rk];
				d = a[s2][0] * ndu[rk][pk];
			}
			if (rk >= -1)	j1 = 1;
			else			j1 = -rk;

			if (r - 1 <= pk)	j2 = k - 1;
			else				j2 = p - r;

			for (int j = j1; j <= j2; ++j)
			{
				a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][rk + j];
				d += a[s2][j] * ndu[rk + j][pk];
			}

			if (r <= pk)
			{
				a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r];
				d += a[s2][k] * ndu[r][pk];
			}
			ders[k][r] = d;
			//Switch rows//
			int j;

			j = s1; 
			s1 = s2; 
			s2 = j;
		}
	}

	//Multiply through by the correct factors//
	double r = p;
	for (int k = 1; k <= n; ++k)
	{
		for (int j = 0; j <= p; ++j) ders[k][j] *= r;
		r *= (p - k);
	}

	for (int q = 0; q < p + 1; ++q) delete[] ndu[q];
	delete[] ndu;
	for (int q = 0; q < 2; ++q) delete[] a[q];
	delete[] a;
	delete[] left;
	delete[] right;

	return ders;
}



vertex* CurveDerivateAlg1(splineInfo spi, size_t p, double u, size_t d)
{
	int n = spi.cpCount;
	double* U = spi.knotVector; 
	vertex* P = spi.controlPoints;

	vertex* CK = new vertex[d + 1];
	memset(CK, 0, sizeof(vertex) * (d + 1));
	double du;
	d <= p ? du = d : du = p;

	vertex nullVx;
	nullVx.pos = vec3(0, 0, 0);
	nullVx.Color = vec4(0, 0, 0, 0);
	nullVx.normal = vec3(0, 0, 0);

	for (int k = p + 1; k < -(int)d; ++k) CK[k] = nullVx;

	int span = FindSpan(n, p, u, U, spi.knotLength);
	double** nders = DersBasisFuns(span, u, p, du, U);
	for (int k = 0; k < du; ++k)
	{
		CK[k] = nullVx;
		for (int j = 0; j <= p; ++j) CK[k] = CK[k] + nders[k][j] * P[span - p + j];
	}

	for (int i = 0; i < du + 1; ++i) delete[] nders[i];
	delete[] nders;

	return CK;
}

void SurfMeshParams(size_t n, size_t m, vertex** Q, double** uk, double** vl)
{
	size_t num = m;
	(*uk) = new double[n];
	(*uk)[0] = 0;
	(*uk)[n - 1] = 1;
	double d = 0;

	for (size_t k = 1; k < n - 1; ++k) (*uk)[k] = 0;
	double total = 0;
	double* cds = new double[n];
	vec3 diff;
	float fDiff;
	for (size_t l = 0; l < m; ++l)
	{
		total = 0;
		for (size_t k = 1; k < n; ++k)
		{
			diff = Q[k][l].pos - Q[k - 1][l].pos;
			cds[k] = sqrtf(DirectX::XMVectorGetX(DirectX::XMVector3Dot(DirectX::XMLoadFloat3(&diff), DirectX::XMLoadFloat3(&diff))));
			total = total + cds[k];
		}
		if (total == 0) num = num - 1;
		else
		{
			d = 0;
			for (size_t k = 1; k < n - 1; ++k)
			{
				d = d + cds[k];
				(*uk)[k] = (*uk)[k] + d / total;
			}
		}
	}
	if (num == 0)
	{
		delete[] cds;
		throw "Error!";
	}
	for (int k = 1; k < (n - 1); k++)
	{
		(*uk)[k] = (*uk)[k] / num;
	}
	delete[] cds;
	//vl
	num = n;
	(*vl) = new double[m];
	(*vl)[0] = 0;
	(*vl)[m - 1] = 1;

	for (size_t k = 1; k < m - 1; ++k) (*vl)[k] = 0;
	total = 0;
	cds = new double[m];

	for (size_t l = 0; l < n; ++l)
	{
		total = 0;
		for (size_t k = 1; k < m; ++k)
		{
			diff = Q[k][l].pos - Q[k - 1][l].pos;
			cds[k] = sqrtf(DirectX::XMVectorGetX(DirectX::XMVector3Dot(DirectX::XMLoadFloat3(&diff), DirectX::XMLoadFloat3(&diff))));
			total = total + cds[k];
		}
		if (total == 0) num = num - 1;
		else
		{
			d = 0;
			for (size_t k = 1; k < m - 1; ++k)
			{
				d = d + cds[k];
				(*vl)[k] = (*vl)[k] + d / total;
			}
		}
	}
	if (num == 0)
	{
		delete[] cds;
		throw "Error!";
	}
	for (size_t k = 1; k < m - 1; k++)	(*vl)[k] = (*vl)[k] / num;

#ifdef LOG_ON
	std::cout << "SurfMeshParams, uk is:\n";
	for (size_t i = 0; i < n; ++i) std::cout << (*uk)[i] << " ";

	std::cout << "\n\tvl is:\n";
	for (size_t i = 0; i < m; ++i) std::cout << (*vl)[i] << " ";
#endif
	delete[] cds;
}

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

float* getVxCoordByPosNum(vertex& vx, size_t num)
{
	switch (num)
	{
	case 0:
		return &vx.pos.x;
	case 1:
		return &vx.pos.y;
	case 2:
		return &vx.pos.z;
	default:
		throw "Uncorrect index!";
	}
}

vertex* interpolateCurve(vertex* vtx, size_t vtx_ct, size_t p, double* uk, double* U, size_t r)
{
	//size_t m = vtx_ct + p;
	if (!uk)
	{
		throw "not implemented!";
	}

	double** A = new double* [vtx_ct];
	for (size_t i = 0; i < vtx_ct; ++i) A[i] = (double*)calloc(vtx_ct, sizeof(double));

	size_t span = 0;
	for (size_t i = 0; i < vtx_ct; ++i)
	{
		span = FindSpan(vtx_ct, p, uk[i], U, vtx_ct + p); //knot size from computeKnotVector function
		BasisFuncs(span, uk[i], p, U, (A[i] + span - p));
	}
	double** Lmx = new double*[vtx_ct];
	double** Umx = new double*[vtx_ct];

	LUDecomposition(A, vtx_ct, &Lmx, &Umx);
	double* rhs = new double[vtx_ct];
	vertex* P = new vertex[vtx_ct];
	for (size_t i = 0; i < r; ++i)
	{
		for (size_t j = 0; j < vtx_ct; ++j)	rhs[j] = *getVxCoordByPosNum(vtx[j], i);
		double* sol = LUForwardBackward(Lmx, Umx, rhs, vtx_ct);
		for (size_t j = 0; j < vtx_ct; ++j) *getVxCoordByPosNum(P[j], i) = sol[j];
	}
	return P;
}

double* computeKnotVector(double* _u, size_t p, size_t n)
{
	size_t m = n + p + 1;
	double* u = (double*)calloc(m, sizeof(double));

	for (size_t i = 0; i <= p; ++i) u[i] = 0;
	for (size_t j = 1; j < n - p; ++j)
	{

		for (size_t i = j; i <= j + p - 1; ++i)
		{
			u[p + j] += _u[i];
		}
		u[p + j] /= p;
	}
	for (size_t i = m - p - 1; i < m; ++i) u[i] = 1;

#ifdef LOG_ON
	std::cout << "\nComputed vector:\n";
	for (size_t i = 0; i < m; ++i) std::cout << u[i] << " ";
	std::cout << "\n";
#endif
	return u;
}

vertex* testSpline()
{
	size_t ct = 5;
	vertex* vtx = new vertex[ct];
	vtx[0].pos = vec3(0, 0, 0);
	vtx[1].pos = vec3(3, 4, 0);
	vtx[2].pos = vec3(-1, 4, 0);
	vtx[3].pos = vec3(-4, 0, 0);
	vtx[4].pos = vec3(-4, -3, 0);
	
	double* uk = new double[ct];
	double* U = (double*) calloc(ct + 3 + 1, sizeof(double));
	uk[0] = 0;
	uk[1] = 5. / 17;
	uk[2] = 9. / 17;
	uk[3] = 14. / 17;
	uk[4] = 1.;

	U[4] = 28. / 51;
	for (size_t i = 5; i < ct + 3 + 1; ++i) U[i] = 1;

	vertex* P = interpolateCurve(vtx, ct, 3, uk, U);

	splineInfo spi;
	spi.controlPoints = P;
	spi.cpCount = ct;
	spi.forwardU = uk;
	spi.knotLength = ct + 3 + 1;
	spi.knotVector = U;

	vertex* res = new vertex[1000];
	for (size_t i = 0; i < 1000; ++i)
	{
		res[i] = CurvePoint(spi, 3, i * 1. / 999);
	}

	return res;
}

vertex** transposeMatrix(vertex** mx, size_t n, size_t m)
{
	vertex** res = new vertex * [m];
	for (size_t i = 0; i < m; ++i)
	{
		res[i] = new vertex[n];
		for (size_t j = 0; j < n; ++j)
		{
			res[i][j] = mx[j][i];
		}
	}
	return res;
}

void saveAsTransponsed(vertex** mx, size_t n, size_t m, size_t idx, vertex* vec)
{
	for (size_t i = 0; i < m; ++i)
	{
		mx[i][idx] = vec[i];
	}
}

surfInfo GenInterpBSplineSurface(size_t n, size_t m, vertex** Q, size_t p, size_t q)
{
	double* _uk, * _vl;
	SurfMeshParams(n, m, Q, &_uk, &_vl); //size of knots: n and m
	double* Uk = computeKnotVector(_uk, p, n);
	double* Vl = computeKnotVector(_vl, q, m);

	//Interpolate by t coordinate
	vertex** R = new vertex * [m];
	for (size_t i = 0; i < m; ++i) R[i] = new vertex[n];

	for (size_t l = 0; l < m; ++l)
	{
		saveAsTransponsed(R, n, m, l, interpolateCurve(Q[l], n, p, _uk, Uk));
	}

	//Interpolate by tau coordinate
	vertex** P = new vertex * [n];
	for (size_t i = 0; i < n; ++i) P[i] = new vertex[m];

	for (size_t l = 0; l < n; ++l)
	{
		saveAsTransponsed(P, m, n, l, interpolateCurve(R[l], m, q, _vl, Vl));
	}

	delete[] _uk;
	delete[] _vl;
	delete[] R;
	surfInfo res(P, n, m, p, q, Uk, Vl);
	return res;
}

vertex SurfacePoint(surfInfo* sfI, double u, double v)
{
	size_t uspan = FindSpan(sfI->n, sfI->p, u, sfI->Uk, sfI->n + sfI->p);
	double* Nu, * Nv;
	Nu = new double[sfI->p + 1];
	BasisFuncs(uspan, u, sfI->p, sfI->Uk, Nu);
	Nv = new double[sfI->q + 1];
	size_t vspan = FindSpan(sfI->m, sfI->q, v, sfI->Vl, sfI->m + sfI->q);
	BasisFuncs(vspan, v, sfI->q, sfI->Vl, Nv);
	size_t uind = uspan - sfI->p;
	
	vertex S = vertex();
	for (size_t l = 0; l <= sfI->q; ++l)
	{
		vertex temp = vertex();
		size_t vind = vspan - sfI->q + l;
		for (size_t k = 0; k <= sfI->p; ++k)
		{
			temp = temp + Nu[k] * sfI->controlPoints[uind + k][vind];
		}
		S += Nv[l] * temp;
	}
	delete[] Nu;
	delete[] Nv;
	return S;
}

vertex** SurfaceDerivsAlg1(surfInfo* sfI, double u, double v, size_t d)
{
	vertex** SKL = new vertex * [d + 1];
	for (size_t i = 0; i < d + 1; ++i) SKL[i] = new vertex[d + 1];

	size_t du = min(sfI->p, d);
	size_t dv = min(sfI->q, d);
	for (size_t k = sfI->p + 1; k <= d; ++k)
		for (size_t l = 0; l <= d - k; ++l) SKL[k][l] = vertex();

	for (size_t l = sfI->q + 1; l <= d; ++l)
		for (size_t k = 0; k <= d - l; ++k) SKL[k][l] = vertex();

	size_t uspan = FindSpan(sfI->n, sfI->p, u, sfI->Uk, sfI->n + sfI->p + 1);
	size_t vspan = FindSpan(sfI->m, sfI->q, v, sfI->Vl, sfI->m + sfI->q + 1);
	
	double** Nu = DersBasisFuns(uspan, u, sfI->p, sfI->n, sfI->Uk);
	double** Nv = DersBasisFuns(vspan, v, sfI->q, sfI->m, sfI->Vl);
	
	vertex* temp = new vertex[sfI->q + 1];
	
	for (size_t k = 0; k <= du; ++k)
	{
		for (size_t s = 0; s <= sfI->q; ++s)
		{
			temp[s] = vertex();
			for (size_t r = 0; r <= sfI->p; ++r)
			{
				temp[s] = temp[s] + Nu[k][r] * sfI->controlPoints[uspan - sfI->p + r][vspan - sfI->q + s];
			}
		}
		size_t dd = min(d - k, dv);
		for (size_t l = 0; l <= dd; ++l)
		{
			SKL[k][l] = vertex();
			for (size_t s = 0; s <= sfI->q; ++s)
			{
				SKL[k][l] = SKL[k][l] + Nv[l][s] * temp[s];
			}
		}
	}
	for (int i = 0; i <= sfI->n; ++i) delete[] Nu[i];
	for (int i = 0; i <= sfI->m; ++i) delete[] Nv[i];
	delete[] Nu;
	delete[] Nv;
	delete[] temp;

	return SKL;
}

double getSplineLen(double _left, double _right, vertex *(*ffunc)(splineInfo, size_t, double, size_t), splineInfo _spi, size_t _p, size_t n)
{
	vec3 v1 = integrate(0, 1, *CurveDerivateAlg1, _spi);
	return sqrt(pow(v1.x, 2) + pow(v1.y, 2) + pow(v1.z, 2));
}

void getKnotForApproximationSurf(double* U, double* ub, int m, int n, int p)
{
	double d = (m + 1) / (n - p + 1);
	double a = 0;
	int i = 0;
	for (int x = 0;  x <= p; ++x)
	{
		U[x] = 0;
	}
	for (int j = 1; j < n - p; ++j)
	{
		i = (int)(j * d);
		a = j * d - i;
		U[p + j] = (1 - a) * ub[i - 1] + a * ub[i];
	}
	for (int j = n-p; j < n + p + 2; ++j)
	{
		U[p + j] = 1;
	}
}

void makeNAproximateionMatrix_new(double*** N, double* u, double* _knot, int knot_size, int n, int m, int p)
{
	//m -= 2;
	(*N) = new double* [m - 2];
	for (int i = 0; i < m - 2; ++i)
	{
		(*N)[i] = new double[n - 2];
		memset((*N)[i], 0, sizeof(double) * (n - 2));
	}
	double* pRes = new double[p + 1];
	for (int i = 1; i < m - 1; ++i)
	{
		int zspan = FindSpan(n, p, u[i], _knot, knot_size);
		BasisFuncs(zspan, u[i], p, _knot, pRes);
		for (int pi = 0; pi < p; ++pi)
		{
			(*N)[i - 1][zspan + pi] = pRes[pi];
		}
	}
	delete[] pRes;
}

void makeNAproximateionMatrix(double*** N, double* u, int n, int m, int p, int u_size)
{
	//m -= 2;
	(*N) = new double* [m - 2];
	for (int i = 0; i < m - 2; ++i)
	{
		(*N)[i] = new double[n - 2];
		memset((*N)[i], 0, sizeof(double) * (n - 2));
	}
	double* pRes = new double[p + 1];
	for (int i = 1; i < m - 1; ++i)
	{
		int zspan = FindSpan(m-2, p, u[i], u, u_size);
		BasisFuncs(zspan, u[i], p, u, pRes);
		for (int pi = 0; pi < p; ++pi)
		{
			(*N)[i - 1][zspan + pi] = pRes[pi];
		}
	}
	delete[] pRes;
}

surfInfo BSplineSurface(int r, int s, vertex** const Q, int p, int q, double** U, double** V, vertex*** P, int u_pt_cnt, int v_pt_cnt)
{
	double* ub, * vb;
	SurfMeshParams(r, s, Q, &ub, &vb); 
	(*U) = computeKnotVector(ub, p, r);
	(*V) = computeKnotVector(vb, q, s);
	

	(*P) = new vertex*[u_pt_cnt];
	for (int i = 0; i < u_pt_cnt; ++i)
	{
		(*P)[i] = new vertex[v_pt_cnt];
	}

	surfInfo surfInf;
	surfInf.controlPoints = Q;
	surfInf.n = r;
	surfInf.m = s;
	surfInf.p = 3;
	surfInf.q = 3;
	surfInf.Uk = *U;
	surfInf.Vl = *V;
	for (int i = 0; i < u_pt_cnt; ++i)
	{
		std::cout << i << "\n";
		for (int j = 0; j < v_pt_cnt; ++j)
		{
			(*P)[i][j] = SurfacePoint(&surfInf, (float)i / (u_pt_cnt - 1), (float)j / (v_pt_cnt - 1));
		}
	}
	delete[] ub, vb;
	return surfInf;
}


void transponateMatrix(double** _N, int _n, int _m, double*** _pRes)
{
	(*_pRes) = new double* [_m];
	for (int i = 0; i < _m; ++i)
	{
		(*_pRes)[i] = new double[_n];
	}
	for (int i = 0; i < _n; ++i)
	{
		for (int j = 0; j < _m; ++j)
		{
			(*_pRes)[j][i] = _N[i][j];
		}
	}
}

//Return multiple of two matrix N and transponated N. N is (m-2) * (n-2) matrix.
void makeNTNuAproxMx(double** _N, double** _Nt, double*** _pRes, int _n, int _m)
{
	(*_pRes) = new double* [_n];
	for (int i = 0; i < _n; ++i)
	{
		(*_pRes)[i] = new double[_n];
		memset((*_pRes)[i], 0, sizeof(double) * _n);
	}
	for (int i = 0; i < _n; ++i)
	{
		for (int j = 0; j < _n; ++j)
		{
			for (int k = 0; k < _m; ++k)
			{
				(*_pRes)[i][j] += _N[i][k] + _Nt[k][j];
			}
		}
	}
}


void GlobalAproximationBSpline(int r, int s, vertex** const Q, int p, int q, int n, int m, double** U, double** V, vertex*** P)
{
	int U_knot_size = n + p + 1;
	int V_knot_size = m + q + 1;
	(*U) = new double[U_knot_size];
	(*V) = new double[V_knot_size];

	double* ub, *vb;
	SurfMeshParams(r, s, Q, &ub, &vb); // Get interpolation params u|, v|
	getKnotForApproximationSurf(*U, ub, n, r, p); //Get knot vector U
	getKnotForApproximationSurf(*V, vb, m, s, q); //Get knot vector V
	
	
	//U direction fits
	double** Nu, ** Nut, ** NNut;
	int Nu_u_size = r - 2;
	int Nu_v_size = n - 2;
	//makeNAproximateionMatrix(&Nu, *U, n, r, p, U_knot_size);
	makeNAproximateionMatrix_new(&Nu, ub, (*U), U_knot_size, n, r, p);
	transponateMatrix(Nu, Nu_u_size, Nu_v_size, &Nut);
	makeNTNuAproxMx(Nut, Nu, &NNut, Nu_v_size, Nu_u_size);
	std::ofstream pFile;
	pFile.open("mx_Nu_log.txt");
	for (int i = 0; i < Nu_u_size; ++i)
	{
		for (int j = 0; j < Nu_v_size; ++j)
		{
			pFile << Nu[i][j] << "\t";
		}
		pFile << "\n";
	}
	pFile.close();
	pFile.open("mx_NNut_log.txt");
	for (int i = 0; i < Nu_v_size; ++i)
	{
		for (int j = 0; j < Nu_v_size; ++j)
		{
			pFile << NNut[i][j] << "\t";
		}
		pFile << "\n";
	}
	pFile.close();

	double** Lar, ** Uar;
	vertex **Temp = new vertex*[n];
	for (int i = 0; i < r; ++i)
	{
		Temp[i] = new vertex[s];
	}

	Lar = new double* [Nu_v_size];
	Uar = new double* [Nu_v_size];
	LUDecomposition(NNut, Nu_v_size, &Lar, &Uar);
	
	//Get R Column
	double* Rx = new double[n - 2];
	double* Ry = new double[n - 2];
	double* Rz = new double[n - 2];
	double* basisN = new double[p + 1];

	for (int j = 0; j < s; ++j)
	{
		Temp[0][j] = Q[0][j];
		Temp[n][j] = Q[r][j];
		for (int i = 1; i < n - 1; ++i)
		{
			BasisFuncs(0, ub[i], p, (*U), basisN);
			Rx[i] = Q[i][j].pos.x - basisN[0] * Q[0][j].pos.x;
			Ry[i] = Q[i][j].pos.y - basisN[0] * Q[0][j].pos.y;
			Rz[i] = Q[i][j].pos.z - basisN[0] * Q[0][j].pos.z;
			BasisFuncs(n, ub[i], p, (*U), basisN);
			Rx[i] -= basisN[0] * Q[r][j].pos.x;
			Ry[i] -= basisN[0] * Q[r][j].pos.y;
			Rz[i] -= basisN[0] * Q[r][j].pos.z;
		}
		double* pResX = LUForwardBackward(Lar, Uar, Rx, Nu_u_size);
		double* pResY = LUForwardBackward(Lar, Uar, Ry, Nu_u_size);
		double* pResZ = LUForwardBackward(Lar, Uar, Rz, Nu_u_size);
		for (int i = 1; i < n - 1; ++i)
		{
			Temp[i][j].pos.x = pResX[i - 1];
			Temp[i][j].pos.y = pResY[i - 1];
			Temp[i][j].pos.z = pResZ[i - 1];
		}
		delete[] pResX,  pResY, pResZ;
	}
	delete[] Rx, Ry, Rz, basisN;
	for (int i = 0; i < Nu_u_size; ++i)
	{
		delete[] Nu[i], Lar[i], Uar[i];
	}
	delete[] Nu, Lar, Uar;
	for (int i = 0; i < Nu_v_size; ++i)
	{
		delete[] NNut[i];
		delete[] Nut[i];
	}
	delete[] NNut;
	delete[] Nut;
	
	/////////////////V directions fits
	double** Nv, ** Nvt, ** NNvt;
	int Nv_u_size = s - 2;
	int Nv_v_size = m - 2;
	makeNAproximateionMatrix(&Nv, *V, m, s, q, V_knot_size);
	transponateMatrix(Nv, Nv_u_size, Nv_v_size, &Nvt);
	makeNTNuAproxMx(Nv, Nvt, &NNvt, Nv_u_size, Nv_v_size);

	vertex** _P = (*P);
	_P = new vertex * [n];
	for (int i = 0; i < n; ++i)
	{
		_P[i] = new vertex[m];
	}

	Rx = new double[m - 2];
	Ry = new double[m - 2];
	Rz = new double[m - 2];
	basisN = new double[q + 1];
	LUDecomposition(NNvt, Nv_v_size, &Lar, &Uar);
	for (int i = 0; i < n; ++i)
	{
		_P[i][0] = Temp[i][0];
		_P[i][m] = Temp[i][s];
		for (int j = 1; j < m - 1; ++j)
		{
			BasisFuncs(0, vb[j], q, (*V), basisN);
			Rx[i] = Q[j][i].pos.x - basisN[0] * Q[j][0].pos.x;
			Ry[i] = Q[j][i].pos.y - basisN[0] * Q[j][0].pos.y;
			Rz[i] = Q[j][i].pos.z - basisN[0] * Q[j][0].pos.z;
			BasisFuncs(m, vb[j], q, (*V), basisN);
			Rx[i] -= basisN[0] * Q[j][s].pos.x;
			Ry[i] -= basisN[0] * Q[j][s].pos.y;
			Rz[i] -= basisN[0] * Q[j][s].pos.z;
		}
		double* pResX = LUForwardBackward(Lar, Uar, Rx, Nv_u_size);
		double* pResY = LUForwardBackward(Lar, Uar, Ry, Nv_u_size);
		double* pResZ = LUForwardBackward(Lar, Uar, Rz, Nv_u_size);
		for (int j = 1; j < m - 1; ++j)
		{
			_P[i][j].pos.x = pResX[j - 1];
			_P[i][j].pos.y = pResY[j - 1];
			_P[i][j].pos.z = pResZ[j - 1];
		}
		delete[] pResX, pResY, pResZ;
	}
	delete[] Rx, Ry, Rz, basisN;
	
	for (int i = 0; i < Nv_u_size; ++i)
	{
		delete[] Nv[i], Lar[i], Uar[i];
	}
	delete[] Nv, Lar, Uar;
	for (int i = 0; i < Nv_v_size; ++i)
	{
		delete[] Nvt[i], NNvt[i];
	}
	delete[] Nvt, NNvt;

	for (int i = 0; i < n; ++i) delete[] Temp[i];
	delete[] Temp;
	delete[] ub, vb;
}