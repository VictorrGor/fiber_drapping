#include "Bspline.h"

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
		throw "Error!";
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
		throw "Error!";
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
		if (!L[i][i]) throw "LUForwardBackward: L[i][i] is null!";
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
	return S;
}

double getSplineLen(double _left, double _right, vertex *(*ffunc)(splineInfo, size_t, double, size_t), splineInfo _spi, size_t _p, size_t n)
{
	vec3 v1 = integrate(0, 1, *CurveDerivateAlg1, _spi);
	return sqrt(pow(v1.x, 2) + pow(v1.y, 2) + pow(v1.z, 2));
}
