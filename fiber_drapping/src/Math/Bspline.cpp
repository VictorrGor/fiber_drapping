// This is a personal academic project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: http://www.viva64.com
#include "../Math/Bspline.h"


d_vertex* getDerivatePoints(size_t _ptCount, d_vertex* vtx)
{
	d_vertex* res = new d_vertex[2];

	res[0] = vtx[0];
	res[0].x = 0;// 2 * res[0].pos.x;

	res[0].y = 0;// 2 * res[0].pos.y;
	res[0].z = 0;//2 * res[0].pos.z;

	res[1] = vtx[_ptCount - 1];

	res[1].x = 0;//2 * res[1].pos.x;
	res[1].y = 0;//2 * res[1].pos.y;
	res[1].z = 0;//2 * res[1].pos.z;
	return res;
}

void getInterpolationKnotVector(size_t pointCount, double* forwardU, double* backwardU, size_t knotVectorSize, d_vertex* vtx, size_t splineDegree)
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

splineInfo addInterpolationSpline(d_vertex* vtx, size_t pointCount, size_t splineDegree)
{
	char* buf = new char[128];

	HRESULT hRes = S_OK;

	double step = 0.05;
	d_vertex* derivatePoints = getDerivatePoints(pointCount, vtx);

	size_t knotVectorSize = (pointCount - 2) + 2 * (splineDegree + 1);
	double* forwardU = new double[pointCount];				
	double* backwardU = new double[knotVectorSize];
	getInterpolationKnotVector(pointCount, forwardU, backwardU, knotVectorSize, vtx, splineDegree);


	size_t n = pointCount - 1;
	d_vertex* R = new d_vertex[pointCount];
	for (size_t i = 3; i < n; ++i) R[i] = vtx[i - 1];

	double* dd = new double[pointCount];
	double pRes[4];
	double den = 0;

	d_vertex* P = new d_vertex[pointCount + 2];
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

d_vertex CurvePoint(splineInfo _spI, size_t _p, double _u)
{
	d_vertex res = { 0., 0., 0. };
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

d_vertex* makeBSpline(size_t _vxCount, size_t _p, splineInfo _spI)
{
	d_vertex* result = new d_vertex[_vxCount];

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

			if (r-1 <= pk)	j2 = k - 1;
			else			j2 = p - r;

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



d_vertex* CurveDerivateAlg1(splineInfo spi, size_t p, double u, size_t d)
{
	int n = spi.cpCount;
	double* U = spi.knotVector;
	d_vertex* P = spi.controlPoints;

	d_vertex* CK = new d_vertex[d + 1];
	memset(CK, 0, sizeof(d_vertex) * (d + 1));
	double du;
	d <= p ? du = d : du = p;

	d_vertex nullVx = {0, 0, 0};

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

void SurfMeshParams(size_t n, size_t m, d_vertex** Q, double** uk, double** vl)
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
	for (size_t l = 0; l < m; ++l)
	{
		total = 0;
		for (size_t k = 1; k < n; ++k)
		{
			diff = subtructAsVec3(Q[k][l], Q[k - 1][l]);
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
			diff = subtructAsVec3(Q[k][l], Q[k - 1][l]);
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

double* getVxCoordByPosNum(d_vertex& vx, size_t num)
{
	switch (num)
	{
	case 0:
		return &vx.x;
	case 1:
		return &vx.y;
	case 2:
		return &vx.z;
	default:
		throw "Uncorrect index!";
	}
}

d_vertex* interpolateCurve(d_vertex* vtx, size_t vtx_ct, size_t p, double* uk, double* U, size_t r)
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
	d_vertex* P = new d_vertex[vtx_ct];
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


void saveAsTransponsed(d_vertex** mx, size_t n, size_t m, size_t idx, d_vertex* vec)
{
	for (size_t i = 0; i < m; ++i)
	{
		mx[i][idx] = vec[i];
	}
}

surfInfo GenInterpBSplineSurface(size_t n, size_t m, d_vertex** Q, size_t p, size_t q)
{
	double* _uk, * _vl;
	SurfMeshParams(n, m, Q, &_uk, &_vl); //size of knots: n and m
	double* Uk = computeKnotVector(_uk, p, n);
	double* Vl = computeKnotVector(_vl, q, m);

	//Interpolate by t coordinate
	d_vertex** R = new d_vertex * [m];
	for (size_t i = 0; i < m; ++i) R[i] = new d_vertex[n];

	for (size_t l = 0; l < m; ++l)
	{
		saveAsTransponsed(R, n, m, l, interpolateCurve(Q[l], n, p, _uk, Uk));
	}

	//Interpolate by tau coordinate
	d_vertex** P = new d_vertex * [n];
	for (size_t i = 0; i < n; ++i) P[i] = new d_vertex[m];

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

d_vertex SurfacePoint(const surfInfo* sfI, double u, double v)
{
	size_t uspan = FindSpan(sfI->n, sfI->p, u, sfI->Uk, sfI->n + sfI->p);
	double* Nu, * Nv;
	Nu = new double[sfI->p + 1];
	BasisFuncs(uspan, u, sfI->p, sfI->Uk, Nu);
	Nv = new double[sfI->q + 1];
	size_t vspan = FindSpan(sfI->m, sfI->q, v, sfI->Vl, sfI->m + sfI->q);
	BasisFuncs(vspan, v, sfI->q, sfI->Vl, Nv);
	size_t uind = uspan - sfI->p;
	
	d_vertex S = d_vertex();
	for (size_t l = 0; l <= sfI->q; ++l)
	{
		d_vertex temp = d_vertex();
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

d_vertex** SurfaceDerivsAlg1(const surfInfo* sfI, double u, double v, size_t d)
{
	d_vertex** SKL = new d_vertex * [d + 1];
	for (size_t i = 0; i < d + 1; ++i) SKL[i] = new d_vertex[d + 1];

	size_t du = min(sfI->p, d);
	size_t dv = min(sfI->q, d);
	for (size_t k = sfI->p + 1; k <= d; ++k)
		for (size_t l = 0; l <= d - k; ++l) SKL[k][l] = d_vertex();

	for (size_t l = sfI->q + 1; l <= d; ++l)
		for (size_t k = 0; k <= d - l; ++k) SKL[k][l] = d_vertex();

	size_t uspan = FindSpan(sfI->n, sfI->p, u, sfI->Uk, sfI->n + sfI->p + 1);
	size_t vspan = FindSpan(sfI->m, sfI->q, v, sfI->Vl, sfI->m + sfI->q + 1);
	
	double** Nu = DersBasisFuns(uspan, u, sfI->p, d, sfI->Uk);
	double** Nv = DersBasisFuns(vspan, v, sfI->q, d, sfI->Vl);
	
	/*std::cout << "Nu:\n";
	for (int i = 0; i <= sfI->p; ++i) std::cout << "\ti:" << i << ": " << Nu[1][i] << ";\n";
	std::cout << "Nv:\n";
	for (int i = 0; i <= sfI->q; ++i) std::cout << "\ti:" << i << ": " << Nv[1][i] << ";\n";*/

	d_vertex* temp = new d_vertex[sfI->q + 1];
	
	for (size_t k = 0; k <= du; ++k)
	{
		for (size_t s = 0; s <= sfI->q; ++s)
		{
			temp[s] = d_vertex();
			for (size_t r = 0; r <= sfI->p; ++r)
			{
				temp[s] = temp[s] + Nu[k][r] * sfI->controlPoints[uspan - sfI->p + r][vspan - sfI->q + s];
			}
		}
		size_t dd = min(d - k, dv);
		for (size_t l = 0; l <= dd; ++l)
		{
			SKL[k][l] = d_vertex();
			for (size_t s = 0; s <= sfI->q; ++s)
			{
				SKL[k][l] = SKL[k][l] + Nv[l][s] * temp[s];
			}
		}
	}
	for (int i = 0; i <= d; ++i) delete[] Nu[i];
	for (int i = 0; i <= d; ++i) delete[] Nv[i];
	delete[] Nu;
	delete[] Nv;
	delete[] temp;

	return SKL;
}

///SKL[d+1][d+1]; d - derivation degree


DersBasisFunsInit* initDersBasisFunsStruct(size_t p, size_t n)
{
	DersBasisFunsInit* res = new DersBasisFunsInit();
	res->ndu = new double* [p + 1];
	res->a = new double* [2];
	res->ders = new double* [n + 1];
	for (size_t i = 0; i < p + 1; ++i) res->ndu[i] = new double[p + 1];
	for (size_t i = 0; i < 2; ++i) res->a[i] = new double[p + 1];
	for (size_t i = 0; i < n + 1; ++i) res->ders[i] = new double[p + 1];

	res->left = new double[p + 1];
	res->right = new double[p + 1];

	return res;
}

void releaseDersBasisFunsStruct(size_t p, size_t n, DersBasisFunsInit* _obj)
{
	for (size_t i = 0; i < p + 1; ++i) delete[] _obj->ndu[i];
	for (size_t i = 0; i < 2; ++i) delete[] _obj->a[i];
	for (size_t i = 0; i < n + 1; ++i) delete[] _obj->ders[i];

	delete[] _obj->a;
	delete[] _obj->ders;
	delete[] _obj->ndu;
	delete[] _obj->left;
	delete[] _obj->right;
	delete _obj;
}

DerivationInit* initDerivationInitStruct(const surfInfo* sfI, size_t der_degree)
{
	DerivationInit* res = new DerivationInit();
	res->der_degree = der_degree;
	res->SKL = new d_vertex * [der_degree + 1];
	for (size_t i = 0; i <= der_degree; ++i) res->SKL[i] = new d_vertex[der_degree + 1];
	res->temp = new d_vertex[sfI->q + 1];
	res->sU = initDersBasisFunsStruct(sfI->p, der_degree);
	res->sV = initDersBasisFunsStruct(sfI->q, der_degree);

	return res;
};

void releaseDerivationInitStruct(const surfInfo* sfI, DerivationInit* _obj)
{
	for (int i = 0; i <= _obj->der_degree; ++i) delete[] _obj->SKL[i];
	releaseDersBasisFunsStruct(sfI->p, _obj->der_degree, _obj->sU);
	releaseDersBasisFunsStruct(sfI->q, _obj->der_degree, _obj->sV);
	delete[] _obj->SKL, 
	delete[] _obj->temp;
	delete _obj;
}

void SurfaceDerivsAlg1(const surfInfo* sfI, double u, double v, DerivationInit* der_init)
{
	size_t d = der_init->der_degree;
	size_t du = min(sfI->p, d);
	size_t dv = min(sfI->q, d);
	for (size_t k = sfI->p + 1; k <= d; ++k)
		for (size_t l = 0; l <= d - k; ++l) der_init->SKL[k][l] = d_vertex();

	for (size_t l = sfI->q + 1; l <= d; ++l)
		for (size_t k = 0; k <= d - l; ++k) der_init->SKL[k][l] = d_vertex();

	size_t uspan = FindSpan(sfI->n, sfI->p, u, sfI->Uk, sfI->n + sfI->p + 1);
	size_t vspan = FindSpan(sfI->m, sfI->q, v, sfI->Vl, sfI->m + sfI->q + 1);

	DersBasisFuns(uspan, u, sfI->p, der_init->der_degree, sfI->Uk, der_init->sU);
	DersBasisFuns(vspan, v, sfI->q, der_init->der_degree, sfI->Vl, der_init->sV);

	double** Nu = der_init->sU->ders;
	double** Nv = der_init->sV->ders;
	
	for (size_t k = 0; k <= du; ++k)
	{
		for (size_t s = 0; s <= sfI->q; ++s)
		{
			der_init->temp[s] = d_vertex();
			for (size_t r = 0; r <= sfI->p; ++r)
			{
				der_init->temp[s] = der_init->temp[s] + Nu[k][r] * sfI->controlPoints[uspan - sfI->p + r][vspan - sfI->q + s];
			}
		}
		size_t dd = min(d - k, dv);
		for (size_t l = 0; l <= dd; ++l)
		{
			der_init->SKL[k][l] = d_vertex();
			for (size_t s = 0; s <= sfI->q; ++s)
			{
				der_init->SKL[k][l] = der_init->SKL[k][l] + Nv[l][s] * der_init->temp[s];
			}
		}
	}
}

//U - knot vector; ders - result; n (n<= p) - derivate degree
void DersBasisFuns(size_t i, double u, int p, int n, double* U, DersBasisFunsInit* der_bf)
{
	double** ders = der_bf->ders;
	double** ndu = der_bf->ndu;
	double** a = der_bf->a;

	double* left = der_bf->left;
	double* right = der_bf->right;
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
			else			j2 = p - r;

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
}
