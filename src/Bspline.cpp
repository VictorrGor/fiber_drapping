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
#ifdef LOG_ON
	Log(vtx, pointCount);
#endif
	vertex* derivatePoints = getDerivatePoints(pointCount, vtx);

	size_t knotVectorSize = (pointCount - 2) + 2 * (splineDegree + 1);
	double* forwardU = new double[pointCount];				//Индексы для уравнения Qk = Sum(0..n)[Ni,p(forwardU)*Pi]
	double* backwardU = new double[knotVectorSize];//Индексы для построения полученного спалйна
	getInterpolationKnotVector(pointCount, forwardU, backwardU, knotVectorSize, vtx);

	//////////////Output//////////////////////////
#ifdef LOG_ON
	Log("Forward vector is: ");
	for (int i = 0; i < pointCount; ++i)
	{
		sprintf(buf, "%lf ", forwardU[i]);
		Log(buf);
	}
	Log("\n");

	Log("Knot vector is: ");
	for (int i = 0; i < knotVectorSize; ++i)
	{
		sprintf(buf, "%lf ", backwardU[i]);
		Log(buf);
	}
	Log("\n");
#endif 
	//////////////Output//////////////////////////

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

#ifdef LOG_ON
	Log(P, pointCount + 2);
#endif

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

#ifdef LOG_ON
	for (size_t i = 0; i < size; ++i)
	{
		file << knots[i] << "\n";
	}
#endif

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
#ifdef LOG_ON
	std::ofstream pFile;
	pFile.open("out.txt");
	if (!pFile.good()) MessageBox(NULL, L"Can't open log file!", L"Error!", 0);
	pFile << "x\ty\tz\n";
#endif
	vertex* result = new vertex[_vxCount];

	double* _knots = _spI.knotVector;
	double dt = _spI.knotVector[_p - 1];
	double step = (_knots[_spI.cpCount] - _knots[_p - 1]) / (_vxCount - 1);
	double Nq = 0;

	for (size_t j = 0; j < _vxCount; ++j, dt += step)
	{
		result[j] = CurvePoint(_spI, _p, dt);
	}
#ifdef LOG_ON
	pFile.close();
#endif
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

double getSplineLen(double _left, double _right, vertex *(*ffunc)(splineInfo, size_t, double, size_t), splineInfo _spi, size_t _p, size_t n)
{
	vec3 v1 = integrate(0, 1, *CurveDerivateAlg1, _spi);
	return sqrt(pow(v1.x, 2) + pow(v1.y, 2) + pow(v1.z, 2));
}
