#include "DataStructures.h"


splineInfo::splineInfo()
{
	controlPoints = nullptr;
	knotVector = nullptr;
	forwardU = nullptr;
	cpCount = 0;
	knotLength = 0;
}

splineInfo::splineInfo(const splineInfo & _spI)
{
	cpCount = _spI.cpCount;
	knotLength = _spI.knotLength;
	controlPoints = new vertex[cpCount];
	knotVector = new double[knotLength];
	forwardU = new double[cpCount];

	for (size_t i = 0; i < cpCount; ++i)
	{
		controlPoints[i] = _spI.controlPoints[i];
		forwardU[i] = _spI.forwardU[i];
	}
	for (size_t i = 0; i < knotLength; ++i)
	{
		knotVector[i] = _spI.knotVector[i];
	}
}

splineInfo& splineInfo::operator=(const splineInfo & _spI)
{

	if (controlPoints) delete[] controlPoints;
	if (knotVector) delete[] knotVector;
	if (forwardU) delete[] forwardU;

	cpCount = _spI.cpCount;
	knotLength = _spI.knotLength;
	controlPoints = new vertex[cpCount];
	knotVector = new double[knotLength];
	forwardU = new double[cpCount];

	for (size_t i = 0; i < cpCount; ++i)
	{
		controlPoints[i] = _spI.controlPoints[i];
		forwardU[i] = _spI.forwardU[i];
	}
	for (size_t i = 0; i < knotLength; ++i)
	{
		knotVector[i] = _spI.knotVector[i];
	}

	return *this;
}

splineInfo::~splineInfo()
{
	if (controlPoints) delete[] controlPoints;
	if (knotVector) delete[] knotVector;
	if (forwardU) delete[] forwardU;
	
	controlPoints = nullptr;
	forwardU = nullptr;
	knotVector = nullptr;

	cpCount = 0;
	knotLength = 0;
}

vec3 operator*(const vec3& _vx, double _m)
{
	return vec3(_vx.x * _m, _vx.y * _m, _vx.z * _m);
}

vec3 operator*(const vec3 & _lvx, const vec3 & _rvx)
{
	return vec3(_lvx.x * _rvx.x, _lvx.y * _rvx.y, _lvx.z * _rvx.z);
}

vec3 operator+(const vec3& _vx1, const vec3& _vx2)
{
	vec3 res = vec3(_vx1.x + _vx2.x, _vx1.y + _vx2.y, _vx1.z + _vx2.z);
	return res;
}

vec3 operator-(const vec3 & _vx1, const vec3 & _vx2)
{
	return vec3(_vx1.x - _vx2.x, _vx1.y - _vx2.y, _vx1.z - _vx2.z);
}

vec3 operator/(const vec3& _vx, double _m)
{
	return vec3(_vx.x * _m, _vx.y * _m, _vx.z * _m);
}

vec3 sqrt(const vec3& _vx)
{
	return vec3(sqrt(_vx.x), sqrt(_vx.y), sqrt(_vx.z));
}

double vertex::getLength()
{
	return sqrt(pow(pos.x, 2) + pow(pos.y, 2) + pow(pos.z, 2));
}

double getLength(const vec3& pos)
{
	return sqrt(pow(pos.x, 2) + pow(pos.y, 2) + pow(pos.z, 2));
}



vertex & vertex::operator=(const vertex & _vx)
{
	pos = _vx.pos;
	Color = _vx.Color;
	normal = _vx.normal;
	return *this;
}

vertex& vertex::operator+=(const vertex& _vx)
{
	pos = pos + _vx.pos;
	return *this;
}

void vertex::makeAbs()
{
	pos.x = fabs(pos.x);
	pos.y = fabs(pos.y);
	pos.z = fabs(pos.z);
}

vertex vertex::getAbs()
{
	vertex res(*this);
	res.makeAbs();
	return res;
}

vertex operator-(const vertex& _vx1, const vertex& _vx2)
{
	vertex result;
	result.Color = _vx1.Color;
	result.pos.x = _vx1.pos.x - _vx2.pos.x;
	result.pos.y = _vx1.pos.y - _vx2.pos.y;
	result.pos.z = _vx1.pos.z - _vx2.pos.z;

	return result;
}
vertex operator+(const vertex& _vx1, const vertex& _vx2)
{
	vertex result;
	result.Color = _vx1.Color;
	result.pos.x = _vx1.pos.x + _vx2.pos.x;
	result.pos.y = _vx1.pos.y + _vx2.pos.y;
	result.pos.z = _vx1.pos.z + _vx2.pos.z;

	return result;
}
vertex operator*(const double& _multiplier, const vertex& _vx)
{
	vertex res = _vx;
	res.pos.x *= _multiplier;
	res.pos.y *= _multiplier;
	res.pos.z *= _multiplier;
	return res;
}
vertex operator*(const vertex& _vx, const double& _multiplier)
{
	vertex res = _vx;
	res.pos.x *= _multiplier;
	res.pos.y *= _multiplier;
	res.pos.z *= _multiplier;
	return res;
}
vertex operator/(const vertex& _vx, const double& _multiplier)
{
	vertex res = _vx;
	res.pos.x /= _multiplier;
	res.pos.y /= _multiplier;
	res.pos.z /= _multiplier;
	return res;
}

vertex* convert2DimArrayTo1(vertex** mx, UINT n, UINT m)
{
	vertex* res = new vertex[n * m];
	for (UINT i = 0; i < n; ++i)
		for (UINT j = 0; j < m; ++j)
			res[i * n + j] = mx[i][j];

	return res;
}


surfInfo::surfInfo()
{
	this->controlPoints = nullptr;
	this->Uk = nullptr;
	this->Vl = nullptr;
}

surfInfo::~surfInfo()
{
	if (controlPoints)
	{
		for (size_t i = 0; i < n; ++i)
		{
			if (controlPoints[i]) delete[] controlPoints[i];
		}
		delete[] controlPoints;
	}
	if (Uk) delete Uk;
	if (Vl) delete Vl;
}

surfInfo::surfInfo(const surfInfo& obj)
{
	this->controlPoints = nullptr;
	this->Uk = nullptr;
	this->Vl = nullptr;

	this->n = obj.n;
	this->m = obj.m;
	this->p = obj.p;
	this->q = obj.q;
	if (obj.controlPoints)
	{
		this->controlPoints = new vertex*[this->n];
		for (size_t i = 0; i < this->n; ++i)
		{
			if (obj.controlPoints[i])
			{
				this->controlPoints[i] = new vertex[this->m];
				for (size_t j = 0; j < this->m; ++j) this->controlPoints[i][j] = obj.controlPoints[i][j];
			}
			else throw "surfaceInfo: incorrect obj!";
		}
		if (obj.Uk)
		{
			this->Uk = new double[this->n + this->p];
			for (size_t i = 0; i < n + p; ++i) this->Uk[i] = obj.Uk[i];
		}
		else throw "surfaceInfo: incorrect obj!";
		if (obj.Vl)
		{
			this->Vl = new double[this->m + this->q];
			for (size_t i = 0; i < this->m + this->q; ++i) this->Vl[i] = obj.Vl[i];
		}
		else throw "surfaceInfo: incorrect obj!";
	}
}

surfInfo& surfInfo::operator=(const surfInfo& obj)
{
	if (this->controlPoints)
	{
		for (size_t i = 0; i < this->n; ++i)
		{
			delete[] this->controlPoints[i];
		}
		delete[] this->controlPoints;
	}
	if (this->Uk) delete[] this->Uk;
	if (this->Vl) delete[] this->Vl;


	this->controlPoints = nullptr;
	this->Uk = nullptr;
	this->Vl = nullptr;
	this->n = obj.n;
	this->m = obj.m;
	this->p = obj.p;
	this->q = obj.q;
	if (obj.controlPoints)
	{
		this->controlPoints = new vertex * [this->n];
		for (size_t i = 0; i < this->n; ++i)
		{
			if (obj.controlPoints[i])
			{
				this->controlPoints[i] = new vertex[this->m];
				for (size_t j = 0; j < this->m; ++j) this->controlPoints[i][j] = obj.controlPoints[i][j];
			}
			else throw "surfaceInfo: incorrect obj!";
		}
		if (obj.Uk)
		{
			this->Uk = new double[this->n + this->p];
			for (size_t i = 0; i < n + p; ++i) this->Uk[i] = obj.Uk[i];
		}
		else throw "surfaceInfo: incorrect obj!";
		if (obj.Vl)
		{
			this->Vl = new double[this->m + this->q];
			for (size_t i = 0; i < this->m + this->q; ++i) this->Vl[i] = obj.Vl[i];
		}
		else throw "surfaceInfo: incorrect obj!";
	}
	return *this;
}

surfInfo::surfInfo(vertex** cp, size_t _n, size_t _m, size_t _p, size_t _q, double* _Uk, double* _Vl) : controlPoints(cp),
	n(_n), m(_m), p(_p), q(_q), Uk(_Uk), Vl(_Vl)
{
}

bSplinePt::bSplinePt()
{
	this->pt = nullptr;
	this->u = -1;
	this->v = -1;
}

bSplinePt::~bSplinePt()
{
	this->pt = nullptr;
}

bSplinePt::bSplinePt(vertex* _pt, double _u, double _v): pt(_pt), u(_u), v(_v)
{
}

bSplinePt& bSplinePt::operator=(const bSplinePt& obj)
{
	this->pt = obj.pt;
	this->u = obj.u;
	this->v = obj.v;
	return *this;
}
