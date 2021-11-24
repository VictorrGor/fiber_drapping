// This is a personal academic project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: http://www.viva64.com
#include "Render/DataStructures.h"


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


vertex* convert2DimArrayTo1(vertex** mx, UINT n, UINT m)
{
	vertex* res = new vertex[n * m];
	for (UINT i = 0; i < n; ++i)
		for (UINT j = 0; j < m; ++j)
			res[i * m + j] = mx[i][j];

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

bSplinePt::bSplinePt(const bSplinePt& _obj)
{
	
	this->pt = _obj.pt;
	this->u = _obj.u;
	this->v = _obj.v;
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
