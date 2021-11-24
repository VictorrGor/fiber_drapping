#pragma once

#include <DirectXMath.h>
#include <Windows.h>
#include "Vertex.h"
#pragma pack(1)
using namespace DirectX;

typedef XMMATRIX mtx;

struct vertex;

struct Light
{
	Light() { ZeroMemory(this, sizeof(Light)); }
	XMFLOAT3 Direction;
	float pad;
	XMFLOAT4 Ambient;
	XMFLOAT4 Diffuse;
};


struct splineInfo
{
	splineInfo();
	splineInfo(const splineInfo& _spI);
	splineInfo& operator=(const splineInfo& _spI);
	~splineInfo();

	vertex* controlPoints;
	size_t cpCount;
	double* knotVector;
	size_t knotLength;
	double* forwardU;
};

struct surfInfo
{	
	vertex** controlPoints;
	size_t n, m; //array size
	size_t p, q; //surface spline degree
	//knot vectors with length: pt_ct + degree
	double* Uk;//n+p
	double* Vl;//m+q

	surfInfo();
	~surfInfo();
	surfInfo(const surfInfo& obj);
	surfInfo& operator=(const surfInfo& obj);
	surfInfo(vertex** cp, size_t _n, size_t _m, size_t _p, size_t _q, double* _Uk, double* _Vl);
};


vertex* convert2DimArrayTo1(vertex** mx, UINT n, UINT m);

//Warp B Spline surface point. Does not release memory
struct bSplinePt
{
	vertex* pt;
	double u;
	double v;

	bSplinePt();
	bSplinePt(const bSplinePt& _obj);
	~bSplinePt();
	bSplinePt(vertex* _pt, double _u, double _v);
	bSplinePt& operator=(const bSplinePt& obj);
};



struct TextTextureInfo
{
	char ch;
	float leftU;
	float rightU;
	UINT pixelWidth;
};