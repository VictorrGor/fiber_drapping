#pragma once

#include <DirectXMath.h>
#include <Render/Vertex.h>
using namespace DirectX;

typedef XMMATRIX mtx;
typedef XMFLOAT3 vec3;

struct d_vertex {
	double x;
	double y;
	double z;
	d_vertex();
	d_vertex(double _x, double _y, double _z);
	d_vertex& operator+=(const d_vertex& _obj);
};

d_vertex operator*(const double& _s, const d_vertex& _obj);
d_vertex operator*(const d_vertex& _obj, const double& _s);
d_vertex operator+(const d_vertex& _obj1, const d_vertex& _obj2);
d_vertex operator-(const d_vertex& _obj1, const d_vertex& _obj2);
d_vertex operator/(const d_vertex& _obj, const double& _s);


double getDistance(const d_vertex& _obj1, const d_vertex& _obj2);


///Save result of substraction as vec3 of float for using fast SSE comands
vec3 subtructAsVec3(const d_vertex& _obj1, const d_vertex& _obj2);

struct splineInfo
{
	splineInfo();
	splineInfo(const splineInfo& _spI);
	splineInfo& operator=(const splineInfo& _spI);
	~splineInfo();

	d_vertex* controlPoints;
	size_t cpCount;
	double* knotVector;
	size_t knotLength;
	double* forwardU;
};

struct surfInfo
{
	d_vertex** controlPoints;
	size_t n, m; //array size
	size_t p, q; //surface spline degree
	//knot vectors with length: pt_ct + degree
	double* Uk;//n+p
	double* Vl;//m+q

	surfInfo();
	~surfInfo();
	surfInfo(const surfInfo& obj);
	surfInfo& operator=(const surfInfo& obj);
	surfInfo(d_vertex** cp, size_t _n, size_t _m, size_t _p, size_t _q, double* _Uk, double* _Vl);
};


//Warp B Spline surface point. Does not release memory
struct bSplinePt
{
	d_vertex* pt;
	double u;
	double v;

	bSplinePt();
	bSplinePt(const bSplinePt& _obj);
	~bSplinePt();
	bSplinePt(d_vertex* _pt, double _u, double _v);
	bSplinePt& operator=(const bSplinePt& obj);
};

