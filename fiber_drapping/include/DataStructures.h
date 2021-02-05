#pragma once

#include <DirectXMath.h>
#include <Windows.h>

#pragma pack(1)
using namespace DirectX;

typedef XMFLOAT3 vec3;
typedef XMFLOAT4 vec4;
typedef XMMATRIX mtx;


struct PointLight
{
	PointLight() { ZeroMemory(this, sizeof(this)); }
	XMFLOAT4 Ambient;
	XMFLOAT4 Diffuse;
	XMFLOAT4 Specular; // Packed into 4D vector: (Position, Range) 
	XMFLOAT3 Position;
	float Range; // Packed into 4D vector: (A0, A1, A2, Pad) 
	XMFLOAT3 Att;
	float Pad; // Pad the last float so we can set an // array of lights if we wanted. 
};

struct Light
{
	Light() { ZeroMemory(this, sizeof(this)); }
	XMFLOAT3 Direction;
	float pad;
	XMFLOAT4 Ambient;
	XMFLOAT4 Diffuse;
};

struct vertex
{
	vec3 pos;
	vec4 Color;
	vec3 normal;

	double getLength();
	vertex()
	{
		pos = vec3(0, 0, 0);
		Color = vec4(0, 0, 0, 0);
		normal = vec3(0, 0, 0);
	}

	vertex& operator=(const vertex& _vx);
	vertex& operator+=(const vertex& _vx);
	//vetrex& ope
	void makeAbs();
	vertex getAbs();
};

vec3 operator*(const vec3& _vx, double _m);
vec3 operator*(const vec3& _lvx, const vec3& _rvx);

vec3 operator+(const vec3& _vx1, const vec3& _vx2);
vec3 operator-(const vec3& _vx1, const vec3& _vx2);
vec3 operator/(const vec3& _vx, double _m);
vec3 sqrt(const vec3& _vx);

vertex operator-(const vertex& _vx1, const vertex& _vx2);
vertex operator+(const vertex& _vx1, const vertex& _vx2);
vertex operator*(const double& _multiplier, const vertex& _vx);
vertex operator*(const vertex& _vx, const double& _multiplier);
vertex operator/(const vertex& _vx, const double& _multiplier);
double getLength(const vec3& pos);

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

struct vertexCB
{
	DirectX::XMMATRIX mWVP;
	DirectX::XMMATRIX mWorld;
};


struct PS_perFrame_CB
{
	Light ll;
};

//Per object
struct PixelShaderCB
{
	XMFLOAT4 Ambient;
	XMFLOAT4 Diffuse;
	XMFLOAT4 Specular;
	XMFLOAT4 Reflect;
	//bool isUsesMaterial;
	bool dummy[16];// size of buffer must by devided by 16

};


vertex* convert2DimArrayTo1(vertex** mx, size_t n, size_t m);