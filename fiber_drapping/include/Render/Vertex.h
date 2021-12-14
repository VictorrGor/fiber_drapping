#pragma once
#include <d3d11.h>
#include <DirectXMath.h>
#include "Math/MathStructures.h"
#pragma pack(1)

using namespace DirectX;

typedef XMFLOAT3 vec3;
typedef XMFLOAT4 vec4;

class d_vertex;

struct vertex
{
	vec3 pos;
	vec4 Color;
	vec3 normal;

	double getLength();
	vertex();
	vertex(const vertex& _obj);
	vertex(const d_vertex& _obj);
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

float getDistance(const vec3& _obj1, const vec3& _obj2);

struct TexturedVertex
{
	vertex vx;
	XMFLOAT2 tcoord;
};