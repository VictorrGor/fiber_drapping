// This is a personal academic project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: http://www.viva64.com
#include "../Render/Vertex.h"


vertex::vertex(const vertex& _obj)
{
	pos = _obj.pos;
	Color = _obj.Color;
	normal = _obj.normal;
}

vertex& vertex::operator=(const vertex& _vx)
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

vec3 operator*(const vec3& _vx, double _m)
{
	return vec3(_vx.x * _m, _vx.y * _m, _vx.z * _m);
}

vec3 operator*(const vec3& _lvx, const vec3& _rvx)
{
	return vec3(_lvx.x * _rvx.x, _lvx.y * _rvx.y, _lvx.z * _rvx.z);
}

vec3 operator+(const vec3& _vx1, const vec3& _vx2)
{
	vec3 res = vec3(_vx1.x + _vx2.x, _vx1.y + _vx2.y, _vx1.z + _vx2.z);
	return res;
}

vec3 operator-(const vec3& _vx1, const vec3& _vx2)
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

float getDistance(const vec3& _obj1, const vec3& _obj2)
{
	return sqrtf(powf((_obj1.x - _obj2.x), 2) + powf((_obj1.y - _obj2.y), 2) + powf((_obj1.z - _obj2.z), 2));
}
