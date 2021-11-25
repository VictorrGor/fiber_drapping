#pragma once

#include <DirectXMath.h>
#include <Windows.h>
#include "Vertex.h"
#pragma pack(1)
using namespace DirectX;


struct vertex;

struct Light
{
	Light() { ZeroMemory(this, sizeof(Light)); }
	XMFLOAT3 Direction;
	float pad;
	XMFLOAT4 Ambient;
	XMFLOAT4 Diffuse;
};

vertex* convert2DimArrayTo1(vertex** mx, UINT n, UINT m);

struct TextTextureInfo
{
	char ch;
	float leftU;
	float rightU;
	UINT pixelWidth;
};