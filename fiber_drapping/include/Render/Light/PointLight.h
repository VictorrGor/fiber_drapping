#pragma once
#include <DirectXMath.h>
#include <Windows.h>
#pragma pack(1)
using namespace DirectX;


struct PointLight
{
	PointLight() { ZeroMemory(this, sizeof(PointLight)); }
	XMFLOAT4 Ambient;
	XMFLOAT4 Diffuse;
	XMFLOAT4 Specular; // Packed into 4D vector: (Position, Range) 
	XMFLOAT3 Position;
	float Range; // Packed into 4D vector: (A0, A1, A2, Pad) 
};
