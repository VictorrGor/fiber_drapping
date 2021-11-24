#pragma once
#pragma pack(1)
#include <d3d11.h>
#include <DirectXMath.h>

#include "Light\PointLight.h"

using namespace DirectX;

struct vertexCB
{
	DirectX::XMMATRIX mWVP;
	DirectX::XMMATRIX mWorld;
};


struct PS_perFrame_CB
{
	PointLight ll;
	XMFLOAT3 eyePos;
	float dummy;
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
