#pragma once
#include <d3d11.h>
#include "Vertex.h"

#include "ConstantBufers.h"

class Object
{
	ID3D11Buffer*	pVertexBuf;
	UINT			vecCount;
	ID3D11Buffer*	pIndexBuf;
	UINT			indexCount;

	ID3D11VertexShader* pVxSh;
	ID3D11PixelShader*	pPxSh;
	D3D_PRIMITIVE_TOPOLOGY	toplology;
	PixelShaderCB*	pCB;
	ID3D11Buffer*	pPS_CB;///depricated
public:
	Object(ID3D11Device* _pDevice, ID3D11VertexShader* _pVxSh, ID3D11PixelShader* _pPxSh, UINT _vertexCount, vertex* vecArr, D3D_PRIMITIVE_TOPOLOGY _toplology,
		ID3D11Buffer* _pVxBuf = nullptr);
	Object(ID3D11Device* _pDevice, ID3D11VertexShader* _pVxSh, ID3D11PixelShader* _pPxSh, UINT _vertexCount, vertex* vecArr, D3D_PRIMITIVE_TOPOLOGY _toplology,
		UINT _indexCount, UINT* indexArray, ID3D11Buffer* _pVxBuf = nullptr, ID3D11Buffer* _pIndexBuf = nullptr);
	~Object();

	void setMaterial(bool isUsed, XMFLOAT4 Ambient = XMFLOAT4(), XMFLOAT4 Diffuse = XMFLOAT4(), XMFLOAT4 Specular = XMFLOAT4(), XMFLOAT4 Reflect = XMFLOAT4());
	void setMaterialState(bool isUsed)
	{
		pCB->dummy[0] = isUsed;
	};
	void Render(ID3D11DeviceContext* _pDeviceContext, ID3D11Buffer* pPS_CB_per_obj, UINT _startLocation);
};



class TexturedObject
{
	ID3D11Buffer* pVertexBuf;
	ID3D11Buffer* pIndexBuf;
	UINT vecCount;
	UINT indexCount;
	D3D_PRIMITIVE_TOPOLOGY	toplology;

	ID3D11Texture2D* texture;
	ID3D11SamplerState* sampleState;

	ID3D11VertexShader* pVxSh;
	ID3D11PixelShader* pPxSh;

	///@todo Depricated. Change this old code to new one, which will use constant buffers
	struct PixelShader_textured
	{
		ID3D11Texture2D* texture;
		ID3D11SamplerState* sampleState;
	};
	PixelShader_textured* pCB;
	ID3D11Buffer* pPSBuff;
public:
	TexturedObject(ID3D11Device* _pDevice, ID3D11VertexShader* _pVxSh, ID3D11PixelShader* _pPxSh, UINT _vertexCount, 
		TexturedVertex* vecArr, D3D_PRIMITIVE_TOPOLOGY _toplology, ID3D11Buffer* _pVxBuf = nullptr);
	TexturedObject(ID3D11Device* _pDevice, ID3D11VertexShader* _pVxSh, ID3D11PixelShader* _pPxSh, UINT _vertexCount, 
		TexturedVertex* vecArr, D3D_PRIMITIVE_TOPOLOGY _toplology, UINT _indexCount, UINT* indexArray, 
		ID3D11Buffer* _pVxBuf = nullptr, ID3D11Buffer* _pIndexBuf = nullptr);

	void Render(ID3D11DeviceContext* _pDeviceContext);
};
