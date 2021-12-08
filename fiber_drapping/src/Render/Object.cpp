// This is a personal academic project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: http://www.viva64.com
#include "Render/Object.h"

Object::Object(ID3D11Device* _pDevice, ID3D11VertexShader* _pVxSh, ID3D11PixelShader* _pPxSh,
	UINT _vertexCount, vertex* vecArr, D3D_PRIMITIVE_TOPOLOGY _toplology,
	ID3D11Buffer* _pVxBuf) : pVertexBuf(_pVxBuf), vecCount(_vertexCount),
	pVxSh(_pVxSh), pPxSh(_pPxSh), toplology(_toplology), pIndexBuf(nullptr), indexCount(0)
{
	if (!_pVxBuf)
	{

		D3D11_BUFFER_DESC bufferDesc;
		bufferDesc.Usage = D3D11_USAGE_DEFAULT;
		bufferDesc.ByteWidth = sizeof(vertex) * _vertexCount;
		bufferDesc.BindFlags = D3D11_BIND_VERTEX_BUFFER;
		bufferDesc.CPUAccessFlags = 0;
		bufferDesc.MiscFlags = 0;

		D3D11_SUBRESOURCE_DATA InitData;
		ZeroMemory(&InitData, sizeof(InitData));
		InitData.pSysMem = vecArr;
		_pDevice->CreateBuffer(&bufferDesc, &InitData, &pVertexBuf);


		D3D11_BUFFER_DESC bd;
		ZeroMemory(&bd, sizeof(bd));
		bd.Usage = D3D11_USAGE_DEFAULT;
		bd.ByteWidth = sizeof(PixelShaderCB);
		bd.BindFlags = D3D11_BIND_CONSTANT_BUFFER;
		bd.CPUAccessFlags = 0;

		HRESULT hRes = _pDevice->CreateBuffer(&bd, NULL, &pPS_CB);

		if (FAILED(hRes))
		{
#ifdef DEBUG_CONSOLE
			std::cout << "[Error] Cann't create Constant Buffer for Object!" << std::endl;
#endif
			pPS_CB = nullptr;
		}
		else
		{
			pCB = new PixelShaderCB();
			setMaterial(false);
		}
	}
}

Object::Object(ID3D11Device* _pDevice, ID3D11VertexShader* _pVxSh, ID3D11PixelShader* _pPxSh,
	UINT _vertexCount, vertex* vecArr, D3D_PRIMITIVE_TOPOLOGY _toplology,
	UINT _indexCount, UINT* indexArray, ID3D11Buffer* _pVxBuf, ID3D11Buffer* _pIndexBuf) : pVertexBuf(_pVxBuf), vecCount(_vertexCount),
	pVxSh(_pVxSh), pPxSh(_pPxSh), toplology(_toplology), pIndexBuf(nullptr), indexCount(_indexCount)
{
	if (!_pVxBuf)
	{

		D3D11_BUFFER_DESC bufferDesc;
		ZeroMemory(&bufferDesc, sizeof(bufferDesc));
		bufferDesc.Usage = D3D11_USAGE_DEFAULT;
		bufferDesc.ByteWidth = sizeof(vertex) * _vertexCount;
		bufferDesc.BindFlags = D3D11_BIND_VERTEX_BUFFER;
		bufferDesc.CPUAccessFlags = 0;
		bufferDesc.MiscFlags = 0;

		D3D11_SUBRESOURCE_DATA InitData;
		ZeroMemory(&InitData, sizeof(InitData));
		InitData.pSysMem = vecArr;
		HRESULT res = _pDevice->CreateBuffer(&bufferDesc, &InitData, &pVertexBuf);


		D3D11_BUFFER_DESC bd;
		ZeroMemory(&bd, sizeof(bd));
		bd.Usage = D3D11_USAGE_DEFAULT;
		bd.ByteWidth = sizeof(PixelShaderCB);
		bd.BindFlags = D3D11_BIND_CONSTANT_BUFFER;
		bd.CPUAccessFlags = 0;

		HRESULT hRes = _pDevice->CreateBuffer(&bd, NULL, &pPS_CB);

		if (FAILED(hRes))
		{
#ifdef DEBUG_CONSOLE
			std::cout << "[Error] Cann't create Constant Buffer for Object!" << std::endl;
#endif
			pPS_CB = nullptr;
		}
		else
		{
			pCB = new PixelShaderCB();
			setMaterial(false);
		}
	}
	if (!_pIndexBuf && _indexCount)
	{
		D3D11_BUFFER_DESC bd;
		ZeroMemory(&bd, sizeof(bd));
		bd.Usage = D3D11_USAGE_DEFAULT;
		bd.ByteWidth = sizeof(UINT) * indexCount;
		bd.BindFlags = D3D11_BIND_INDEX_BUFFER;
		bd.CPUAccessFlags = 0;
		bd.MiscFlags = 0;

		D3D11_SUBRESOURCE_DATA InitData;
		ZeroMemory(&InitData, sizeof(InitData));
		InitData.pSysMem = indexArray;
		_pDevice->CreateBuffer(&bd, &InitData, &pIndexBuf);
	}
}

Object::~Object()
{
	if (pVertexBuf)
	{
		pVertexBuf->Release();
		pVertexBuf = nullptr;
	}
	if (pIndexBuf)
	{
		pIndexBuf->Release();
		pIndexBuf = nullptr;
	}
	if (pPS_CB)
	{
		pPS_CB->Release();
		pPS_CB = nullptr;
	}
	if (pCB)
	{
		delete pCB;
		pCB = nullptr;
	}
}

void Object::setMaterial(bool isUsed, XMFLOAT4 Ambient, XMFLOAT4 Diffuse, XMFLOAT4 Specular, XMFLOAT4 Reflect)
{
	pCB->dummy[0] = isUsed;
	pCB->Ambient = Ambient;
	pCB->Diffuse = Diffuse;
	pCB->Specular = Specular;
}

void Object::Render(ID3D11DeviceContext* _pDeviceContext, ID3D11Buffer* pPS_CB_per_obj, UINT _startLocation)
{
	UINT stride = sizeof(vertex);
	UINT offset = 0;
	_startLocation *= 3;

	_pDeviceContext->IASetPrimitiveTopology(toplology);
	_pDeviceContext->IASetVertexBuffers(0, 1, &pVertexBuf, &stride, &offset);
	_pDeviceContext->VSSetShader(pVxSh, NULL, 0);
	_pDeviceContext->PSSetShader(pPxSh, NULL, 0);
	_pDeviceContext->UpdateSubresource(pPS_CB_per_obj, 0, NULL, pCB, 0, 0);
	_pDeviceContext->PSSetConstantBuffers(1, 1, &pPS_CB_per_obj);

	if (_startLocation > vecCount) _startLocation = vecCount;
	if (pIndexBuf)
	{
		_pDeviceContext->IASetIndexBuffer(pIndexBuf, DXGI_FORMAT_R32_UINT, offset);
		_pDeviceContext->DrawIndexed(indexCount - _startLocation, _startLocation, 0);
	}
	else
	{
		//_pDeviceContext->Draw(vecCount - _startLocation, _startLocation);
		if (!_startLocation)
		{
			_pDeviceContext->Draw(vecCount - _startLocation, _startLocation);
		}
		else
		{
			_pDeviceContext->Draw(_startLocation, 0);
		}
	}
}



TexturedObject::TexturedObject(ID3D11Device* _pDevice, ID3D11VertexShader* _pVxSh, ID3D11PixelShader* _pPxSh, 
	UINT _vertexCount, TexturedVertex* vecArr, D3D_PRIMITIVE_TOPOLOGY _toplology, ID3D11Buffer* _pVxBuf) : pVertexBuf(_pVxBuf), vecCount(_vertexCount),
	pVxSh(_pVxSh), pPxSh(_pPxSh), toplology(_toplology)
{
	if (!_pVxBuf)
	{

		D3D11_BUFFER_DESC bufferDesc;
		bufferDesc.Usage = D3D11_USAGE_DEFAULT;
		bufferDesc.ByteWidth = sizeof(TexturedVertex) * _vertexCount;
		bufferDesc.BindFlags = D3D11_BIND_VERTEX_BUFFER;
		bufferDesc.CPUAccessFlags = 0;
		bufferDesc.MiscFlags = 0;

		D3D11_SUBRESOURCE_DATA InitData;
		ZeroMemory(&InitData, sizeof(InitData));
		InitData.pSysMem = vecArr;
		_pDevice->CreateBuffer(&bufferDesc, &InitData, &pVertexBuf);
	}
}

TexturedObject::TexturedObject(ID3D11Device* _pDevice, ID3D11VertexShader* _pVxSh, ID3D11PixelShader* _pPxSh, 
	UINT _vertexCount, TexturedVertex* vecArr, D3D_PRIMITIVE_TOPOLOGY _toplology, UINT _indexCount, 
	UINT* _indexArray, ID3D11Buffer* _pVxBuf, ID3D11Buffer* _pIndexBuf) : pVertexBuf(_pVxBuf), vecCount(_vertexCount),
	pVxSh(_pVxSh), pPxSh(_pPxSh), toplology(_toplology), indexCount(_indexCount)
{
	if (!_pVxBuf)
	{

		D3D11_BUFFER_DESC bufferDesc;
		bufferDesc.Usage = D3D11_USAGE_DEFAULT;
		bufferDesc.ByteWidth = sizeof(TexturedVertex) * _vertexCount;
		bufferDesc.BindFlags = D3D11_BIND_VERTEX_BUFFER;
		bufferDesc.CPUAccessFlags = 0;
		bufferDesc.MiscFlags = 0;

		D3D11_SUBRESOURCE_DATA InitData;
		ZeroMemory(&InitData, sizeof(InitData));
		InitData.pSysMem = vecArr;
		_pDevice->CreateBuffer(&bufferDesc, &InitData, &pVertexBuf);


		if (!pIndexBuf && _indexArray && indexCount)
		{
			D3D11_BUFFER_DESC bd_indexes;
			ZeroMemory(&bd_indexes, sizeof(bd_indexes));
			bd_indexes.Usage = D3D11_USAGE_DEFAULT;
			bd_indexes.ByteWidth = sizeof(UINT) * indexCount;
			bd_indexes.BindFlags = D3D11_BIND_INDEX_BUFFER;
			bd_indexes.CPUAccessFlags = 0;

			D3D11_SUBRESOURCE_DATA InitData;
			ZeroMemory(&InitData, sizeof(InitData));
			InitData.pSysMem = _indexArray;
			//_pDevice->CreateBuffer(&bd_indexes, &InitData, &pIndexBuf);
			HRESULT hRes = _pDevice->CreateBuffer(&bd_indexes, &InitData, &pIndexBuf);

			if (FAILED(hRes))
			{
#ifdef DEBUG_CONSOLE
				std::cout << "[Error] Cann't create Index Buffer for Object!" << std::endl;
#endif
				pIndexBuf = nullptr;
			}

		}
	}
}

TexturedObject::TexturedObject()
{
	pVertexBuf = nullptr;
	pIndexBuf = nullptr;
	texture = nullptr;
	sampleState = nullptr;
}

TexturedObject::~TexturedObject()
{
	if (pVertexBuf)
	{
		pVertexBuf->Release();
		pVertexBuf = nullptr;
	}
	if (pIndexBuf)
	{
		pIndexBuf->Release();
		pIndexBuf = nullptr;
	}
}

void TexturedObject::Render(ID3D11DeviceContext* _pDeviceContext)
{
	UINT stride = sizeof(TexturedVertex);
	UINT offset = 0;
	_pDeviceContext->IASetPrimitiveTopology(toplology);
	_pDeviceContext->IASetVertexBuffers(0, 1, &pVertexBuf, &stride, &offset);

	if (pIndexBuf)
	{
		_pDeviceContext->IASetIndexBuffer(pIndexBuf, DXGI_FORMAT_R32_UINT, offset);
		_pDeviceContext->DrawIndexed(indexCount, 0, 0);
	}
	else
		_pDeviceContext->Draw(vecCount, 0);
}
