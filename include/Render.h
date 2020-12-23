#pragma once
#include <Windows.h>
#include <d3d11.h>
#include <d3dcompiler.h>
#include <wchar.h>
#include <vector>
#include <winerror.h>
#include <iostream>
#include <fstream>
#include "MathLib.h"
#include "Bspline.h"
#include "DataStructures.h"
#include "Drapping.h"



extern std::ofstream file;

//TODO
//������� � addSpline ������� ������!!
// ������� � addInterpolationSpline ������� ������!!
//�������� ������� ��������� ������� �����

using namespace DirectX;

XMVECTOR ComputeNormal(FXMVECTOR p0, FXMVECTOR p1, FXMVECTOR p2);

class RenderSys;



void Log(const char* _str);
void Log(const double* _arr, size_t _size);
void Log(const vertex* _arr, size_t _size);



//Left-handed^ (xOz surface) and Y-up

class Object
{
	ID3D11Buffer*			pVertexBuf;
	size_t vecCount;

	ID3D11Buffer*			pIndexBuf;
	size_t					indexCount;

	ID3D11VertexShader*		pVxSh;
	ID3D11PixelShader*		pPxSh;
	D3D_PRIMITIVE_TOPOLOGY	toplology;
	PixelShaderCB*			pCB;

	ID3D11Buffer*			pPS_CB;///depricated

	RenderSys* rs;
public:
	friend RenderSys;
	Object(RenderSys* _rs, ID3D11VertexShader* _pVxSh, ID3D11PixelShader* _pPxSh, size_t _vertexCount, vertex* vecArr, D3D_PRIMITIVE_TOPOLOGY _toplology,
		ID3D11Buffer* _pVxBuf = nullptr);
	Object(RenderSys* _rs, ID3D11VertexShader* _pVxSh, ID3D11PixelShader* _pPxSh, size_t _vertexCount, vertex* vecArr, D3D_PRIMITIVE_TOPOLOGY _toplology,
		size_t _indexCount, size_t* indexArray, ID3D11Buffer* _pVxBuf = nullptr, ID3D11Buffer* _pIndexBuf = nullptr);
	~Object();

	void setMaterial(bool isUsed, XMFLOAT4 Ambient = XMFLOAT4(), XMFLOAT4 Diffuse = XMFLOAT4(), XMFLOAT4 Specular = XMFLOAT4(), XMFLOAT4 Reflect = XMFLOAT4());
};



class Mouse
{
	bool isLeftKeyPressed;
	bool isRightKeyPressed;
	int wheel_pos;

	POINT mousePos;
	POINT savedMPos;
public:
	friend RenderSys;

	Mouse();
	void updWheelPos(int newPos);
	void updLK(bool isPressed);
	void updRK(bool isPressed);
	void updMousePos(POINT mPos);
};

class RenderSys
{
	D3D_DRIVER_TYPE         g_driverType;
	D3D_FEATURE_LEVEL       g_featureLevel;
	ID3D11Device*           g_pd3dDevice;
	ID3D11DeviceContext*    g_pImmediateContext;
	IDXGISwapChain*         g_pSwapChain;
	ID3D11RenderTargetView* g_pRenderTargetView;

	ID3D11VertexShader*		pVxSh;
	ID3D11PixelShader*		pPxSh;

	ID3D11InputLayout*		pVtxLayout;
	ID3D11Buffer*			g_pVertexBuffer;
	ID3D11Buffer*			g_pIndexBuffer;

	ID3D11Buffer*			g_pConstantBuffer;
	ID3D11Buffer*			g_pPS_CB;
	ID3D11Buffer*			pPS_CB_per_obj;

	DirectX::XMMATRIX		g_World;
	DirectX::XMMATRIX		g_View;
	DirectX::XMMATRIX		g_Projection;

	std::vector<Object*>	objects;

	struct MemoryPool
	{
		vertex* mem;
		size_t objCount;
		size_t maxSize = 3000000;
		MemoryPool* nextPool;

		MemoryPool()
		{
			objCount = 0;
			mem = new vertex[maxSize];
			memset(mem, 0, sizeof(vertex) *  maxSize);
			nextPool = nullptr;
		}
		~MemoryPool()
		{
			delete[] mem;
			mem = nullptr;
			objCount = 0;
			if (nextPool)
			{
				delete nextPool;
				nextPool = nullptr;
			}
		}
		void addNew(vertex* _pt)
		{
			if (objCount == maxSize)
			{
				if (!nextPool) nextPool = new MemoryPool();
				nextPool->addNew(_pt);
			}
			else
			{
				mem[objCount] = *_pt;
				++objCount;
			}
		}
		void addNew(vertex* _pt, size_t _count)
		{
			if (objCount + _count >= maxSize)
			{
				if (!nextPool) nextPool = new MemoryPool();
				
				size_t usedThis = maxSize - objCount - 1;
				nextPool->addNew((_pt + usedThis), _count - usedThis);

				for (size_t i = 0; i < usedThis; ++i)
				{
					mem[objCount] = *(_pt+i);
					++objCount;
				}
			}
			else
			{
				for (size_t i = 0; i < _count; ++i)
				{
					mem[objCount] = *(_pt+i);
					++objCount;
				}
			}
		}
		void clear()
		{
			memset(mem, 0, sizeof(vertex) *  maxSize);
			objCount = 0;
			if (nextPool) nextPool->clear();
		}
		
		void transferToRender(RenderSys *rs)
		{
			/*XMVECTOR xm;
			for (int i = 0; i < objCount / 3; ++i)
			{
				xm = ComputeNormal(XMLoadFloat3(&mem[i*3].pos), XMLoadFloat3(&mem[i * 3 + 1].pos), XMLoadFloat3(&mem[i * 3 + 2].pos));
				for (int j = 0; j < 3; ++j) XMStoreFloat3(&mem[i * 3 + j].normal, xm);
			}*/


			Object* obj = new Object(rs, rs->pVxSh, rs->pPxSh, objCount, mem, D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST);
			obj->setMaterial(true,  vec4(0.2, 0.2, 0.2, 1), vec4(0.3, 0.5, 0.3, 1), 
				vec4(0.2, 0.2, 0.2, 1),	vec4(0.2, 0.2, 0.2, 1)); 
			rs->objects.push_back(obj);

			if (nextPool) nextPool->transferToRender(rs);
		}
	};
	MemoryPool mp;

	Mouse* mouse;
public:
	friend Object;

	RenderSys();
	// Render the frame
	void Render();
	// Clean up the objects we've created
	void CleanupDevice();

	HRESULT InitShaders();

	HRESULT InitConstantBuffers(UINT width, UINT height);

	HRESULT InitDevice(HWND* hWnd);
	
	void testObject()
	{
		//vertex* vx = new vertex[9];
		//vertex vvx;
		//vvx.Color = vec4(0, 1.0, 0, 1.0);
		//vvx.pos = vec3(0.5, 0, -0.5);
		//vx[0] = vvx;
		//vvx.pos = vec3(0, 0.5, 0);
		//vx[1] = vvx;
		//vvx.pos = vec3(0.5, 0, 0.5);
		//vx[2] = vvx;

		////vvx.Color = vec4(1, 0.0, 0, 1.0);
		//vvx.pos = vec3(0.5, 0, 0.5);
		//vx[3] = vvx;
		//vvx.pos = vec3(0, 0.5, 0);
		//vx[4] = vvx;
		//vvx.pos = vec3(0, 0, 0);
		//vx[5] = vvx;

		////vvx.Color = vec4(0, 0, 1, 1.0);
		//vvx.pos = vec3(0, 0, 0);
		//vx[6] = vvx;
		//vvx.pos = vec3(0, 0.5, 0);
		//vx[7] = vvx;
		//vvx.pos = vec3(0.5, 0, -0.5);
		//vx[8] = vvx;

		//mp.addNew(vx, 9);

		vertex* sphere = new vertex[10000];
		size_t ptCt = 100;
		float teta = 0;
		float fi = 0;
		float dTeta = DirectX::XM_PIDIV2 / (ptCt - 1);
		float dFi = DirectX::XM_2PI / (ptCt - 1) ;

		for (size_t i = 0; i < ptCt; i++)
		{
			teta = 0;
			for (size_t j = 0; j < ptCt; j++)
			{
				sphere[j + i * ptCt].pos.x = 1 * sin(teta) * cos(fi);
				sphere[j + i * ptCt].pos.y = 1 * cos(teta);
				sphere[j + i * ptCt].pos.z = 1 * sin(teta) * sin(fi);
				sphere[j + i * ptCt].Color = vec4(0.2, 0.2, 0.2, 1);
				teta += dTeta;
			}
			fi += dFi;
		}
		
		size_t arS = (ptCt - 1) * 2 * 3 * ptCt;
		//size_t* indices = new size_t[arS];
		vertex* vs = new vertex[arS];
		//memset(indices, 0, arS);
		size_t actualIndex = 0;

		//size_t i = 0;
		vec3 pos;
		for (size_t i = 0; i < ptCt - 1; i++)
		{
			for (size_t j = 0; j < ptCt - 1; j++)
			{
				vs[actualIndex] = sphere[i * ptCt + j];
				vs[actualIndex + 1] = sphere[(i + 1) * ptCt + j];
				vs[actualIndex + 2] = sphere[i * ptCt + j + 1];
				actualIndex += 3;
			}
		}
		for (size_t i = 0; i < ptCt - 1; i++)
		{
			for (size_t j = 0; j < ptCt - 1; j++)
			{
				vs[actualIndex + 1] = sphere[i * ptCt + j + 1];
				vs[actualIndex + 2] = sphere[(i + 1) * ptCt + j];
				vs[actualIndex ] = sphere[(i+1) * ptCt + j + 1];
				actualIndex += 3;
			}
		}

		//for (size_t j = 0; j < ptCt - 1; j++)
		//{
		//	vs[actualIndex] = (ptCt - 1) * ptCt + j;
		//	vs[actualIndex + 1] = (ptCt - 1) * ptCt + j;
		//	vs[actualIndex + 2] = (ptCt - 1) * ptCt + j + 1;
		//	actualIndex += 3;
		//}


		Object* obj = new Object(this, pVxSh, pPxSh, actualIndex, vs, D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST);// , actualIndex, indices);
		obj->setMaterial(true, vec4(1.2, 0.2, 0.2, 1), vec4(0.3, 0.5, 0.3, 1),
			vec4(0.2, 0.2, 0.2, 1), vec4(0.2, 0.2, 0.2, 1));
		objects.push_back(obj);
	}

	void drawDrappingPoints(vertex** points);

	HRESULT InitObjects();
	
	void coonsLines_new();
	
	HRESULT drawSpline(splineInfo _info);

	HRESULT drawTriangle(vertex* _pt);

	Mouse* getMouse();
};