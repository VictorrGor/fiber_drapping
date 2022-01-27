#pragma once
#include <Windows.h>
#include <winerror.h>
#include <d3d11.h>
#include <d3dcompiler.h>
#include <iostream>
#include <fstream>
#include <wchar.h>
#include <vector>

#include "Math/MathLib.h"
#include "Math/Bspline.h"
#include "Math/Drapping.h"

#include "Input/Mouse.h"

#include "RenderStructures.h"
#include "DDSTextureLoader.h"
#include "Object.h"
#include "Interface.h"
#include "Scene.h"

extern std::ofstream file;

using namespace DirectX;

XMVECTOR ComputeNormal(FXMVECTOR p0, FXMVECTOR p1, FXMVECTOR p2);
vec3 calculateTriangleNormal(const vertex& v0, const vertex& v1, const vertex& v2);

class Mouse;
class RenderSys;


//Left-handed^ (xOz surface) and Y-up

class Camera
{
	vec3 lookAt;
	vec3 position;
	vec3 up;
	float speed; //move speed
	float rotationSpeed;
public:
	friend RenderSys;
	Camera() : lookAt({ 1,0,0 }), position({ 0,0,0 }), up({0, 1, 0}), speed(0.05), rotationSpeed(0.01) {};
	void moveForward();
	void moveBackward();
	void moveLeft();
	void moveRight();
	void moveUp();
	void moveDown();

	void rotateCamera(float XoZAngle, float YoZAngle);
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

	ID3D11Buffer*			g_pConstantBuffer;	//Vertex shader constant buffer
	ID3D11Buffer*			g_pPS_CB;			//Pixel shader per frame constant buffer
	ID3D11Buffer*			pPS_CB_per_obj;		//Pixel shader per object constant buffer
	ID3D11DepthStencilView* DepthBuffer;


	///@todo Word matrix might be in object description
	DirectX::XMMATRIX		g_World;
	DirectX::XMMATRIX		g_View;
	DirectX::XMMATRIX		g_Projection;

	std::vector<Object*>	objects;
	std::vector<PointLight*> pointLightObjects;


	///@todo depricated usage
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
				
				nextPool->addNew(_pt, _count);
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


			Object* obj = new Object(rs->g_pd3dDevice, rs->pVxSh, rs->pPxSh, objCount, mem, D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST);
			//obj->setMaterial(true,  vec4(0.2, 0.2, 0.2, 1), vec4(0.3, 0.5, 0.3, 1), 
			//	vec4(0.2, 0.2, 0.2, 1),	vec4(0.2, 0.2, 0.2, 1)); 
			rs->objects.push_back(obj);

			if (nextPool) nextPool->transferToRender(rs);
		}
	};
	MemoryPool mp;

	Mouse* mouse;
	Interface* pInterface;
	
	int renderCount;

public:
	UINT width;
	UINT height;
	Camera mCamera;

	RenderSys();
	~RenderSys();
	// Render the frame
	void Render();
	// Clean up the objects we've created
	void CleanupDevice();
	HRESULT InitShaders();
	HRESULT InitConstantBuffers(UINT width, UINT height);
	HRESULT InitDevice(HWND* hWnd);
	HRESULT InitObjects();

	HRESULT drawTriangle(vertex* _pt);
	Mouse* getMouse();
	void disableMaterials();
	void OnResize(UINT _width, UINT _height);
	void pushObject(Object* _obj);
	
	Interface* getInterface();
	ID3D11Device* getDevice();
	
	//u, v - parametric coordinates. isU - is true if u - const coordinate. 
	void drawLineOnBSplineSurface(surfInfo* sfi, double u, double v, bool isU);

	//Debug functions, wich provide draw only given count of object.
	int getRenderCount();

	void incrementRenderObjCount();
	void decrementRenderObjCount();
};
