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
//Вернуть в addSpline уборщик мусора!!
// Вернуть в addInterpolationSpline уборщик мусора!!
//Починить функцию рисования больших точек

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
	UINT vecCount;

	ID3D11Buffer*			pIndexBuf;
	UINT					indexCount;

	ID3D11VertexShader*		pVxSh;
	ID3D11PixelShader*		pPxSh;
	D3D_PRIMITIVE_TOPOLOGY	toplology;
	PixelShaderCB*			pCB;

	ID3D11Buffer*			pPS_CB;///depricated

	RenderSys* rs;
public:
	friend RenderSys;
	Object(RenderSys* _rs, ID3D11VertexShader* _pVxSh, ID3D11PixelShader* _pPxSh, UINT _vertexCount, vertex* vecArr, D3D_PRIMITIVE_TOPOLOGY _toplology,
		ID3D11Buffer* _pVxBuf = nullptr);
	Object(RenderSys* _rs, ID3D11VertexShader* _pVxSh, ID3D11PixelShader* _pPxSh, UINT _vertexCount, vertex* vecArr, D3D_PRIMITIVE_TOPOLOGY _toplology,
		UINT _indexCount, UINT* indexArray, ID3D11Buffer* _pVxBuf = nullptr, ID3D11Buffer* _pIndexBuf = nullptr);
	~Object();

	void setMaterial(bool isUsed, XMFLOAT4 Ambient = XMFLOAT4(), XMFLOAT4 Diffuse = XMFLOAT4(), XMFLOAT4 Specular = XMFLOAT4(), XMFLOAT4 Reflect = XMFLOAT4());
	void setMaterialState(bool isUsed)
	{
		pCB->dummy[0] = isUsed;
	};
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
	std::vector<PointLight*> pointLightObjects;

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


			Object* obj = new Object(rs, rs->pVxSh, rs->pPxSh, objCount, mem, D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST);
			//obj->setMaterial(true,  vec4(0.2, 0.2, 0.2, 1), vec4(0.3, 0.5, 0.3, 1), 
			//	vec4(0.2, 0.2, 0.2, 1),	vec4(0.2, 0.2, 0.2, 1)); 
			rs->objects.push_back(obj);

			if (nextPool) nextPool->transferToRender(rs);
		}
	};
	MemoryPool mp;

	Mouse* mouse;
public:
	friend Object;

	RenderSys();
	~RenderSys();
	// Render the frame
	void Render();
	// Clean up the objects we've created
	void CleanupDevice();

	HRESULT InitShaders();

	HRESULT InitConstantBuffers(UINT width, UINT height);

	HRESULT InitDevice(HWND* hWnd);

	//All lines must start at once point
	void drapping_part(surfInfo* sfI, double u1, double v1, bool isU1, 
									 double u2, double v2, bool isU2)
	{
		UINT size = GIRD_SIZE;
		bSplinePt** P = new bSplinePt * [size]; //points warper
		
		vertex** Q = new vertex * [size];
		for (UINT i = 0; i < size; ++i)
		{
			P[i] = new bSplinePt[size];

			Q[i] = new vertex[size];
			for (UINT j = 0; j < size; ++j)
			{
				Q[i][j].Color = vec4(0, 1, 0, 1);
				P[i][j].pt = &Q[i][j];
			}
		}

		UINT dim = 2;
		double** W = new double* [dim]; //Jakobian
		double** invW = new double* [dim]; //inverse
		for (UINT i = 0; i < dim; ++i)
		{
			W[i] = new double[dim];
			invW[i] = new double[dim];
		}
		UINT err_ct = 0;

		double cycle_step = 1. / (size - 1);

		//Generating initial lines
		for (UINT i = 0; i < size; ++i)
		{
			if (isU1)
			{
				P[0][i].u = u1;
				P[0][i].v = i * cycle_step;
			}
			else
			{
				P[0][i].u = i * cycle_step;
				P[0][i].v = v1;
			}
			if (isU2)
			{ 
				P[i][0].u = u2;
				P[i][0].v = i * cycle_step;
			}
			else
			{
				P[i][0].u = i * cycle_step;
				P[i][0].v = v2;
			}
			Q[0][i] = SurfacePoint(sfI, P[0][i].u, P[0][i].v);
			Q[i][0] = SurfacePoint(sfI, P[i][0].u, P[i][0].v);
		}

		double delta_u = 0.01;
		
		for (UINT i = 0; i <  size - 1; ++i)
		{
			for (UINT j = 0; j < size - 1/*size - 1*/; ++j)
			{
				if ((P[i][j].u < 0) || (P[i][j + 1].u < 0) || (P[i + 1][j].u < 0))
					continue;

				bSplinePt* ptIJ = &P[i + 1][j + 1];
				bSplinePt* ptIm1J = &P[i][j + 1];
				bSplinePt* ptIJm1 = &P[i + 1][j];

				ptIJ->u = ptIJm1->u;
				ptIJ->v = ptIm1J->v;// ptIJm1->v - 2* delta_u;
				if ((ptIJ->u == ptIJm1->u) && (ptIJ->v == ptIJm1->v) || (ptIJ->u == ptIm1J->u) && (ptIJ->v == ptIm1J->v))
				{
					ptIJ->u = (ptIJm1->u + ptIm1J->u) / 2;
				}
				if (ptIJ->u > 1) ptIJ->u -= 2 * delta_u;
				if (ptIJ->v > 1) ptIJ->v -= 2 * delta_u;
				if (ptIJ->u < 0) ptIJ->u = (ptIJm1->u + ptIm1J->u) /2;//0;
				if (ptIJ->v < 0) ptIJ->v = min(ptIJm1->v, ptIm1J->v);// 0;// 1 + ptIJ->v;

				(*ptIJ->pt) = SurfacePoint(sfI, ptIJ->u, ptIJ->v);

#ifdef _DEBUG
				std::cout << "\ti:" << i << "; j:" << j << "\n";
#endif

				if ((i == 0) && (j == 1))
				{
					vertex* trace;
					int trace_size;
					if (!getBSplineDrapPoint_with_trace(W, invW, ptIJ, ptIm1J, ptIJm1, sfI, &trace, trace_size))
					{
						ptIJ->u = -1;
						ptIJ->v = -1;
						++err_ct;
						std::cout << err_ct << "\n";
						std::cout << "\ti:" << i << "; j:" << j << "\n";
					}
					for (int i = 0; i < trace_size; ++i)
					{
						trace[i].Color = vec4(0, 1, i * (1. / (trace_size - 1)), 1.);
					}
					Object* trace_obj = new Object(this, pVxSh, pPxSh, trace_size, trace, D3D11_PRIMITIVE_TOPOLOGY_LINESTRIP);
					this->objects.push_back(trace_obj);
				}
				if (!getBSplineDrapPoint(W, invW, ptIJ, ptIm1J, ptIJm1, sfI))
				{
					ptIJ->u = -1;
					ptIJ->v = -1;
					++err_ct;
					std::cout << err_ct << "\n";
					std::cout << "\ti:" << i << "; j:" << j << "\n";
				}
#ifdef _DEBUG
				std::cout << "\tu:" << ptIJ->u << "; v:" << ptIJ->v << "\n";
#endif
				ptIJ->pt->Color = vec4(0, 1, 0, 1);
			}
		}
		std::cout << "\nDrapping errors: " << err_ct <<"\n";
		Object* obj = new Object(this, pVxSh, pPxSh, size * size /*size * size*/, convert2DimArrayTo1(Q, size, size), D3D11_PRIMITIVE_TOPOLOGY_POINTLIST);
		//objects.push_back(obj);
		Object* obj1 = new Object(this, pVxSh, pPxSh, size * size /*size * size*/, convert2DimArrayTo1(Q, size, size), D3D11_PRIMITIVE_TOPOLOGY_LINESTRIP);
		//objects.push_back(obj1);


		vertex* triangle = new vertex[3];
		for (UINT i = 0; i < size - 1; ++i)
		{
			for (UINT j = 0; j < size - 1; ++j)
			{
				if ((P[i][j].u < 0) || (P[i][j + 1].u < 0) || (P[i + 1][j].u < 0))
					continue;
				float angle = getAngle(Q, i, j, i + 1, j, i, j + 1);
				float red, green, blue, coeff;
				coeff = 1;
				if (angle < 90)
					red = 1 - angle * coeff / 90;
				else
					red = 0;
				if (angle > 90)
				{
					green = 1 - (angle - 90) * coeff / 90;
					blue = (angle - 90) * coeff / 90;
				}
				else
				{
					green = angle * coeff / 90;
					blue = 0;
				}

				triangle[0] = Q[i][j];
				triangle[1] = Q[i][j + 1];
				triangle[2] = Q[i + 1][j];
				triangle[0].normal = calculateTriangleNormal(Q[i][j], Q[i][j + 1], Q[i + 1][j]);
				triangle[1].normal = triangle[0].normal;
				triangle[2].normal = triangle[0].normal;

				triangle[0].Color = vec4(red, green, blue, 1);
				triangle[1].Color = vec4(red, green, blue, 1);
				triangle[2].Color = vec4(red, green, blue, 1);
				drawTriangle(triangle);
			}
		}
		for (UINT i = 0; i < size - 1; ++i)
		{
			for (UINT j = 0; j < size - 1; ++j)
			{
				float angle = getAngle(Q, i + 1, j + 1, i + 1, j, i, j + 1);
				float red, green, blue, coeff;
				coeff = 1;
				if (angle < 90)
					red = 1 - angle * coeff / 90;
				else
					red = 0;
				if (angle > 90)
				{
					green = 1 - (angle - 90) * coeff / 90;
					blue = (angle - 90) * coeff / 90;
				}
				else
				{
					green = angle * coeff / 90;
					blue = 0;
				}


				vertex* triangle = new vertex[3];
				if ((P[i+1][j+1].u < 0) || (P[i][j + 1].u < 0) || (P[i + 1][j].u < 0))
					continue;
				triangle[0] = Q[i+1][j];
				triangle[1] = Q[i][j + 1];
				triangle[2] = Q[i + 1][j + 1];
				triangle[0].normal = calculateTriangleNormal(Q[i+1][j], Q[i][j + 1], Q[i + 1][j + 1]);
				triangle[1].normal = triangle[0].normal;
				triangle[2].normal = triangle[0].normal;

				triangle[0].Color = vec4(red, green, blue, 1);
				triangle[1].Color = vec4(red, green, blue, 1);
				triangle[2].Color = vec4(red, green, blue, 1);
				drawTriangle(triangle);
			}
		}


		for (UINT i = 0; i < size; ++i)
		{
			delete[] Q[i];
			delete[] P[i];
		}
		for (UINT i = 0; i < dim; ++i)
		{
			delete W[i];
			delete invW[i];
		}
		delete[] Q;
		delete[] P;
		delete[] W;
		delete[] invW;
	}

	void test_surface()
	{
		UINT n = 10;
		UINT m = 10;
		vertex** Q = new vertex * [n];
		for (UINT i = 0; i < n; ++i) Q[i] = new vertex[m];

		double step_fi = XM_PI * 2 / (n - 1);
		double step_teta = XM_PIDIV2 / (m - 1);
		double R = 1;

		for (UINT i = 0; i < n; ++i)
			for (UINT j = 0; j < m; ++j)
				Q[i][j].pos = vec3(R * cos(i * step_fi) * cos(step_teta * j), R * sin(step_teta * j) , R * cos(j * step_teta) * sin(step_fi * i));


		surfInfo sfI = GenInterpBSplineSurface(n, m, Q, 3, 3);


		UINT size = 100;
		vertex* res = new vertex[size * size];
		

		for(UINT i = 0; i < size; ++i)
			for (UINT j = 0; j < size; ++j)
			{
				res[i * size + j] = SurfacePoint(&sfI, 1. / (size - 1) * i, 1. / (size - 1) * j);
				res[i * size + j].Color = vec4(0, 0, 0, 1);
			}
		
		Object* obj = new Object(this, pVxSh, pPxSh, size * size, res, D3D11_PRIMITIVE_TOPOLOGY_POINTLIST);
		//Object* obj = new Object(this, pVxSh, pPxSh, n*m, convert2DimArrayTo1(Q, n, n), D3D11_PRIMITIVE_TOPOLOGY_POINTLIST);
		//objects.push_back(obj);

		//Drapping part
		drawLineOnBSplineSurface(&sfI, 0, 0, false);
		drawLineOnBSplineSurface(&sfI, 0.25, 0, true);
		//drawLineOnBSplineSurface(&sfI, 0.5, 0, true);
		//drawLineOnBSplineSurface(&sfI, 0.75, 0, true);
		
		drapping_part(&sfI, 0, 0, true, 0.25, 0, true);
		drapping_part(&sfI, 0.25, 0, true, 0.5, 0, true);
		drapping_part(&sfI, 0.5, 0, true, 0.75, 0, true);
		drapping_part(&sfI, 0.75, 0, true, 1, 0, true);
	}

	void generateSphere();

	HRESULT InitObjects();
	
	void coonsLines_new();
	
	void lighting_test();

	HRESULT drawSpline(splineInfo _info);

	HRESULT drawTriangle(vertex* _pt);

	Mouse* getMouse();
	vertex** drawSinSurf();

	//u, v - parametric coordinates. isU - is true if u - const coordinate. 
	void drawLineOnBSplineSurface(surfInfo* sfi, double u, double v, bool isU);


	vec3 calculateTriangleNormal(const vertex& v0, const vertex& v1, const vertex& v2)
	{
		XMFLOAT3 u = v1.pos - v0.pos;
		XMFLOAT3 v = v2.pos - v0.pos;
		XMVECTOR u_sse = XMLoadFloat3(&u);
		XMVECTOR v_sse = XMLoadFloat3(&v);
		XMFLOAT3 normal;
		XMStoreFloat3(&normal, XMVector3Cross(u_sse, v_sse));
		return normal;
	};

	void disableMaterials()
	{
		for (std::vector<Object*>::iterator it = objects.begin(); it < objects.end(); ++it)
		{ 
			(*it)->setMaterialState(false);
		}
	}
};
