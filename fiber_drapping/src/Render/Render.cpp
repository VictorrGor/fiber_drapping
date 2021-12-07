// This is a personal academic project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: http://www.viva64.com
#include "Render/Render.h"

std::ofstream file;


XMVECTOR ComputeNormal(FXMVECTOR p0, FXMVECTOR p1, FXMVECTOR p2)
{
	XMVECTOR u = p1 - p0;
	XMVECTOR v = p2 - p0;
	return XMVector3Normalize(XMVector3Cross(u, v));
}

RenderSys::RenderSys()
{
	g_driverType = D3D_DRIVER_TYPE_NULL;
	g_featureLevel = D3D_FEATURE_LEVEL_11_0;
	g_pd3dDevice = NULL;
	g_pImmediateContext = NULL;
	g_pSwapChain = NULL;
	g_pRenderTargetView = NULL;
	mouse = nullptr;
	pInterface = nullptr;
}

///@todo Release ibject like ClenupDevice
RenderSys::~RenderSys()
{
	for (std::vector<Object*>::iterator it = objects.begin(); it < objects.end(); ++it)
	{
		delete (*it);
	}
	for (std::vector<PointLight*>::iterator it = pointLightObjects.begin(); it < pointLightObjects.end(); ++it)
	{
		delete (*it);
	}
}

void RenderSys::Render()
{
	//Draw objects
	g_pImmediateContext->IASetInputLayout(this->pVtxLayout);
	UINT stride = sizeof(vertex);
	UINT offset = 0;

	static float speed = 0;

	float R = 1. * mouse->getWheelPos() * 0.1;
	static float h = 1.3;

	PS_perFrame_CB perFrame;
	if (this->pointLightObjects.size())
	{
		perFrame.ll = *(this->pointLightObjects.back());
	}

	static float t = 0.0f;
	{
		/* 
		static DWORD dwTimeStart = 0;
		DWORD dwTimeCur = GetTickCount();
		if (dwTimeStart == 0) dwTimeStart = dwTimeCur;
		t = (dwTimeCur - dwTimeStart) / 1000.0f;
		speed = t / 2;
		*/

		static POINT deltaMouse;
		deltaMouse = mouse->updateSavedPos();
		static float speedCoeff = 0.01;

		speed += deltaMouse.x * speedCoeff;
		h += deltaMouse.y * speedCoeff;

		DirectX::XMVECTOR Eye = DirectX::XMVectorSet(R * cos(speed), h, R * sin(speed), 0.0f);
		DirectX::XMVECTOR At = DirectX::XMVectorSet(0.0f, 0.0f, 0.0f, 0.0f);
		DirectX::XMVECTOR Up = DirectX::XMVectorSet(0.0f, 1.0f, 0.0f, 0.0f);
		g_View = DirectX::XMMatrixLookAtLH(Eye, At, Up);
		XMStoreFloat3(&perFrame.eyePos, Eye);
	}

	g_pImmediateContext->UpdateSubresource(g_pPS_CB, 0, NULL, &perFrame, 0, 0);
	g_pImmediateContext->PSSetConstantBuffers(0, 1, &g_pPS_CB);

	//
	// Clear the back buffer
	//
	float ClearColor[4] = { 1.0f, 1.f, 1.f, 1.0f }; // red,green,blue,alpha
	g_pImmediateContext->ClearRenderTargetView(g_pRenderTargetView, ClearColor);

	vertexCB cb0;
	cb0.mWorld = DirectX::XMMatrixTranspose(g_World);
	cb0.mWVP = XMMatrixTranspose(g_World * g_View * g_Projection);

	g_pImmediateContext->VSSetShader(pVxSh, NULL, 0);
	g_pImmediateContext->UpdateSubresource(g_pConstantBuffer, 0, NULL, &cb0, 0, 0);
	g_pImmediateContext->VSSetConstantBuffers(0, 1, &g_pConstantBuffer);


	Object* pObj = nullptr;
	for (std::vector<Object*>::iterator it = objects.begin(); it < objects.end() - 1; ++it)
	{
		pObj = *it;
		pObj->Render(g_pImmediateContext, pPS_CB_per_obj, 0);
	}
	if (renderCount > 0)
	{
		(*(objects.end() - 1))->Render(g_pImmediateContext, pPS_CB_per_obj, renderCount);
	}
	else
	{
		(*(objects.end() - 1))->Render(g_pImmediateContext, pPS_CB_per_obj, 0);
	}
	//Draw Interface
	pInterface->Render(g_pImmediateContext);

	g_pSwapChain->Present(1, 0);
	g_pImmediateContext->ClearDepthStencilView(DepthBuffer, D3D11_CLEAR_DEPTH | D3D11_CLEAR_STENCIL, 1.0f, 0);
}

void RenderSys::CleanupDevice()
{
	if (g_pImmediateContext) g_pImmediateContext->ClearState();
	if (g_pConstantBuffer) g_pConstantBuffer->Release();
	if (g_pPS_CB) g_pPS_CB->Release();
	if (g_pRenderTargetView) g_pRenderTargetView->Release();
	if (g_pSwapChain) g_pSwapChain->Release();
	if (g_pImmediateContext) g_pImmediateContext->Release();
	if (g_pd3dDevice) g_pd3dDevice->Release();
	if (pPS_CB_per_obj) pPS_CB_per_obj->Release();
	if (mouse) delete mouse;
}

HRESULT RenderSys::InitShaders()
{
	HRESULT hRes = S_OK;

	////////////////////////////////
	//Compile shaders
	////////////////////////////////
	ID3DBlob* pBlobVertex = NULL;
	ID3DBlob* pPSerr = NULL;
	hRes = D3DCompileFromFile(L"shaders/VertexShader.hlsl", NULL, D3D_COMPILE_STANDARD_FILE_INCLUDE, "main", "vs_4_0", NULL, NULL, &pBlobVertex, &pPSerr);
	if (FAILED(hRes))
	{
		if (pPSerr) pPSerr->Release();
		MessageBox(NULL,
			L"The FX file cannot be compiled.  Please run this executable from the directory that contains the FX file.", L"Error", MB_OK);
		return hRes;
	}

	hRes = g_pd3dDevice->CreateVertexShader(pBlobVertex->GetBufferPointer(), pBlobVertex->GetBufferSize(), NULL, &pVxSh);


	// Define the input layout
	D3D11_INPUT_ELEMENT_DESC layout[] =
	{
		{ "POSITION", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 0, D3D11_INPUT_PER_VERTEX_DATA, 0 },
		{ "COLOR", 0, DXGI_FORMAT_R32G32B32A32_FLOAT, 0, 12, D3D11_INPUT_PER_VERTEX_DATA, 0 },
		{ "NORMAL", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 28, D3D11_INPUT_PER_VERTEX_DATA, 0 }

	};
	UINT numElements = ARRAYSIZE(layout);

	// Create the input layout
	hRes = g_pd3dDevice->CreateInputLayout(layout, numElements, pBlobVertex->GetBufferPointer(),
		pBlobVertex->GetBufferSize(), &pVtxLayout);
	pBlobVertex->Release();
	if (FAILED(hRes))
		return hRes;

	// Set the input layout
	g_pImmediateContext->IASetInputLayout(pVtxLayout);

	// Compile the pixel shader
	ID3DBlob* pPSBlob = NULL;
	hRes = D3DCompileFromFile(L"shaders/PixelShader.hlsl", NULL, D3D_COMPILE_STANDARD_FILE_INCLUDE, "main", "ps_4_0", NULL, NULL, &pPSBlob, &pPSerr);
	if (FAILED(hRes))
	{
		MessageBox(NULL,
			L"The FX file cannot be compiled.  Please run this executable from the directory that contains the FX file.", L"Error", MB_OK);
		return hRes;
	}

	// Create the pixel shader
	hRes = g_pd3dDevice->CreatePixelShader(pPSBlob->GetBufferPointer(), pPSBlob->GetBufferSize(), NULL, &pPxSh);
	pPSBlob->Release();
	if (FAILED(hRes))
		return hRes;

	if (pPSerr) pPSerr->Release();

	return hRes;
}

HRESULT RenderSys::InitConstantBuffers(UINT width, UINT height)
{
	HRESULT hRes = S_OK;
	//Initlize constant buffers
	D3D11_BUFFER_DESC bd;
	ZeroMemory(&bd, sizeof(bd));

	bd.Usage = D3D11_USAGE_DEFAULT;
	bd.ByteWidth = sizeof(vertexCB);
	bd.BindFlags = D3D11_BIND_CONSTANT_BUFFER;
	bd.CPUAccessFlags = 0;
	hRes = g_pd3dDevice->CreateBuffer(&bd, NULL, &g_pConstantBuffer);
	if (FAILED(hRes))
		return hRes;

	// Initialize the world matrix
	g_World = DirectX::XMMatrixIdentity();
	// Initialize the projection matrix
	g_Projection = DirectX::XMMatrixPerspectiveFovLH(DirectX::XM_PIDIV2, width / (FLOAT)height, 0.01f, 100.0f);

	D3D11_BUFFER_DESC bd1;
	ZeroMemory(&bd1, sizeof(bd1));

	bd1.Usage = D3D11_USAGE_DEFAULT;
	bd1.ByteWidth = sizeof(PS_perFrame_CB);
	int a = sizeof(XMFLOAT3);
	int b = sizeof(Light);
	bd1.BindFlags = D3D11_BIND_CONSTANT_BUFFER;
	bd1.CPUAccessFlags = 0;
	hRes = g_pd3dDevice->CreateBuffer(&bd1, NULL, &g_pPS_CB);
	if (FAILED(hRes))
		return hRes;

	Light light;
	light.Direction = XMFLOAT3(0.25f, 1, 0.50f);
	light.Ambient = XMFLOAT4(0.2f, 0.2f, 0.2f, 1.0f);
	light.Diffuse = XMFLOAT4(1.0f, 1.0f, 1.0f, 1.0f);

	//g_pImmediateContext->UpdateSubresource(g_pPS_CB, 0, NULL, &light, 0, 0);
	//g_pImmediateContext->PSSetConstantBuffers(0, 1, &g_pPS_CB);
	///////Per obj
	ZeroMemory(&bd1, sizeof(bd1));

	bd1.Usage = D3D11_USAGE_DEFAULT;
	bd1.ByteWidth = sizeof(PixelShaderCB);
	bd1.BindFlags = D3D11_BIND_CONSTANT_BUFFER;
	bd1.CPUAccessFlags = 0;
	hRes = g_pd3dDevice->CreateBuffer(&bd1, NULL, &pPS_CB_per_obj);
	if (FAILED(hRes))
		return hRes;

	return hRes;
}

HRESULT RenderSys::InitDevice(HWND* hWnd)
{
	int sizeP = sizeof(POINT);
	int xs = sizeof(LONG);
	mouse = new Mouse();
	int sss = sizeof(Mouse);
	HRESULT hRes = S_OK;

	RECT rc;
	GetClientRect(*hWnd, &rc);

	UINT width = rc.right - rc.left;
	UINT height = rc.bottom - rc.top;

	UINT createDeviceFlags = 0;

#ifdef _DEBUG
	createDeviceFlags |= D3D11_CREATE_DEVICE_DEBUG;
#endif


	D3D_DRIVER_TYPE driverTypes[] =
	{
		D3D_DRIVER_TYPE_HARDWARE,
		D3D_DRIVER_TYPE_WARP,
		D3D_DRIVER_TYPE_REFERENCE,
	};
	UINT numDriverTypes = ARRAYSIZE(driverTypes);

	D3D_FEATURE_LEVEL featureLevels[] =
	{
		D3D_FEATURE_LEVEL_11_0,
		D3D_FEATURE_LEVEL_10_1,
		D3D_FEATURE_LEVEL_10_0,
	};
	UINT numFeatureLevels = ARRAYSIZE(featureLevels);


	DXGI_SWAP_CHAIN_DESC sd;
	ZeroMemory(&sd, sizeof(sd));
	sd.BufferCount = 2;
	sd.BufferDesc.Width = width;
	sd.BufferDesc.Height = height;
	sd.BufferDesc.Format = DXGI_FORMAT_R8G8B8A8_UNORM;
	sd.BufferDesc.RefreshRate.Numerator = 60;
	sd.BufferDesc.RefreshRate.Denominator = 1;
	sd.BufferUsage = DXGI_USAGE_RENDER_TARGET_OUTPUT;
	sd.OutputWindow = *hWnd;
	sd.SampleDesc.Count = 1;
	sd.SampleDesc.Quality = 0;
	sd.Windowed = TRUE;


	for (UINT driverTypeIndex = 0; driverTypeIndex < numDriverTypes; driverTypeIndex++)
	{
		g_driverType = driverTypes[driverTypeIndex];
		hRes = D3D11CreateDeviceAndSwapChain(NULL, g_driverType, NULL, createDeviceFlags, featureLevels, numFeatureLevels,
			D3D11_SDK_VERSION, &sd, &g_pSwapChain, &g_pd3dDevice, &g_featureLevel, &g_pImmediateContext);
		if (SUCCEEDED(hRes))
			break;
	}
	if (FAILED(hRes))
		return hRes;

	// Create a render target view
	ID3D11Texture2D* pBackBuffer = NULL;
	hRes = g_pSwapChain->GetBuffer(0, __uuidof(ID3D11Texture2D), (LPVOID*)&pBackBuffer);
	if (FAILED(hRes))
		return hRes;

	hRes = g_pd3dDevice->CreateRenderTargetView(pBackBuffer, NULL, &g_pRenderTargetView);
	pBackBuffer->Release();
	if (FAILED(hRes))
		return hRes;

	D3D11_TEXTURE2D_DESC depthTextureDesc;
	ZeroMemory(&depthTextureDesc, sizeof(depthTextureDesc));
	depthTextureDesc.Width = width;
	depthTextureDesc.Height = height;
	depthTextureDesc.MipLevels = 1;
	depthTextureDesc.ArraySize = 1;
	depthTextureDesc.SampleDesc.Count = 1;
	depthTextureDesc.Format = DXGI_FORMAT_D24_UNORM_S8_UINT;
	depthTextureDesc.BindFlags = D3D11_BIND_DEPTH_STENCIL;

	ID3D11Texture2D* DepthStencilTexture;
	hRes = g_pd3dDevice->CreateTexture2D(&depthTextureDesc, NULL, &DepthStencilTexture);

	if (FAILED(hRes))
		return hRes;

	D3D11_DEPTH_STENCIL_VIEW_DESC dsvDesc;
	ZeroMemory(&dsvDesc, sizeof(dsvDesc));
	dsvDesc.Format = depthTextureDesc.Format;
	dsvDesc.ViewDimension = D3D11_DSV_DIMENSION_TEXTURE2DMS;

	hRes = g_pd3dDevice->CreateDepthStencilView(DepthStencilTexture, &dsvDesc, &DepthBuffer);
	DepthStencilTexture->Release();

	if (FAILED(hRes))
		return hRes;

	g_pImmediateContext->OMSetRenderTargets(1, &g_pRenderTargetView, DepthBuffer);

	// Setup the viewport
	D3D11_VIEWPORT vp;
	vp.Width = (FLOAT)width;
	vp.Height = (FLOAT)height;
	vp.MinDepth = 0.0f;
	vp.MaxDepth = 1.0f;
	vp.TopLeftX = 0;
	vp.TopLeftY = 0;
	g_pImmediateContext->RSSetViewports(1, &vp);

	//Init Interface
	this->pInterface = new Interface(g_pd3dDevice, width, height);

	hRes = InitShaders();
	if (FAILED(hRes))
		return hRes;
	//Create vertex bufffers for B-spline and Big points
	hRes = InitObjects();
	if (FAILED(hRes))
		return hRes;
	hRes = InitConstantBuffers(width, height);

	return hRes;
}


void RenderSys::pushObject(Object* _obj)
{
	this->objects.push_back(_obj);
}



HRESULT RenderSys::InitObjects()
{
	HRESULT hRes = S_OK;
	//Generate Light points
	PointLight* pl = new PointLight();
	pl->Position = vec3(0.0, 3, 0);
	pl->Range = 20.;
	pl->Ambient = vec4(0.5, 0.5, 0.5, 1);
	pl->Diffuse = vec4(1., 1., 1., 1);
	pl->Specular = vec4(1, 1., 1, 1);
	this->pointLightObjects.push_back(pl);


	//testPlane(this, this->g_pd3dDevice, pVxSh, pPxSh);
	generateSurfaceByBSpline(this, this->g_pd3dDevice, pVxSh, pPxSh);
	//drawSinSurf(this, this->g_pd3dDevice, pVxSh, pPxSh);
	//drawSinSurf(this, g_pd3dDevice, pVxSh, pPxSh);
	//generateSphere(this, g_pd3dDevice, pVxSh, pPxSh);
	//lighting_test(this, g_pd3dDevice, pVxSh, pPxSh);
	mp.transferToRender(this);
	return hRes;
}

HRESULT RenderSys::drawTriangle(vertex* _pt)
{
	HRESULT hRes = S_OK;
	mp.addNew(_pt, 3);
	return hRes;

}



void RenderSys::drawLineOnBSplineSurface(surfInfo* sfi, double u, double v, bool isU)
{
	UINT size = 500;
	vertex* line = new vertex[size];
	for (UINT i = 0; i < size; ++i)
	{
		if (isU) line[i] = SurfacePoint(sfi, u, i * 1. / (size - 1));
		else line[i] = SurfacePoint(sfi, i * 1. / (size - 1), v);
		line[i].Color = vec4(1, 0, 0, 1);
	}
	Object* obj = new Object(g_pd3dDevice, pVxSh, pPxSh, size, line, D3D11_PRIMITIVE_TOPOLOGY_LINESTRIP);
	objects.push_back(obj);
}

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

void RenderSys::disableMaterials()
{
	for (std::vector<Object*>::iterator it = objects.begin(); it < objects.end(); ++it)
	{
		(*it)->setMaterialState(false);
	}
}

void RenderSys::OnResize(UINT _width, UINT _height)
{
	width = _width;
	height = _height;
	//@TODO Call resize for all object, that need resizing
	if (pInterface) pInterface->Resize(_width, _height);
}


Mouse* RenderSys::getMouse()
{
	return mouse;
}


Interface* RenderSys::getInterface()
{
	return pInterface;
}

ID3D11Device* RenderSys::getDevice()
{
	return g_pd3dDevice;
}

int RenderSys::getRenderCount()
{
	return renderCount;
}


void RenderSys::incrementRenderObjCount()
{
	renderCount++;
}
void RenderSys::decrementRenderObjCount()
{
	renderCount--;
}