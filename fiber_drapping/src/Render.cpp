#include "Render.h"

std::ofstream file;


XMVECTOR ComputeNormal(FXMVECTOR p0, FXMVECTOR p1, FXMVECTOR p2)
{
	XMVECTOR u = p1 - p0;
	XMVECTOR v = p2 - p0;
	return XMVector3Normalize(XMVector3Cross(u, v));
}



void Log(const char * _str)
{
	file << _str;// << "\n";
}

void Log(const double* _arr, size_t _size)
{
	for (size_t i = 0; i < _size; ++i)
	{
		file << _arr[i] << " ";
	}
	file << "\n";
}
void Log(const vertex* _arr, size_t _size)
{
	for (size_t i = 0; i < _size; ++i)
	{
		file << _arr[i].pos.x << " " << _arr[i].pos.y << " " << _arr[i].pos.z << "\n";
	}
	file << "\n";
}

Object::Object(RenderSys * _rs, ID3D11VertexShader * _pVxSh, ID3D11PixelShader * _pPxSh, 
	size_t _vertexCount, vertex * vecArr, D3D_PRIMITIVE_TOPOLOGY _toplology,
				ID3D11Buffer * _pVxBuf) : pVertexBuf(_pVxBuf), vecCount(_vertexCount),
				pVxSh(_pVxSh), pPxSh(_pPxSh), toplology(_toplology), pIndexBuf(nullptr), indexCount(0), rs(_rs)
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
		_rs->g_pd3dDevice->CreateBuffer(&bufferDesc, &InitData, &pVertexBuf);


		D3D11_BUFFER_DESC bd;
		ZeroMemory(&bd, sizeof(bd));
		bd.Usage = D3D11_USAGE_DEFAULT;
		bd.ByteWidth = sizeof(PixelShaderCB);
		bd.BindFlags = D3D11_BIND_CONSTANT_BUFFER;
		bd.CPUAccessFlags = 0;

		HRESULT hRes = rs->g_pd3dDevice->CreateBuffer(&bd, NULL, &pPS_CB);

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

Object::Object(RenderSys* _rs, ID3D11VertexShader* _pVxSh, ID3D11PixelShader* _pPxSh, size_t _vertexCount, vertex* vecArr, D3D_PRIMITIVE_TOPOLOGY _toplology,
	size_t _indexCount, size_t* indexArray, ID3D11Buffer* _pVxBuf, ID3D11Buffer* _pIndexBuf) : pVertexBuf(_pVxBuf), vecCount(_vertexCount),
	pVxSh(_pVxSh), pPxSh(_pPxSh), toplology(_toplology), pIndexBuf(nullptr), indexCount(_indexCount), rs(_rs)
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
		_rs->g_pd3dDevice->CreateBuffer(&bufferDesc, &InitData, &pVertexBuf);


		D3D11_BUFFER_DESC bd;
		ZeroMemory(&bd, sizeof(bd));
		bd.Usage = D3D11_USAGE_DEFAULT;
		bd.ByteWidth = sizeof(PixelShaderCB);
		bd.BindFlags = D3D11_BIND_CONSTANT_BUFFER;
		bd.CPUAccessFlags = 0;

		HRESULT hRes = rs->g_pd3dDevice->CreateBuffer(&bd, NULL, &pPS_CB);

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
		bd.ByteWidth = sizeof(size_t) * indexCount;
		bd.BindFlags = D3D11_BIND_INDEX_BUFFER;
		bd.CPUAccessFlags = 0;
		bd.MiscFlags = 0;

		D3D11_SUBRESOURCE_DATA InitData;
		ZeroMemory(&InitData, sizeof(InitData));
		InitData.pSysMem = indexArray;
		_rs->g_pd3dDevice->CreateBuffer(&bd, &InitData, &pIndexBuf);
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
	/*ps_pO.dummy[0] = true;
	ps_pO.Ambient = vec4(0.2, 0.2, 0.2, 1);
	ps_pO.Diffuse = vec4(0.3, 0.5, 0.3, 1);
	ps_pO.Specular = vec4(0.2, 0.2, 0.2, 1);
	ps_pO.Reflect = vec4(0.2, 0.2, 0.2, 1);

	ps_pO.Ambient = vec4(1, 0., 0., 1);
	ps_pO.Diffuse = vec4(0., 1, 0., 1);
	ps_pO.Specular = vec4(0., 0, 1, 1);
	ps_pO.Reflect = vec4(1, 0., 0., 1);*/

}

RenderSys::RenderSys()
{
	D3D_DRIVER_TYPE         g_driverType = D3D_DRIVER_TYPE_NULL;
	D3D_FEATURE_LEVEL       g_featureLevel = D3D_FEATURE_LEVEL_11_0;
	ID3D11Device* g_pd3dDevice = NULL;
	ID3D11DeviceContext* g_pImmediateContext = NULL;
	IDXGISwapChain* g_pSwapChain = NULL;
	ID3D11RenderTargetView* g_pRenderTargetView = NULL;
}

void RenderSys::Render()
{
	UINT stride = sizeof(vertex);
	UINT offset = 0;

	static float speed = 0;

	float R = 1. * mouse->wheel_pos * 0.1;
	static float h = 1.3;

	static float t = 0.0f;
	{
		//Для поворота по времени
		/* 
		static DWORD dwTimeStart = 0;
		DWORD dwTimeCur = GetTickCount();
		if (dwTimeStart == 0) dwTimeStart = dwTimeCur;
		t = (dwTimeCur - dwTimeStart) / 1000.0f;
		speed = t / 2;
		*/

		static float xPos = 0;
		static float yPos = 0;
		
		//Производим перерасчёт смешения камеры, только если нажата кнопка
		xPos = (mouse->mousePos.x - mouse->savedMPos.x) * 0.01;
		yPos = (mouse->mousePos.y - mouse->savedMPos.y) * 0.01;
		speed += xPos;
		h += yPos;
		mouse->savedMPos.x = mouse->mousePos.x;
		mouse->savedMPos.y = mouse->mousePos.y;

		DirectX::XMVECTOR Eye = DirectX::XMVectorSet(R * cos(speed), h, R * sin(speed), 0.0f);
		DirectX::XMVECTOR At = DirectX::XMVectorSet(0.0f, 0.0f, 0.0f, 0.0f);
		DirectX::XMVECTOR Up = DirectX::XMVectorSet(0.0f, 1.0f, 0.0f, 0.0f);
		g_View = DirectX::XMMatrixLookAtLH(Eye, At, Up);
	}


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
	for (std::vector<Object*>::iterator it = objects.begin(); it < objects.end(); ++it)
	{
		pObj = *it;
		offset = 0;

		g_pImmediateContext->IASetPrimitiveTopology(pObj->toplology);
		g_pImmediateContext->IASetVertexBuffers(0, 1, &pObj->pVertexBuf, &stride, &offset);
		//g_pImmediateContext->VSSetShader(pObj->pVxSh, NULL, 0);

		g_pImmediateContext->PSSetShader(pObj->pPxSh, NULL, 0);
		g_pImmediateContext->UpdateSubresource(pPS_CB_per_obj, 0, NULL, pObj->pCB, 0, 0);
		g_pImmediateContext->PSSetConstantBuffers(1, 1, &pPS_CB_per_obj);

		if (pObj->pIndexBuf)
		{
			g_pImmediateContext->IASetIndexBuffer(pObj->pIndexBuf, DXGI_FORMAT_R32_UINT, offset);
			g_pImmediateContext->DrawIndexed(pObj->indexCount, 0, 0);
			g_pImmediateContext->Draw(pObj->vecCount, 0);
		}
		else
		{
			g_pImmediateContext->Draw(pObj->vecCount, 0);
		}
	}

	g_pSwapChain->Present(1, 0);
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

	//INIT PS_CB

	D3D11_BUFFER_DESC bd1;
	ZeroMemory(&bd1, sizeof(bd1));

	bd1.Usage = D3D11_USAGE_DEFAULT;
	bd1.ByteWidth = sizeof(PS_perFrame_CB);
	bd1.BindFlags = D3D11_BIND_CONSTANT_BUFFER;
	bd1.CPUAccessFlags = 0;
	hRes = g_pd3dDevice->CreateBuffer(&bd1, NULL, &g_pPS_CB);
	if (FAILED(hRes))
		return hRes;

	Light light;
	light.Direction = XMFLOAT3(0.25f, 1, 0.50f);
	light.Ambient = XMFLOAT4(0.2f, 0.2f, 0.2f, 1.0f);
	light.Diffuse = XMFLOAT4(1.0f, 1.0f, 1.0f, 1.0f);

	g_pImmediateContext->UpdateSubresource(g_pPS_CB, 0, NULL, &light, 0, 0);
	g_pImmediateContext->PSSetConstantBuffers(0, 1, &g_pPS_CB);
	///////Per obj
	ZeroMemory(&bd1, sizeof(bd1));

	bd1.Usage = D3D11_USAGE_DEFAULT;
	bd1.ByteWidth = sizeof(PixelShaderCB);
	bd1.BindFlags = D3D11_BIND_CONSTANT_BUFFER;
	bd1.CPUAccessFlags = 0;
	hRes = g_pd3dDevice->CreateBuffer(&bd1, NULL, &pPS_CB_per_obj);
	if (FAILED(hRes))
		return hRes;

}

HRESULT RenderSys::InitDevice(HWND* hWnd)
{
	mouse = new Mouse();

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

	g_pImmediateContext->OMSetRenderTargets(1, &g_pRenderTargetView, NULL);

	// Setup the viewport
	D3D11_VIEWPORT vp;
	vp.Width = (FLOAT)width;
	vp.Height = (FLOAT)height;
	vp.MinDepth = 0.0f;
	vp.MaxDepth = 1.0f;
	vp.TopLeftX = 0;
	vp.TopLeftY = 0;
	g_pImmediateContext->RSSetViewports(1, &vp);


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

void RenderSys::drawDrappingPoints(vertex** points)
{
	//Переформирование из двумерной сетки в одномерный массив
	vertex* dummyArray = new vertex[GIRD_SIZE * GIRD_SIZE];//Массив-копия для визуализации соединения линиями

	for (size_t i = 0; i < GIRD_SIZE; ++i)
	{
		for (size_t q = 0; q < GIRD_SIZE; q++)
		{
			dummyArray[q + i * GIRD_SIZE] = points[q][i];
		}

	}
	//Вывод всей сетки на экран в виде точек
	Object* obj = new Object(this, pVxSh, pPxSh, GIRD_SIZE*GIRD_SIZE, dummyArray, D3D11_PRIMITIVE_TOPOLOGY_POINTLIST);
	objects.push_back(obj);
	

	//Формирование данных для отрисовки волокон межде узлами
	size_t arraySize = (pow((GIRD_SIZE - 1) / 2, 2) * 2 + (GIRD_SIZE - 1)) * 2;
	size_t* indexAr = new size_t[arraySize];
	memset(indexAr, 0, sizeof(size_t) * arraySize);


	vertex* lineDummy = new vertex[arraySize];//pow((GIRD_SIZE - 1) / 2 + 1, 2)];
	size_t actualIndex = 0;


	for (int a_index = (GIRD_SIZE - 1) / 2; a_index >= 0; --a_index)
	{
		for (int b_index = (GIRD_SIZE - 1) / 2; b_index >= 0; --b_index)
		{

			if (b_index - 1 >= 0)
			{
				lineDummy[actualIndex] = points[b_index][a_index];
				actualIndex++;
				lineDummy[actualIndex] = points[b_index - 1][a_index];
				actualIndex++;
			}
			if (a_index - 1 >= 0)
			{
				lineDummy[actualIndex] = points[b_index][a_index];
				actualIndex++;
				lineDummy[actualIndex] = points[b_index][a_index-1];
				actualIndex++;
			}
		}
	}
	obj = new Object(this, pVxSh, pPxSh, actualIndex, lineDummy, D3D11_PRIMITIVE_TOPOLOGY_LINELIST);
	objects.push_back(obj);

	//Анализ углов и визуализация результатов



	vertex* triangleDummy = new vertex[(GIRD_SIZE - 1) * (GIRD_SIZE - 1) * 3];
	actualIndex = 0;
	size_t i, j, p, q, s, h;

	
	float angle, blue, red, green;
	
	float coeff = 0.7;
	for (int a_index = (GIRD_SIZE - 1) / 2; a_index >= 1; --a_index)
	{
		for (int b_index = (GIRD_SIZE - 1) / 2; b_index >= 1; --b_index)
		{
			i = b_index;
			j = a_index;
			p = i - 1;
			q = j;
			s = i;
			h = j - 1;
			angle = getAngle(points, i, j, p, q, s, h);

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


			triangleDummy[actualIndex] = points[b_index][a_index];
			actualIndex++;
			triangleDummy[actualIndex] = points[b_index][a_index - 1];
			actualIndex++;
			triangleDummy[actualIndex] = points[b_index - 1][a_index];
			actualIndex++;

			triangleDummy[actualIndex] = points[b_index - 1][a_index - 1];
			actualIndex++;
			triangleDummy[actualIndex] = points[b_index - 1][a_index];
			actualIndex++;
			triangleDummy[actualIndex] = points[b_index][a_index - 1];
			actualIndex++;

			for (size_t ctr = 1; ctr <= 6; ctr++)
			{
				triangleDummy[actualIndex - ctr].Color = vec4(red, green, blue, 1.0);
			}


		}
	}
	obj = new Object(this, pVxSh, pPxSh, actualIndex, triangleDummy, D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST);
	objects.push_back(obj);
}

HRESULT RenderSys::InitObjects()
{
	HRESULT hRes = S_OK;
	//coonsLines_new();
	drawDrappingPoints(makeGird());
	//testObject();
	mp.transferToRender(this);
	return hRes;
}

void RenderSys::coonsLines_new()
{
	size_t subSplineCtM = 10;
	size_t subSplineCtN = 10;

	std::cout << sizeof(vertex) << std::endl;
	std::cout << "vec3: " << sizeof(vec3) << "; vec4: " << sizeof(vec4) << "; required: " << sizeof(vec4) + sizeof(vec3) * 2 << std::endl;

	size_t N_ctr = 10;
	size_t M_ctr = 10;

	//Step
	double teta = DirectX::XM_PIDIV2 / (M_ctr - 1);
	double delTeta = teta / subSplineCtM;
	double fi = DirectX::XM_2PI / (N_ctr - 1);
	double R = 0.98;

	//vertex* sphere = new vertex[N_ctr * subSplineCtN * M_ctr * subSplineCtM];

	int ss = N_ctr * (M_ctr - 1) * (subSplineCtM + 1) + (N_ctr - 1) * M_ctr * (subSplineCtN + 1);

	vertex* sphere = new vertex[N_ctr * (M_ctr - 1) * (subSplineCtM + 1) + (N_ctr - 1) * M_ctr * (subSplineCtN + 1)];
	vertex* line = new vertex[3];
	line[0].pos = vec3(-2, 0, 0);
	line[0].Color = vec4(1, 0, 0, 1);
	line[1].pos = vec3(0, 0, 0);
	line[1].Color = vec4(1, 0, 0, 1);
	line[2].pos = vec3(2, 0, 0);
	line[2].Color = vec4(1, 0, 0, 1);
	drawSpline(addInterpolationSpline(line, 3));

	//Vec for calculate normals
	vec3 dpdFi;
	vec3 dpdTeta;

	for (size_t t = 0; t < N_ctr; ++t)
	{
		for (size_t tau = 0; tau < M_ctr - 1; ++tau)
		{
			for (size_t tM = 0; tM <= subSplineCtM; ++tM)
			{
				vertex vtx;
				vtx.Color = vec4(0, 0, 0, 1);


				vtx.pos.x = R * cos(fi * t) * sin(teta * tau + tM * delTeta);
				vtx.pos.z = R * sin(fi * t) * sin(teta * tau + tM * delTeta);
				vtx.pos.y = R * cos(teta * tau + tM * delTeta);
				dpdFi = vec3(-R * sin(fi * t) * sin(teta * tau + tM * delTeta), 0,
					R * cos(fi * t) * sin(teta * tau + tM * delTeta));
				dpdTeta = vec3(R * cos(fi * t) * cos(teta * tau + tM * delTeta),
					-R * sin(teta * tau + tM * delTeta),
					R * sin(fi * t) * cos(teta * tau + tM * delTeta));

				XMStoreFloat3(&vtx.normal, XMVector3Cross(XMLoadFloat3(&dpdFi), XMLoadFloat3(&dpdTeta)));

				sphere[t * (subSplineCtM + 1) * (M_ctr - 1) + tau * (subSplineCtM + 1) + tM] = vtx;
			}
		}
	}


	size_t start = N_ctr * (subSplineCtM + 1) * (M_ctr - 1);//Начало записи точек для параллелей

	teta = DirectX::XM_PIDIV2 / (M_ctr - 1);
	fi = DirectX::XM_2PI / (N_ctr - 1);
	double delFi = fi / subSplineCtN;

	for (size_t tau = 0; tau < M_ctr; ++tau)
	{
		for (size_t t = 0; t < N_ctr - 1; ++t)
		{
			for (size_t tN = 0; tN <= subSplineCtN; ++tN)
			{
				vertex vtx;
				vtx.Color = vec4(0, 0, 0, 1);

				vtx.pos.x = R * cos(fi * t + tN * delFi) * sin(teta * tau);
				vtx.pos.z = R * sin(fi * t + tN * delFi) * sin(teta * tau);
				vtx.pos.y = R * cos(teta * tau);

				dpdFi = vec3(-R * sin(fi * t + tN * delFi) * sin(teta * tau), 0,
					R * cos(fi * t + tN * delFi) * sin(teta * tau));

				dpdTeta = vec3(R * cos(fi * t + tN * delFi) * cos(teta * tau),
					-R * sin(teta * tau),
					R * sin(fi * t + tN * delFi) * cos(teta * tau));

				XMStoreFloat3(&vtx.normal, XMVector3Cross(XMLoadFloat3(&dpdFi), XMLoadFloat3(&dpdTeta)));



				sphere[start + tau * (subSplineCtN + 1) * (N_ctr - 1) + t * (subSplineCtN + 1) + tN] = vtx;
			}
		}
	}

	splineInfo* spI = new splineInfo[N_ctr * (M_ctr - 1) + (N_ctr - 1) * M_ctr];// Набор сплайнов без изменений(для рисования сферы)
	splineInfo* patchSpI = new splineInfo[(N_ctr - 1) * M_ctr];// Набор сплайнов для изменений(для рисования лоскутов) только для параллелейы

	char bu1f[100];
	size_t bias = N_ctr * (M_ctr - 1);

	for (size_t t = 0; t < N_ctr; ++t)
	{
		for (size_t j = 0; j < M_ctr - 1; ++j)
		{
			spI[t * (M_ctr - 1) + j] = addInterpolationSpline((sphere + t * (subSplineCtM + 1) * (M_ctr - 1) + (subSplineCtM + 1) * j), subSplineCtM + 1);
			drawSpline(spI[t * (M_ctr - 1) + j]);

		}
	}

	for (size_t j = 0; j < M_ctr; ++j)
	{
		for (size_t t = 0; t < N_ctr - 1; ++t)
		{
			//++splineNct;
			spI[bias + j * (N_ctr - 1) + t] = addInterpolationSpline((sphere + start + j * (subSplineCtN + 1) * (N_ctr - 1) + (subSplineCtN + 1) * t), subSplineCtN + 1);

			drawSpline(spI[bias + j * (N_ctr - 1) + t]);
		}
	}
	//drawBigPoints(sphere, N_ctr * (M_ctr - 1) * (subSplineCtM + 1) + (N_ctr - 1) * M_ctr * (subSplineCtN + 1), 0.01);

	for (size_t i = 0; i < (N_ctr - 1) * M_ctr; ++i)
	{
		patchSpI[i] = spI[bias + i];
	}


	////coons
	//
	size_t spIBias = N_ctr * (M_ctr - 1);
	size_t omegaSize = 5;
	size_t ksiSize = 5;
	vertex** surface = new vertex * [(N_ctr - 1) * (M_ctr - 1)];



	for (size_t i = 0; i < (N_ctr - 1) * (M_ctr - 1); ++i)
	{
		surface[i] = new vertex[omegaSize * ksiSize];
	}

	double tStep = 1.0 / (omegaSize - 1);
	double tauStep = 1.0 / (ksiSize - 1);

	double omega = 0;
	double ksi = 0;


	//Makeing changes in spline length
	/*
	double eps = 0.01;
	bool isGrow = true;


	for (size_t s = 0; s < N_ctr - 1; ++s)
	{
		for (size_t j = 1; j < M_ctr - 2; ++j)
		{

			splineInfo* cc1 = &spI[spIBias + (M_ctr - j) * (N_ctr - 1) + s];//bottom
			splineInfo* cc2 = &patchSpI[(M_ctr - j - 1) * (N_ctr - 1) + s];//up

			double botomLen = getSplineLen(0, 1, *CurveDerivateAlg1, *cc1);
			double upLen = getSplineLen(0, 1, *CurveDerivateAlg1, *cc2);
			printf("up len = %lf\nbottom len = %lf\n", upLen, botomLen);

			vertex* vx = &cc2->controlPoints[subSplineCtM / 2];
			vx->normal = vec3(2 * vx->pos.x, 2 * vx->pos.y, 2 * vx->pos.z);

			size_t counter = 0;
			while (fabs(botomLen - upLen) > eps)
			{
				counter++;
				if (upLen < botomLen)
				{
					vx->pos = vx->pos + vx->normal * 0.0001;
				}
				else
				{
					vx->pos = vx->pos + vx->normal * (-0.0001);
				}
				botomLen = getSplineLen(0, 1, *CurveDerivateAlg1, *cc1);
				upLen = getSplineLen(0, 1, *CurveDerivateAlg1, *cc2);
			}
			printf("count is %d: \n", counter);
		}
	}
	for (size_t i = 0; i < M_ctr; ++i)
	{
		splineInfo c1 = spI[bias + i * (N_ctr - 1)];

		vec3 i1 = integrate(0, 1, *CurveDerivateAlg1, c1);
		double dii1 = sqrt(pow(i1.x, 2) + pow(i1.y, 2) + pow(i1.z, 2));
		//std::cout << dii1 << std::endl;
		printf("%d : %lf\n", i, dii1);
	}*/

	//Make patches
	for (size_t i = 0; i < (M_ctr - 1); ++i)
		for (size_t j = 0; j < (N_ctr - 1); ++j)
		{
			splineInfo c1 = spI[i * (M_ctr - 1) + j];
			splineInfo c2 = spI[(i + 1) * (M_ctr - 1) + j];
			splineInfo d1 = patchSpI[j * (N_ctr - 1) + i];
			splineInfo d2 = spI[spIBias + (j + 1) * (N_ctr - 1) + i];



			double t = 0;
			double tau = 0;

			for (size_t o = 0; o < omegaSize; ++o)
			{
				tau = 0;
				omega = t;

				for (size_t k = 0; k < ksiSize; ++k)
				{
					ksi = tau;

					surface[i * (N_ctr - 1) + j][o * ksiSize + k] = (1 - ksi) * CurvePoint(c2, 3, t) + ksi * CurvePoint(c1, 3, t);
					surface[i * (N_ctr - 1) + j][o * ksiSize + k] = surface[i * (N_ctr - 1) + j][o * ksiSize + k] + (1 - omega) * CurvePoint(d2, 3, tau) + omega * CurvePoint(d1, 3, tau);
					surface[i * (N_ctr - 1) + j][o * ksiSize + k] = surface[i * (N_ctr - 1) + j][o * ksiSize + k] - (CurvePoint(c2, 3, 0) * (1 - omega) * (1 - ksi) + CurvePoint(c2, 3, 1) * omega * (1 - ksi) + CurvePoint(c1, 3, 0) * (1 - omega) * ksi + CurvePoint(c1, 3, 1) * omega * ksi);
					surface[i * (N_ctr - 1) + j][o * ksiSize + k].Color = vec4(0, 1, 0, 1);

					tau += tauStep;
					if (tau > 1) tau = 1;
				}
				t += tStep;
				if (t > 1) t = 1;
			}
		}


	vertex* pPt;
	for (size_t i = 0; i < (M_ctr - 1) * (N_ctr - 1); ++i)
	{
		for (size_t j = 0; j < omegaSize * ksiSize; ++j)
		{
			pPt = &surface[i][j];
			pPt->normal = vec3(2 * pPt->pos.x, 2 * pPt->pos.y, 2 * pPt->pos.z);
		}
	}

	vertex test[3];

	//size_t j = 0;
	for (size_t s = 0; s < (M_ctr - 1); ++s)
		for (size_t j = 0; j < (N_ctr - 1); ++j)
		{
			for (size_t o = 0; o < omegaSize - 1; ++o)
			{
				for (size_t i = 0; i < omegaSize - 1; ++i)
				{
					test[1] = surface[s * (N_ctr - 1) + j][o * ksiSize + i];
					test[0] = surface[s * (N_ctr - 1) + j][o * ksiSize + i + 1];
					test[2] = surface[s * (N_ctr - 1) + j][o * ksiSize + ksiSize + i];
					drawTriangle(test);

					test[1] = surface[s * (N_ctr - 1) + j][o * ksiSize + ksiSize + i + 1]; //Выход за границу!!
					test[2] = surface[s * (N_ctr - 1) + j][o * ksiSize + i + 1];
					test[0] = surface[s * (N_ctr - 1) + j][o * ksiSize + ksiSize + i];
					drawTriangle(test);
				}
			}
		}


	{
		splineInfo c1 = spI[0 * (M_ctr - 1) + 0];
		splineInfo c2 = spI[(0 + 1) * (M_ctr - 1) + 0];
		splineInfo d1 = patchSpI[0 * (N_ctr - 1) + 0];
		splineInfo d2 = spI[spIBias + (0 + 1) * (N_ctr - 1) + 0];

		splineInfo ar_sp[4];
		ar_sp[0] = c1;
		ar_sp[1] = c2;
		ar_sp[2] = d1;
		ar_sp[3] = d2;

		vertex* ar[4];
		ar[0] = new vertex[100];
		ar[1] = new vertex[100];
		ar[2] = new vertex[100];
		ar[3] = new vertex[100];
		for (int i = 0; i < 4; ++i)
		{
			for (int j = 0; j < 100; ++j)
			{
				ar[i][j] = CurvePoint(ar_sp[i], 3, (double)(j) / 99.0);
				ar[i][j].Color = vec4(1, 0, 0, 1);
			}



			Object* obj = new Object(this, pVxSh, pPxSh, 100, ar[i], D3D11_PRIMITIVE_TOPOLOGY_POINTLIST);
			objects.push_back(obj);
		}
		std::cout << "Hello: " << getArea(&c1, &c2, &d1, &d2) << "\n";
	}


	for (size_t i = 0; i < (N_ctr - 1) * (M_ctr - 1); ++i)
	{
		delete[] surface[i];
	}

	delete[] surface;
}

HRESULT RenderSys::drawSpline(splineInfo _info)
{
	HRESULT hRes = S_OK;
	UINT offset = 0;
	UINT stride = sizeof(vertex);

	vertex vertices[100];
	splineInfo info = _info;
	vertex* points = info.controlPoints;

	double* knots = info.knotVector;

	vertex* pt = makeBSpline(100, 3, info);

	for (int i = 0; i < 100; ++i)
	{
		vertices[i] = pt[i];
	}

	Object* obj = new Object(this, pVxSh, pPxSh, 100, vertices, D3D11_PRIMITIVE_TOPOLOGY_LINESTRIP);
	objects.push_back(obj);
	//hRes = drawBigPoints(points, info.cpCount, 0.01);
/*
		delete[] pt;
		delete[] knots;
		delete[] points;*/

	return hRes;
}

HRESULT RenderSys::drawTriangle(vertex* _pt)
{
	HRESULT hRes = S_OK;
	mp.addNew(_pt, 3);
	return hRes;

}

Mouse* RenderSys::getMouse()
{
	return mouse;
}


Mouse::Mouse()
{
	isLeftKeyPressed = false;
	isRightKeyPressed = false;
	wheel_pos = 1;
}

void Mouse::updWheelPos(int newPos)
{
	wheel_pos = newPos;
}

void Mouse::updLK(bool isPressed)
{
	isLeftKeyPressed = isPressed;
	savedMPos = mousePos;
}

void Mouse::updRK(bool isPressed)
{
	isRightKeyPressed = isPressed;
}

void Mouse::updMousePos(POINT mPos)
{
	if (!isLeftKeyPressed)
	{
		savedMPos = mPos;
	}
	mousePos = mPos;
}
