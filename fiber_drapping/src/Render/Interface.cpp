// This is a personal academic project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: http://www.viva64.com
#include "../Render/Interface.h"


HRESULT Interface::initCharMap(const char* filename)///@todo needs testing
{
	FILE* f;
	f = fopen(filename, "rb");
	if (!f) return E_FAIL;


	UINT charactersNum = 0;
	fseek(f, -((int)sizeof(UINT)), SEEK_END);
	fread(&charactersNum, sizeof(UINT), 1, f);
	fseek(f, 0, SEEK_SET);

	TextTextureInfo* buf = textureInfoBuf;
	buf = new TextTextureInfo[charactersNum];
	fread(&buf[0], sizeof(TextTextureInfo), charactersNum, f);
	auto a = buf[0];
	for (UINT i = 0; i < charactersNum; ++i)
	{
		chars.insert(std::pair<char, TextTextureInfo*>(buf[i].ch, &buf[i]));
		std::cout << buf[i].ch << ": LeftU: " << chars[buf[i].ch]->leftU << "; RightU: " << chars[buf[i].ch]->rightU << "; Width: " << chars[buf[i].ch]->pixelWidth << "\n";
	}
	fclose(f);
	return S_OK;
};

Interface::Interface(ID3D11Device* _pDevice, UINT _screenWidth, UINT _screenHeight) : screenHeight(_screenHeight), screenWidth(_screenWidth)
{
	HRESULT hRes = S_OK;

	///////////////////
	//Compiling shaders
	///////////////////

	ID3DBlob* pBlobVertex = NULL;
	ID3DBlob* pPSerr = NULL;
	UINT compile_flag = D3DCOMPILE_DEBUG || D3DCOMPILE_SKIP_OPTIMIZATION;


	hRes = D3DCompileFromFile(L"shaders/InterfaceVertexShader.hlsl", NULL, D3D_COMPILE_STANDARD_FILE_INCLUDE, "main", "vs_4_0", compile_flag, NULL, &pBlobVertex, &pPSerr);
	if (FAILED(hRes))
	{
		if (pPSerr) pPSerr->Release();
		MessageBox(NULL,
			L"The FX file cannot be compiled.  Please run this executable from the directory that contains the FX file.", L"Error", MB_OK);
		throw "Cann't init Interface!";
	}

	hRes = _pDevice->CreateVertexShader(pBlobVertex->GetBufferPointer(), pBlobVertex->GetBufferSize(), NULL, &pInterfaceVxSh);


	// Define the input layout
	D3D11_INPUT_ELEMENT_DESC layout[] =
	{
		{ "POSITION", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, D3D11_APPEND_ALIGNED_ELEMENT , D3D11_INPUT_PER_VERTEX_DATA, 0 },
		{ "COLOR", 0, DXGI_FORMAT_R32G32B32A32_FLOAT, 0, D3D11_APPEND_ALIGNED_ELEMENT , D3D11_INPUT_PER_VERTEX_DATA, 0 },
		{ "NORMAL", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, D3D11_APPEND_ALIGNED_ELEMENT , D3D11_INPUT_PER_VERTEX_DATA, 0 },
		{ "TEXCOORD", 0, DXGI_FORMAT_R32G32_FLOAT, 0, D3D11_APPEND_ALIGNED_ELEMENT , D3D11_INPUT_PER_VERTEX_DATA, 0 }

	};
	UINT numElements = ARRAYSIZE(layout);

	// Create the input layout
	hRes = _pDevice->CreateInputLayout(layout, numElements, pBlobVertex->GetBufferPointer(),
		pBlobVertex->GetBufferSize(), &pInterfaceVtxLayout);
	pBlobVertex->Release();
	if (FAILED(hRes))
		throw "Cann't init Interface!";

	// Compile the pixel shader
	ID3DBlob* pPSBlob = NULL;
	hRes = D3DCompileFromFile(L"shaders/InterfacePixelShader.hlsl", NULL, D3D_COMPILE_STANDARD_FILE_INCLUDE, "main", "ps_4_0", NULL, NULL, &pPSBlob, &pPSerr);
	if (FAILED(hRes))
	{
		MessageBox(NULL,
			L"The FX file cannot be compiled.  Please run this executable from the directory that contains the FX file.", L"Error", MB_OK);
		throw "Cann't init Interface!";
	}

	// Create the pixel shader
	hRes = _pDevice->CreatePixelShader(pPSBlob->GetBufferPointer(), pPSBlob->GetBufferSize(), NULL, &pInterfacePxSh);
	pPSBlob->Release();
	if (FAILED(hRes))
		throw "Cann't init Interface!";

	if (pPSerr) pPSerr->Release();

	//////////////////////////
	//Init Texture Resources//
	//////////////////////////
	hRes = CreateDDSTextureFromFile(_pDevice, L"resources/texture/font.dds", (ID3D11Resource**)&pFont, &pFontRV);

	if (FAILED(hRes))
	{
		MessageBox(NULL,
			L"Some problems with Creating texture.", L"Error", MB_OK);
		throw "Failed to create texture! Cann't init Interface!";
	}

	initCharMap("resources/texture/ASCI_desc.binary");
}

Interface::~Interface()
{
	if (pInterfaceVxSh)
	{
		pInterfaceVxSh->Release();
	}
	if (pInterfacePxSh)
	{
		pInterfacePxSh->Release();
	}
	if (pInterfaceVtxLayout)
	{
		pInterfaceVtxLayout->Release();
	}
	if (pFont)
	{
		pFont->Release();
	}
	if (pFontRV)
	{
		pFontRV->Release();
	}
	for (std::vector<TexturedObject*>::iterator it = objects.begin(); it < objects.end(); ++it)
	{
		if ((*it))
		{
			delete (*it);
			(*it) = nullptr;
		}
	}
	objects.clear();
	for (std::vector<ActiveElement*>::iterator it = activeObjs.begin(); it < activeObjs.end(); ++it)
	{
		if ((*it))
		{
			delete (*it);
			(*it) = nullptr;
		}
	}
	activeObjs.clear();
}

TexturedObject* Interface::makeWord(const char* _text, float x_pos, float y_pos, float _fontSize, ID3D11Device* _pDevice) {
	UINT size = strlen(_text);
	TexturedVertex* vec = new TexturedVertex[size * 4];
	float x_act = x_pos;
	TextTextureInfo* tti = nullptr;
	UINT counter = 0;
	float x_step = 0;
	float y_step = _fontSize / screenWidth * 20;

	UINT* indices = new UINT[size * 6];
	UINT indices_ctr = 0;
	for (int i = 0; i < size; ++i)
	{
		tti = chars[_text[i]];

		x_step = (float)tti->pixelWidth / screenWidth * _fontSize;

		vec[counter].vx.pos = vec3(x_act, y_pos, 0);
		vec[counter].tcoord.x = tti->leftU;
		vec[counter].tcoord.y = 0;

		vec[counter + 1].vx.pos = vec3(x_act + x_step, y_pos, 0);
		vec[counter + 1].tcoord.x = tti->rightU;
		vec[counter + 1].tcoord.y = 0;

		vec[counter + 2].vx.pos = vec3(x_act, y_pos - y_step, 0);
		vec[counter + 2].tcoord.x = tti->leftU;
		vec[counter + 2].tcoord.y = 1;

		vec[counter + 3].vx.pos = vec3(x_act + x_step, y_pos - y_step, 0);
		vec[counter + 3].tcoord.x = tti->rightU;
		vec[counter + 3].tcoord.y = 1;

		indices[indices_ctr] = counter;
		indices[indices_ctr + 1] = counter + 1;
		indices[indices_ctr + 2] = counter + 2;
		indices[indices_ctr + 3] = counter + 1;
		indices[indices_ctr + 4] = counter + 3;
		indices[indices_ctr + 5] = counter + 2;

		indices_ctr += 6;
		x_act += x_step;
		counter += 4;
	}
	TexturedObject* res = new TexturedObject(_pDevice, this->pInterfaceVxSh, this->pInterfacePxSh, size * 4, (vec), D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST, indices_ctr, indices);
	delete[] vec;
	delete[] indices;
	return res;
}

void Interface::pushTexturedObject(TexturedObject* _obj)
{
	objects.push_back(_obj);
}

void Interface::pushTexturedObject(const char* _name, TexturedObject* _obj)
{
	mObjects.insert(std::pair<const char*, TexturedObject*>(_name, _obj));
}

bool Interface::removeTexturedObject(const char* _name)
{
	auto res = mObjects.find(_name);
	if (res == mObjects.end()) return false;
	mObjects.erase(res);
	return true;
}

void Interface::pushActiveElement(ActiveElement* _obj)
{
	activeObjs.push_back(_obj);
}

void Interface::deleteObject(TexturedObject* _obj)
{
	for (std::vector<TexturedObject*>::iterator it = objects.begin(); it != objects.end(); ++it)
	{
		if (_obj == (*it))
		{
			delete (*it);
			objects.erase(it);
			break;
		}
	}
}

void Interface::Render(ID3D11DeviceContext* _pDeviceContext)
{
	_pDeviceContext->IASetInputLayout(pInterfaceVtxLayout);
	_pDeviceContext->VSSetShader(pInterfaceVxSh, NULL, 0);
	_pDeviceContext->PSSetShader(pInterfacePxSh, NULL, 0);
	_pDeviceContext->PSSetShaderResources(0, 1, &pFontRV);

	UINT stride = sizeof(TexturedVertex);
	UINT offset = 0;
	for (std::vector<TexturedObject*>::iterator it = objects.begin(); it < objects.end(); ++it)
	{
		(*it)->Render(_pDeviceContext);
	}

	for (std::vector<ActiveElement*>::iterator it = activeObjs.begin(); it < activeObjs.end(); ++it)
	{
		(*it)->Render(_pDeviceContext);
	}
}

void Interface::Resize(UINT _screenWidth, UINT _screenHeight)
{
	screenWidth  = _screenWidth;
	screenHeight = _screenHeight;
}

void Interface::SendClick(POINT _pt)
{
	for (std::vector<ActiveElement*>::iterator it = activeObjs.begin(); it != activeObjs.end(); ++it)
	{
		if ((*it)->isClicked(_pt))
		{
			break;
		}
	}
}

void Interface::makeTextObject(const char* _text, float x_pos, float y_pos, float _fontSize, ID3D11Device* _pDevice)
{
	objects.push_back(makeWord(_text, x_pos, y_pos,  _fontSize, _pDevice));
}

ActiveElement::ActiveElement()
{
	obj = nullptr;
	payload = nullptr;
}

ActiveElement::~ActiveElement()
{
	if (obj)
	{
		delete obj;
		obj = nullptr;
	}
}


bool ActiveElement::isClicked(POINT _pt)
{
	if ((_pt.x >= x) && (_pt.x <= x + boundWidth) && (_pt.y >= y) && (_pt.y <= y + boundHeight))
	{
		if (payload) payload();
		return true;
	}
	return false;
}

void ActiveElement::setPayload(void(*_payload)())
{
	payload = _payload;
}

void ActiveElement::setObj(TexturedObject* _obj)
{
	if (obj) throw "Button::setObj: Memory leak! Rewriting pointer";
	obj = _obj;
}

void ActiveElement::setCoordinates(float _screenX, float _screenY, UINT _boundWidth, UINT _boundHeight, UINT _screenWidth, UINT _screenHeight)
{
	float absX = (_screenX + 1.) / 2.;
	float absY = (-_screenY + 1.) / 2.;
	x = _screenWidth  * absX;
	y = _screenHeight * absY;
	boundWidth = _boundWidth;
	boundHeight = _boundHeight;
}

void ActiveElement::Render(ID3D11DeviceContext* _pDeviceContext)
{
	obj->Render(_pDeviceContext);
}
