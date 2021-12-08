#pragma once
#include <d3d11.h>
#include <DirectXMath.h>
#include "DDSTextureLoader.h"
#include <d3dcompiler.h>
#include <map>
#include <vector>
#include <iostream>
#include <fstream>

#include "RenderStructures.h"
#include "Object.h"

using namespace DirectX;

class ActiveElement
{
protected:
	void (*payload)();
	TexturedObject* obj;
	UINT boundWidth;
	UINT boundHeight;
	float x;
	float y;
public:
	ActiveElement();
	virtual ~ActiveElement();
	bool isClicked(POINT _pt);
	void setPayload(void (*_payload)());
	void setObj(TexturedObject* _obj);
	void setCoordinates(float _screenX, float _screenY, UINT _boundWidth, UINT _boundHeight, UINT _screenWidth, UINT _screenHeight);
	void Render(ID3D11DeviceContext* _pDeviceContext);
};

class Button: public ActiveElement
{
};

class Interface
{
	std::map<char, TextTextureInfo*> chars;
	std::map<const char*, TexturedObject*> mObjects;
	TextTextureInfo* textureInfoBuf;
	//File format:
	//Binary file. First element has UINT type and stores count of characters.
	//Further written character's structers.
	HRESULT initCharMap(const char* filename);
	std::vector<TexturedObject*> objects;
	std::vector<ActiveElement*> activeObjs;

	ID3D11VertexShader* pInterfaceVxSh;
	ID3D11PixelShader* pInterfacePxSh;
	ID3D11InputLayout* pInterfaceVtxLayout;

	ID3D11Texture2D* pFont;
	ID3D11ShaderResourceView* pFontRV;

	UINT screenWidth;
	UINT screenHeight;
public:
	Interface(ID3D11Device* _pDevice, UINT _screenWidth, UINT _screenHeight);
	~Interface();
	void makeTextObject(const char* _text, float x_pos, float y_pos, float _fontSize, ID3D11Device* _pDevice);

	//@todo Height and width are need to be used
	TexturedObject* makeWord(const char* _text, float x_pos, float y_pos, float _fontSize, ID3D11Device* _pDevice);
	void pushTexturedObject(TexturedObject* _obj);
	
	void pushTexturedObject(const char* _name, TexturedObject* _obj);
	bool removeTexturedObject(const char* _name); // if deletion was succesful return true. False otherwise

	void pushActiveElement(ActiveElement* _obj);
	void deleteObject(TexturedObject* _obj);

	void Render(ID3D11DeviceContext* _pDeviceContext);
	///@todo resize all objects too
	void Resize(UINT _screenWidth, UINT _screenHeight);

	void SendClick(POINT _pt);
};
