#pragma once
#include <d3d11.h>
#include <DirectXMath.h>

#include "Render/Object.h"
#include "Render/Render.h"
#include "../Math/Drapping.h"

using namespace DirectX;

class RenderSys;

void lighting_test(RenderSys* _rs, ID3D11Device* _pDevice, ID3D11VertexShader* _pVxSh, ID3D11PixelShader* _pPxSh);

void generateSphere(RenderSys* _rs, ID3D11Device* _pDevice, ID3D11VertexShader* _pVxSh, ID3D11PixelShader* _pPxSh);

//void generateSphereByBSpline(RenderSys* _rs, ID3D11Device* _pDevice, ID3D11VertexShader* _pVxSh, ID3D11PixelShader* _pPxSh);
//
void generateSurfaceByBSpline(RenderSys* _rs, ID3D11Device* _pDevice, ID3D11VertexShader* _pVxSh, ID3D11PixelShader* _pPxSh);

void testPlane(RenderSys* _rs, ID3D11Device* _pDevice, ID3D11VertexShader* _pVxSh, ID3D11PixelShader* _pPxSh);

void cornerDrapping(RenderSys* _rs, ID3D11Device* _pDevice, ID3D11VertexShader* _pVxSh, ID3D11PixelShader* _pPxSh);