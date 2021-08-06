
struct VS_OUTPUT
{
	float4 Pos : SV_POSITION;
	float4 Color : COLOR;
	float3 Normal : NORMAL;
	float2 Tex: TEXCOORD;
};


Texture2D gDiffuseMap;

SamplerState samAnisotropic
{
	Filter = ANISOTROPIC;
	MaxAnisotropy = 4;
	AddressU = WRAP;
	AddressV = WRAP;
};

float4 main(VS_OUTPUT input) : SV_TARGET
{
	return gDiffuseMap.Sample(samAnisotropic, input.Tex);
}