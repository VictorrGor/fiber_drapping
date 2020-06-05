cbuffer vertexCB //: register(b0)
{
	float4x4 wvp;
	float4x4 world;
}
//--------------------------------------------------------------------------------------
struct VS_OUTPUT
{
	float4 Pos : SV_POSITION;
	float4 WorldPos : POSITION;
	float4 Color : COLOR;
	float3 Normal : NORMAL;
};
struct VS_INPUT
{
	float3 Pos : POSITION;
	float4 Color : COLOR;
	float3 Normal : NORMAL;
};

//--------------------------------------------------------------------------------------
// Vertex Shader
//--------------------------------------------------------------------------------------
VS_OUTPUT main(VS_INPUT inp)
{
	VS_OUTPUT output = (VS_OUTPUT)0;

	output.Pos = mul(float4(inp.Pos, 1.0f), wvp);
	output.WorldPos = mul(float4(inp.Pos, 1.0f), world);
	output.Color = inp.Color;
	output.Normal = mul(inp.Normal, world);// obr transp 8ая лекция
	return output;
}
